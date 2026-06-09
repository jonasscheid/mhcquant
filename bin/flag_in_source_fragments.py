#!/usr/bin/env python
"""Flag in-source fragments (ISF) in rescored MHCquant identifications.

In-source dissociation cleaves intact peptides in the ESI source before MS2 isolation. The
resulting b/y-type truncations co-elute with their intact parent but a sequence-based RT
predictor (DeepLC) predicts them to elute elsewhere. This script reads the FDR-filtered idXML
(the confident parents) and the 100% FDR rescored idXML (the ISF candidate source), tags every
output PeptideHit with an ``is_isf`` UserParam, and appends the flagged ISF candidates to the
filtered set so they reach quantification.

Matching is done per origin run (``id_merge_index``); RT is only comparable within a run.

Usage:
    flag_in_source_fragments.py --filtered <fdr_filtered.idXML> \
        --rescored <pout.idXML> --out <out.idXML>
"""

import argparse
import logging
import os
import sys

import numpy as np
import pyopenms as oms

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
log = logging.getLogger("flag_isf")

# Hardcoded detection defaults (not user-facing; see design spec).
MIN_OVERLAP_LEN = 6        # minimum fragment length to consider
DELTA_OBS_TOL = 5.0        # s, co-elution tolerance (observed RT)
PRED_TOL_PERCENTILE = 95.0  # percentile of parent DeepLC residuals -> Δpred tolerance
PRED_TOL_FLOOR = 60.0      # s, fallback Δpred tolerance when too few parents


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Flag in-source fragments in rescored idXML.")
    p.add_argument("--filtered", required=True, help="FDR-filtered idXML (confident parents)")
    p.add_argument("--rescored", required=True, help="100%% FDR rescored idXML (ISF candidates)")
    p.add_argument("--out", required=True,
                   help="output filtered idXML (whitelist + ID export): filtered hits + flagged ISF")
    p.add_argument("--out-rescored", default=None,
                   help="output rescored idXML with is_isf annotated on every PSM (fed to "
                        "quantification so the flag reaches the final table)")
    p.add_argument("--min-len", type=int, default=MIN_OVERLAP_LEN)
    p.add_argument("--delta-obs", type=float, default=DELTA_OBS_TOL)
    p.add_argument("--delta-pred", type=float, default=None,
                   help="Δpredicted-RT tolerance in s (default: data-driven per run)")
    p.add_argument("--pred-percentile", type=float, default=PRED_TOL_PERCENTILE)
    return p.parse_args(argv)


def load_idxml(path):
    prot_ids, pep_ids = [], []
    oms.IdXMLFile().load(path, prot_ids, pep_ids)
    return prot_ids, pep_ids


def spectrum_ref(pid):
    try:
        ref = pid.getSpectrumReference()
        if ref:
            return ref
    except Exception:
        pass
    if pid.metaValueExists("spectrum_reference"):
        return pid.getMetaValue("spectrum_reference")
    return None


def spectra_data_map(prot_ids):
    """Map each ProteinIdentification identifier -> list of run names (spectra_data basenames).

    The merged idXML stores the ordered list of source runs as the `spectra_data` MetaValue on
    the ProteinIdentification; each PeptideIdentification points into it via `id_merge_index`.
    """
    out = {}
    for pr in prot_ids:
        runs = []
        if pr.metaValueExists("spectra_data"):
            for s in pr.getMetaValue("spectra_data"):
                s = s.decode() if isinstance(s, (bytes, bytearray)) else str(s)
                runs.append(os.path.basename(s))
        out[pr.getIdentifier()] = runs
    return out


def resolve_run(pid, sd_map):
    """Resolve the origin run NAME of a PeptideIdentification via spectra_data[id_merge_index].

    Using the run name (not the raw index) makes parent/candidate matching robust even if the
    two idXMLs order their spectra_data differently.
    """
    runs = sd_map.get(pid.getIdentifier(), [])
    if pid.metaValueExists("id_merge_index"):
        idx = int(pid.getMetaValue("id_merge_index"))
        if 0 <= idx < len(runs):
            return runs[idx]
        return "idx%d" % idx
    return runs[0] if runs else "run0"


def _meta_float(hit, key):
    if hit.metaValueExists(key):
        try:
            return float(hit.getMetaValue(key))
        except (TypeError, ValueError):
            return None
    return None


def best_hit_record(pid, run_name):
    """Return a dict describing the rank-1 hit of a PeptideIdentification, or None.

    Matching is strictly per origin run: ``run_name`` is the scope for RT comparison and part
    of the unique identity key.
    """
    hits = pid.getHits()
    if not hits:
        return None
    h = hits[0]
    seq = h.getSequence().toUnmodifiedString()
    obs = _meta_float(h, "observed_retention_time_best")
    if obs is None:
        obs = pid.getRT()
    sref = spectrum_ref(pid)
    return {
        "run": run_name,                    # origin run: scope for RT comparison
        "key": (run_name, sref),            # identity: origin run + spectrum (unique)
        "sequence": seq,
        "rt_obs": obs,
        "rt_pred": _meta_float(h, "predicted_retention_time_best"),
        "rt_diff_best": _meta_float(h, "rt_diff_best"),
        "target_decoy": h.getMetaValue("target_decoy") if h.metaValueExists("target_decoy") else "target",
    }


def terminal_substrings(seq, min_len):
    """Yield every contiguous prefix/suffix of seq with length in [min_len, len(seq)-1]."""
    n = len(seq)
    for k in range(min_len, n):
        yield seq[:k]
        yield seq[n - k:]


def is_prefix_or_suffix(parent_seq, frag_seq):
    return (len(frag_seq) < len(parent_seq)
            and (parent_seq.startswith(frag_seq) or parent_seq.endswith(frag_seq)))


def derive_pred_tol(parents_by_run, percentile):
    resid = [p["rt_diff_best"] for run in parents_by_run.values()
             for p in run if p["rt_diff_best"] is not None]
    if len(resid) >= 5:
        return float(np.percentile(resid, percentile))
    return PRED_TOL_FLOOR


def build_parent_index(parents, min_len):
    """Map each prefix/suffix of every parent (in one run) to the parent records."""
    index = {}
    for p in parents:
        seq = p["sequence"]
        if seq is None or len(seq) <= min_len or p["rt_obs"] is None or p["rt_pred"] is None:
            continue
        for sub in set(terminal_substrings(seq, min_len)):
            index.setdefault(sub, []).append(p)
    return index


def flag_candidate(cand, parent_index, delta_obs, delta_pred):
    """Return the best matching parent dict if cand is an ISF, else None."""
    if (cand["sequence"] is None or len(cand["sequence"]) < MIN_OVERLAP_LEN
            or cand["rt_obs"] is None or cand["rt_pred"] is None):
        return None
    candidates = parent_index.get(cand["sequence"])
    if not candidates:
        return None
    best = None
    for p in candidates:
        if p["key"] == cand["key"]:
            continue  # never match a PSM to itself
        if not is_prefix_or_suffix(p["sequence"], cand["sequence"]):
            continue
        d_obs = cand["rt_obs"] - p["rt_obs"]
        d_pred = cand["rt_pred"] - p["rt_pred"]
        if abs(d_obs) <= delta_obs and abs(d_pred) >= delta_pred:
            if best is None or abs(d_obs) < best[0]:
                best = (abs(d_obs), p)
    return best[1] if best else None


def main(argv=None):
    args = parse_args(argv)

    log.info("Loading filtered idXML (parents): %s", args.filtered)
    filt_prot, filt_pids = load_idxml(args.filtered)
    log.info("Loading rescored idXML (candidates): %s", args.rescored)
    resc_prot, resc_pids = load_idxml(args.rescored)

    filt_sd = spectra_data_map(filt_prot)
    resc_sd = spectra_data_map(resc_prot)

    # Parents (confident) grouped strictly per origin run; identity tracked by unique key.
    parents_by_run = {}
    filtered_keys = set()
    for pid in filt_pids:
        rec = best_hit_record(pid, resolve_run(pid, filt_sd))
        if rec is None:
            continue
        filtered_keys.add(rec["key"])
        if rec["target_decoy"] != "decoy":
            parents_by_run.setdefault(rec["run"], []).append(rec)

    delta_pred = args.delta_pred if args.delta_pred is not None \
        else derive_pred_tol(parents_by_run, args.pred_percentile)
    log.info("Parents: %d in %d run(s); delta_obs<=%.1fs, delta_pred>=%.1fs",
             len(filtered_keys), len(parents_by_run), args.delta_obs, delta_pred)

    parent_index_by_run = {run: build_parent_index(ps, args.min_len)
                           for run, ps in parents_by_run.items()}

    # Walk rescored identifications once. Annotate is_isf on EVERY rescored PSM (so the flag
    # survives into quantification, which rips/quantifies from the rescored runs), and collect
    # the augmented filtered subset (filtered hits + flagged ISF candidates) used as the
    # whitelist / identification export.
    rescored_out = []
    filtered_out = []
    n_isf = n_isf_subthreshold = 0
    for pid in resc_pids:
        rec = best_hit_record(pid, resolve_run(pid, resc_sd))
        if rec is None:
            continue
        in_filtered = rec["key"] in filtered_keys
        isf_parent = None
        if rec["target_decoy"] != "decoy":
            isf_parent = flag_candidate(rec, parent_index_by_run.get(rec["run"], {}),
                                        args.delta_obs, delta_pred)
        is_isf = isf_parent is not None
        # annotate the rank-1 hit and persist
        hits = pid.getHits()
        hits[0].setMetaValue("is_isf", "true" if is_isf else "false")
        if is_isf:
            hits[0].setMetaValue("isf_parent_sequence", isf_parent["sequence"])
        pid.setHits(hits)
        rescored_out.append(pid)
        if in_filtered or is_isf:
            filtered_out.append(pid)
        if is_isf:
            n_isf += 1
            if not in_filtered:
                n_isf_subthreshold += 1

    log.info("is_isf=true: %d (of which %d sub-threshold appended to filtered); "
             "filtered out=%d, rescored out=%d",
             n_isf, n_isf_subthreshold, len(filtered_out), len(rescored_out))

    # Use the rescored ProteinIdentification (superset) so all protein_refs resolve.
    oms.IdXMLFile().store(args.out, resc_prot, filtered_out)
    log.info("Wrote filtered (whitelist + ID export): %s", args.out)
    if args.out_rescored:
        oms.IdXMLFile().store(args.out_rescored, resc_prot, rescored_out)
        log.info("Wrote rescored (is_isf-annotated, for quantification): %s", args.out_rescored)
    return 0


if __name__ == "__main__":
    sys.exit(main())
