#!/usr/bin/env python3
"""
Problem C2: somatic candidate filtration benchmark.

Usage:
    python3 src/problemC.py <input_dir> <out_dir>

<input_dir> must contain:
    tumor_pileup.tsv, normal_pileup.tsv, germline_resource.tsv, panel_of_normals.tsv,
    repeat_annotation.tsv, target_intervals.bed, sample_meta.json

Writes to <out_dir>:
    problemC_calls.tsv, problemC_summary.json

Design rules (from the assignment):
  * Deterministic: standard library only, no randomness, no network, no dependence on
    file row order or column order (columns are read by lowercase name).
  * No manual edits: every correction is done in code and is reported on stderr.
  * Fails loudly on ambiguous input (conflicting duplicate rows, missing columns)
    instead of guessing.
"""
import bisect
import csv
import json
import statistics
import sys
from pathlib import Path

# ----------------------------------------------------------------------------
# Filter hierarchy: first matching rule wins; rank = 1-based position here.
# ----------------------------------------------------------------------------
FILTERS = [
    "off_target", "low_depth", "min_alt", "germline", "panel_of_normals",
    "normal_artifact", "contamination", "orientation_bias", "strand_bias",
    "mapping_quality", "base_quality", "short_tandem_repeat", "read_position", "PASS",
]
KEY = ("chrom", "pos", "ref", "alt")

TUMOR_COLS = KEY + ("depth", "alt_count", "alt_fwd", "alt_rev",
                    "mean_alt_bq", "mean_alt_mq", "alt_read_pos_median")
NORMAL_COLS = KEY + ("depth", "alt_count")
GERM_COLS = KEY + ("af",)
PON_COLS = KEY + ("n_normals_with_alt",)
REPEAT_COLS = ("chrom", "pos", "in_tandem_repeat", "repeat_unit_length", "repeat_n_units")

WARNINGS = []


def warn(msg):
    WARNINGS.append(msg)


# ----------------------------------------------------------------------------
# Reading helpers
# ----------------------------------------------------------------------------
def read_table(path, required):
    """Read a TSV by column name (order does not matter). utf-8-sig drops a BOM."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path.name}: file is empty")
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path.name}: missing required column(s) {missing}; "
                             f"found {reader.fieldnames}")
        rows = []
        for raw in reader:
            if raw is None or all((v or "").strip() == "" for v in raw.values()):
                continue                      # blank line
            rows.append({k: (v or "").strip() for k, v in raw.items() if k is not None})
        return rows


def num(value):
    """float, or None for NA / empty / non-numeric."""
    if value is None:
        return None
    v = value.strip()
    if v == "" or v.upper() in ("NA", "NAN", "NULL", "."):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def make_key(row):
    return (row["chrom"], int(float(row["pos"])), row["ref"].upper(), row["alt"].upper())


def index_rows(rows, key_fn, name):
    """Dict by key. Identical duplicates are collapsed; conflicting duplicates are fatal."""
    out = {}
    for r in rows:
        k = key_fn(r)
        if k in out:
            if out[k] != r:
                raise ValueError(f"{name}: conflicting duplicate rows for key {k}")
            warn(f"{name}: exact duplicate row for {k} collapsed")
            continue
        out[k] = r
    return out


# ----------------------------------------------------------------------------
# BED: 0-based, half-open. A 1-based position p is inside iff start <= p-1 < end.
# ----------------------------------------------------------------------------
def load_bed(path):
    by_chrom = {}
    with open(path, encoding="utf-8-sig") as fh:
        for n, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            parts = s.split()
            if len(parts) < 3:
                raise ValueError(f"{path.name} line {n}: need chrom, start, end")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end <= start:
                raise ValueError(f"{path.name} line {n}: end <= start")
            by_chrom.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, ivs in by_chrom.items():          # merge overlaps so lookup is exact
        ivs.sort()
        cur = [ivs[0][0], ivs[0][1]]
        out = []
        for s, e in ivs[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                out.append(tuple(cur))
                cur = [s, e]
        out.append(tuple(cur))
        merged[chrom] = out
    return merged


def make_on_target(merged):
    starts = {c: [s for s, _ in ivs] for c, ivs in merged.items()}

    def on_target(chrom, pos1):
        ivs = merged.get(chrom)
        if not ivs:
            return False
        p0 = pos1 - 1
        i = bisect.bisect_right(starts[chrom], p0) - 1
        return i >= 0 and ivs[i][0] <= p0 < ivs[i][1]
    return on_target


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    inp, out = Path(argv[1]), Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)

    meta = json.loads((inp / "sample_meta.json").read_text(encoding="utf-8-sig"))
    for k in ("af_of_alleles_not_in_resource", "contamination_fraction"):
        if k not in meta:
            raise ValueError(f"sample_meta.json: missing '{k}'")
    af_default = float(meta["af_of_alleles_not_in_resource"])
    contam_frac = float(meta["contamination_fraction"])

    tumor = read_table(inp / "tumor_pileup.tsv", TUMOR_COLS)
    normal = index_rows(read_table(inp / "normal_pileup.tsv", NORMAL_COLS), make_key, "normal_pileup")
    germ = index_rows(read_table(inp / "germline_resource.tsv", GERM_COLS), make_key, "germline_resource")
    pon = index_rows(read_table(inp / "panel_of_normals.tsv", PON_COLS), make_key, "panel_of_normals")
    # repeat annotation is keyed on (chrom, pos) only, as supplied
    repeat = index_rows(read_table(inp / "repeat_annotation.tsv", REPEAT_COLS),
                        lambda r: (r["chrom"], int(float(r["pos"]))), "repeat_annotation")
    merged_bed = load_bed(inp / "target_intervals.bed")
    on_target = make_on_target(merged_bed)

    if "callable_bases" in meta:
        callable_bases = float(meta["callable_bases"])
    else:
        callable_bases = float(sum(e - s for ivs in merged_bed.values() for s, e in ivs))
        warn("sample_meta.json has no callable_bases: using total BED length for TMB")
    if callable_bases <= 0:
        raise ValueError("callable_bases must be > 0")

    tumor_keys = set()
    rows = []
    for t in tumor:
        k = make_key(t)
        if k in tumor_keys:
            raise ValueError(f"tumor_pileup: duplicate site {k}")
        tumor_keys.add(k)
        chrom, pos, ref, alt = k

        depth = int(float(t["depth"]))
        alt_count = int(float(t["alt_count"]))
        alt_fwd, alt_rev = num(t["alt_fwd"]), num(t["alt_rev"])
        if alt_fwd is not None and alt_rev is not None and alt_fwd + alt_rev != alt_count:
            warn(f"{k}: alt_fwd + alt_rev != alt_count")

        n = normal.get(k)                         # unmatched site -> normal depth 0
        n_depth = int(float(n["depth"])) if n else 0
        n_alt = int(float(n["alt_count"])) if n else 0
        n_pon = int(float(pon[k]["n_normals_with_alt"])) if k in pon else 0

        vaf = alt_count / depth if depth > 0 else 0.0
        normal_vaf = n_alt / n_depth if n_depth > 0 else 0.0

        if k in germ:
            af_used = float(germ[k]["af"])        # exact (chrom,pos,ref,alt) match only
        else:
            af_used = af_default
            if any(g[:2] == k[:2] for g in germ):
                warn(f"{k}: resource has a different alt at this position -> default AF used")

        exp_contam = contam_frac * depth * af_used
        if alt_fwd is not None and alt_rev is not None and alt_count > 0:
            strand_balance = min(alt_fwd, alt_rev) / alt_count
        else:
            strand_balance = None
        is_indel = len(ref) != len(alt)

        rep = repeat.get((chrom, pos))
        rep_yes = bool(rep) and rep["in_tandem_repeat"].lower() == "yes"
        unit = num(rep["repeat_unit_length"]) if rep else None
        n_units = num(rep["repeat_n_units"]) if rep else None

        mq, bq, rp = num(t["mean_alt_mq"]), num(t["mean_alt_bq"]), num(t["alt_read_pos_median"])

        # A condition that cannot be evaluated (missing number) is treated as not matching.
        def lt(x, limit):
            return x is not None and x < limit

        conds = [
            not on_target(chrom, pos),                                           # 1 off_target
            depth < 50 or n_depth < 20,                                          # 2 low_depth
            alt_count < 5 or vaf < 0.02,                                         # 3 min_alt
            af_used >= 1e-5 and normal_vaf >= 0.20,                              # 4 germline
            n_pon >= 2,                                                          # 5 panel_of_normals
            n_alt >= 3 and 0.02 <= normal_vaf < 0.20,                            # 6 normal_artifact
            alt_count < 2 * exp_contam,                                          # 7 contamination
            (ref, alt) in (("C", "A"), ("G", "T")) and lt(strand_balance, 0.20), # 8 orientation_bias
            alt_count >= 10 and lt(strand_balance, 0.10),                        # 9 strand_bias
            lt(mq, 30),                                                          # 10 mapping_quality
            lt(bq, 20),                                                          # 11 base_quality
            is_indel and rep_yes and unit is not None and n_units is not None
                and unit <= 2 and n_units >= 6,                                  # 12 short_tandem_repeat
            lt(rp, 10),                                                          # 13 read_position
            True,                                                                # 14 PASS
        ]
        idx = next(i for i, c in enumerate(conds) if c)
        rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "depth": depth,
                     "alt_count": alt_count, "vaf": vaf, "normal_vaf": normal_vaf,
                     "af_used": af_used, "expected_contam_alt": exp_contam,
                     "filter": FILTERS[idx], "filter_rank": idx + 1})

    if not rows:
        raise ValueError("tumor_pileup has no sites")
    unseen = sorted({r["chrom"] for r in rows} - set(merged_bed))
    if unseen:
        warn(f"chromosomes in tumor_pileup absent from BED (naming mismatch?): {unseen}")

    # ---- calls.tsv: one row per input site, input order preserved ----
    with open(out / "problemC_calls.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chrom", "pos", "ref", "alt", "depth", "alt_count", "vaf", "normal_vaf",
                    "af_used", "expected_contam_alt", "filter", "filter_rank"])
        for r in rows:
            w.writerow([r["chrom"], r["pos"], r["ref"], r["alt"], r["depth"], r["alt_count"],
                        f"{r['vaf']:.6f}", f"{r['normal_vaf']:.6f}", f"{r['af_used']:.8f}",
                        f"{r['expected_contam_alt']:.6f}", r["filter"], r["filter_rank"]])

    # ---- summary.json ----
    passed = [r for r in rows if r["filter"] == "PASS"]
    counts = {f: sum(1 for r in rows if r["filter"] == f) for f in FILTERS}
    summary = {
        "n_sites_evaluated": len(rows),
        "n_pass": len(passed),
        "filter_counts": {f: c for f, c in counts.items() if c > 0},
        "tmb_per_mb": round(len(passed) / (callable_bases / 1e6), 3),
        "median_pass_vaf": round(statistics.median(r["vaf"] for r in passed), 6) if passed else 0.0,
    }
    with open(out / "problemC_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary))

    for msg in sorted(set(WARNINGS)):             # sorted so stderr is deterministic too
        print(f"WARNING: {msg}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
  
