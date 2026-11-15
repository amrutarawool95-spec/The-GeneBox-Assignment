# Problem C: sample run (TGB-SOM-014)

## Command
```
python3 src/problemC.py problemC_somatic out
```
Standard library only, no randomness or network. Two runs gave byte-identical output (checked with `diff -r`).

## Summary metrics
```
{"n_sites_evaluated": 48, "n_pass": 23,
 "filter_counts": {"off_target": 6, "low_depth": 2, "min_alt": 2, "germline": 3, "panel_of_normals": 2,
   "normal_artifact": 1, "contamination": 1, "orientation_bias": 2, "strand_bias": 1, "mapping_quality": 1,
   "base_quality": 1, "short_tandem_repeat": 2, "read_position": 1, "PASS": 23},
 "tmb_per_mb": 19.386, "median_pass_vaf": 0.146341}
```
Callable territory = 1,186,440 bases (sample_meta.json), not the 2.01 Mb BED length.

## Result
Of 48 candidates, 23 pass and 25 are removed. Six sites are off-target: the BED is 0-based half-open, so five sites that sit exactly on a BED start (chr1:1000000, chr3:2500000, chr7:5000000, chr12:8000000, chr17:7570000) fall outside, and so does chr1:1400001, one past a BED end. The rest are removed by germline (3), low_depth (2), min_alt (2), panel_of_normals (2), orientation_bias (2), short_tandem_repeat (2), and one each for normal_artifact, contamination, strand_bias, mapping_quality, base_quality and read_position. Median PASS VAF is 0.146.

## Is TMB = 19.386 mut/Mb suitable for clinical reporting? No.
- **Small denominator.** Only 1.19 Mb is callable, and the count is 23. A Poisson 95% interval on 23 events is 12.3 to 29.1 mut/Mb. That whole range is above the commonly used 10 mut/Mb cutoff, but it is wide, so the number is imprecise and a reported value would carry about 8 mut/Mb of uncertainty from counting alone. The BED length (2.01 Mb) would give 11.4, so the denominator choice alone moves the value by 41%.
- **Artifacts still inside n_pass.** chr3:2800119 (9 forward, 0 reverse alt reads, in a repeat) passes only because strand_bias needs alt_count >= 10. Low-VAF sites at 2 to 4% (chr17:7577538, chr17:7578406, chr3:2560300, chr12:8120400, chr3:2680400, chr12:8333120) sit near the detection limit with 6 to 11 alt reads.
- **Germline leakage.** Absent alleles get af_used = 1e-06, so germline filtering relies on the normal alone, and the normal is tumour-adjacent tissue (see FINDINGS.md).
- **No UMIs**, so PCR and oxidation errors cannot be separated from real variants. The assay has not been validated against a TMB reference such as WES.
- **Not exome-scale.** A 1.19 Mb hybrid-capture panel can give a different TMB than whole-exome sequencing (panels are often enriched for cancer genes), and no calibration against WES is provided.
