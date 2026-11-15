# FINDINGS (Problem C)

## 1. Data defects and traps
- **Multiallelic loci.** The germline resource has a different alt at the same position for several sites: chr17:7577538 (resource G>C, pileup G>A), chr17:7578406 (C>G vs C>A), chr7:5050606 (C>G vs C>T), chr12:8700000 (A>T vs A>G), chr3:2680400 (C>G vs C>T), chr12:8333120 (G>T vs G>A), chr3:2560300 (A>C vs A>G), chr12:8120400 (T>A vs T>C). Matching on chrom+pos only would assign the wrong AF; the exact 4-field match falls back to 1e-06. Detection: key set difference between pileup and resource.
- **Duplicate position, two alts.** chr7:5350700 has C>T and C>G. repeat_annotation.tsv is keyed on (chrom, pos) only, so it has two rows with the same key; a dict build silently keeps one. Harmless here (both "no") but a failure mode.
- **BED coordinates.** BED is 0-based half-open; pileup positions are 1-based. Naive `start <= pos <= end` wrongly keeps chr1:1000000, chr3:2500000, chr7:5000000, chr12:8000000 and chr17:7570000 (all sit exactly on a BED start). Correct test: start <= pos-1 < end.
- **Denominator for TMB.** BED length 2.01 Mb vs callable_bases 1.18644 Mb; the spec says callable territory.
- **Depth mismatch.** Normal depth is much lower than tumour depth (e.g. chr1:1061877 normal 18, chr1:1050120 tumour 42), so low_depth catches sites before PoN or other filters can. chr1:1050120 is also in the panel of normals (3 of 50), but filter 2 fires first.
- Pileup arithmetic (ref+alt = depth, fwd+rev = alt) is consistent for all 48 rows; no duplicate keys.

## 2. normal_type and filter 6
`normal_type` = "tumour-adjacent tissue, same surgical specimen". Adjacent tissue usually contains infiltrating tumour cells, so true somatic variants can show normal alt reads. Filter 6 (normal alt >= 3 and 0.02 <= normal_vaf < 0.20) assumes a clean normal where such reads mean an artifact, so it is more likely to reject real somatic variants here, e.g. chr12:8203311 (normal VAF 0.05, 4 alt reads, af_used 2e-5). Filter 4 can reject real somatic variants the same way: chr7:5350700 C>G (normal VAF 0.25) is removed as germline. Tradeoff: loosening filter 6 or using alt >= 3 relative to the tumour VAF (a ratio test, e.g. normal_vaf < 0.1 x tumour VAF) keeps these but lets real low-level germline and normal-borne artifacts through; a blood normal would remove the problem and allow a stricter rule.

## 3. Threshold evaluation
- **Admits an invalid variant:** filter 9 requires alt_count >= 10. chr3:2800119 T>C has 9 alt reads, all forward (strand_balance 0), sits in an annotated repeat (unit 4, 3 units), and still passes. It is an SNV, so filter 12 (indels only) never applies. Lowering the gate to alt_count >= 5, or using a binomial strand test, would catch it.
- **Rejects a plausible real somatic variant:** filter 7 on chr17:7577120 C>T. The germline resource gives af 0.47, so expected_contam_alt = 0.045 x 320 x 0.47 = 6.77 and the 9 alt reads fall below 2 x 6.77. But the normal has 0 of 90 alt reads, which excludes a germline het, and a TP53 hotspot at 0.47 AF is implausible; the resource row is likely wrong or mis-keyed (identical 0.47 also on chr17:7577538 G>C). A normal-aware contamination test would rescue it. chr7:5350700 C>G (filter 4) is the second candidate; see section 2.
- 
