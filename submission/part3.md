# Part 3: Mathematical Formulations

All final numbers were checked with a short Python script (scipy), and cross-checked against an independent solution.

---

## Question 1 (Q3.1): PGx assay at low prevalence

**Given:** prevalence π = 1/10,000 = 1e-4, sensitivity Se = 0.99, specificity Sp = 0.999. Target PPV = 0.95 for part (b).

**Formula:**
- PPV = TP / (TP + FP), with TP = Se·π and FP = (1 − Sp)(1 − π)
- For a target PPV: FPR = TP·(1/PPV − 1) / (1 − π), where FPR = 1 − Sp

**Derivation:**

(a) Per person tested:
- TP = 0.99 × 1e-4 = 9.9e-5
- FP = 0.001 × 0.9999 = 9.999e-4
- PPV = 9.9e-5 / (9.9e-5 + 9.999e-4) = 9.9e-5 / 1.0989e-3 = 0.09009

So about 91% of positive calls are false, even though the assay looks very good.

(b) Set PPV = 0.95 and keep TP fixed at 9.9e-5:
- Allowed FP = 9.9e-5 × (1/0.95 − 1) = 9.9e-5 × 0.052632 = 5.2105e-6
- FPR = 5.2105e-6 / 0.9999 = 5.2111e-6
- Required specificity = 1 − 5.2111e-6 = 0.99999479
- n = 1 / FPR = 0.949905 / 4.95e-6 = 191,900

**Final Answer:**
(a) PPV ≈ **0.0901 (9.0%)**
(b) **One false positive per 191,900 negatives** (specificity ≈ 0.9999948)

---

## Question 2 (Q3.2): Whole-genome sequencing at 5×, Poisson

**Given:** genome size G = 3.1 × 10⁹ bp, mean coverage λ = 5, depth D ~ Poisson(5). At a heterozygous site each read is allele A or B with probability 0.5.

**Formula:**
- P(D = 0) = e^(−λ), so expected unsequenced bases = G · e^(−λ)
- P(both alleles | depth d) = 1 − 2·(0.5)^d for d ≥ 1 (this is 0 at d = 1)
- P(both) = Σ over d ≥ 1 of [1 − 2·(0.5)^d] · e^(−5) 5^d / d!

**Derivation:**

(a) e^(−5) = 0.006738, so unsequenced bases = 3.1e9 × 0.006738 = 2.0888e7.

(b) Split the sum into two parts, using Σ_{d≥1} x^d/d! = e^x − 1:
- Σ_{d≥1} Poisson(d; 5) = 1 − e^(−5)
- Σ_{d≥1} 2·(0.5)^d · Poisson(d; 5) = 2·e^(−5)·(e^(2.5) − 1) = 2·(e^(−2.5) − e^(−5))

P(both) = (1 − e^(−5)) − 2·(e^(−2.5) − e^(−5)) = 1 − 2e^(−2.5) + e^(−5)

= 1 − 2(0.082085) + 0.006738 = 0.842568

Check by Poisson thinning: reads of allele A ~ Poisson(2.5) and reads of allele B ~ Poisson(2.5), independent, so P(both) = (1 − e^(−2.5))² = 0.917915² = 0.842568. ✓

**Final Answer:**
(a) **≈ 2.09 × 10⁷ bases** (about 20.9 Mb, 0.67% of the genome)
(b) **≈ 0.8426 (84.3%)**. Using a fixed depth of 5 would wrongly give 0.9375.

---

## Question 3 (Q3.3): Relative vs absolute abundance

**Given:** Sample A: relative abundance 4.0%, total load 1.0 × 10¹¹ cells/g. Sample B: relative abundance 2.0%, total load 2.5 × 10¹¹ cells/g.

**Formula:** absolute abundance = relative abundance × total load; fold change = B / A.

**Derivation:**
- A = 0.04 × 1.0e11 = 4.0e9 cells/g
- B = 0.02 × 2.5e11 = 5.0e9 cells/g
- Fold change = 5.0e9 / 4.0e9 = 1.25 (difference = +1.0e9 cells/g)
- Check: relative abundance halved (×0.5) but total load rose ×2.5, and 0.5 × 2.5 = 1.25.

**Final Answer:** **1.25-fold increase in sample B** (4.0e9 to 5.0e9 cells/g), even though the relative abundance fell from 4% to 2%.

---

## Question 4 (Q3.4): "Always hom-ref" baseline caller

**Given:** the caller always predicts 0/0. Hardy–Weinberg equilibrium holds, minor allele frequency q.

**Formula:**
- Concordance at one locus = P(true genotype is 0/0) = (1 − q)²
- (a) average over q ~ Uniform(0, 0.5], density 1/0.5 = 2
- (b) solve (1 − q)² = 0.90

**Derivation:**

(a) E[(1 − q)²] = 2 ∫₀^0.5 (1 − 2q + q²) dq = 2 [q − q² + q³/3] from 0 to 0.5
= 2 (0.5 − 0.25 + 0.125/3) = 2 (7/24) = 7/12

(b) 1 − q = √0.90 = 0.948683, so q = 1 − 0.948683 = 0.051317.

**Final Answer:**
(a) **7/12 ≈ 0.5833**
(b) **q ≈ 0.0513 (about 5.1%)**

---

## Question 5 (Q3.5): Low-abundance organism in a metagenome

**Given:** genome size 3 Mb = 3 × 10⁶ bp, relative abundance f = 0.01% = 1e-4, total yield 5 Gb for part (b). Assumption: abundance equals the share of sequenced bases (no genome-size correction).

**Formula:**
- Mean coverage = (total yield × f) / genome size
- Fraction covered by ≥ 1 read = 1 − e^(−λ) (Poisson)

**Derivation:**

(a) For 1× coverage we need 3e6 bases from this organism. Total yield = 3e6 / 1e-4 = 3e10 bases.

(b) Bases from the organism = 5e9 × 1e-4 = 5e5. λ = 5e5 / 3e6 = 0.1667. Fraction covered = 1 − e^(−0.1667) = 1 − 0.84648 = 0.15352.

**Final Answer:**
(a) **30 Gb**
(b) **≈ 0.1535 (15.4%)** of the genome is covered at least once.

---

## Question 6 (Q3.6): Somatic caller, ≥ 3 alt reads, 5% VAF

**Given:** detection needs ≥ 3 alt reads, VAF = 0.05, depth 60×, 20 subclonal variants, target power 95%.

**Formula:**
- Alt reads X ~ Binomial(n, 0.05)
- Power = P(X ≥ 3) = 1 − [P(0) + P(1) + P(2)]
- Expected misses = N × (1 − power)

**Derivation:**

(a) At n = 60:
- P(0) = 0.95⁶⁰ = 0.04607
- P(1) = 60 × 0.05 × 0.95⁵⁹ = 0.14548
- P(2) = C(60,2) × 0.05² × 0.95⁵⁸ = 0.22588
- Sum = 0.41744, so power = 1 − 0.41744 = 0.5826

(b) Each variant is missed with probability 0.41744. For 20 variants: 20 × 0.41744 = 8.35.

(c) Find the smallest n where P(X ≤ 2) ≤ 0.05:
- n = 123: P(X ≤ 2) = 0.05142 (power 0.9486, not enough)
- n = 124: P(X ≤ 2) = 0.04953 (power 0.9505, enough)

**Final Answer:**
(a) Power **≈ 0.5826 (58.3%)**
(b) **≈ 8.35 missed variants** out of 20 (about 11.65 detected)
(c) **Minimum depth = 124×**

---

## Question 7 (Q3.7): "Reads" vs gigabases for 2 × 150 bp

**Given:** paired-end 2 × 150 bp. Options "5 million reads" advertised as 1.5 Gb and "18 million reads" advertised as 5.4 Gb.

**Formula:**
- If a "read" means a read pair: yield = N × 2 × 150 bp
- If a "read" means a single read: yield = N × 150 bp

**Derivation:**
- Read pairs: 5e6 × 300 = 1.5 Gb, and 18e6 × 300 = 5.4 Gb (matches the advertised numbers exactly)
- Single reads: 5e6 × 150 = 0.75 Gb, and 18e6 × 150 = 2.7 Gb (does not match)

**Final Answer:** The definition that matches is **read pairs**. Reading the same labels as single reads gives half the yield, an error factor of **exactly 2×** (0.75 Gb instead of 1.5 Gb, and 2.7 Gb instead of 5.4 Gb).

---

## Question 8 (Q3.8): Counting diplotypes with 9 alleles

**Given:** 9 non-deletion star alleles, copy number 2 and copy number 3. A diplotype is an unordered multiset of alleles.

**Formula:**
- Number of multisets of size n from 9 types = C(9 + n − 1, n)
- Orderings of a multiset = n! / (product of the factorials of the repeat counts)

**Derivation:**

(a)
- Copy number 2: C(10, 2) = 45
- Copy number 3: C(11, 3) = 165

(b) 9ⁿ counts ordered tuples (which copy carries which allele). A diplotype ignores order, so tuples with the same alleles in a different order collapse into one.
- n = 2: 81 = 9 (same allele twice, 1 ordering each) + 36 pairs × 2 orderings. Unordered count = 9 + 36 = 45.
- n = 3: 729 = 9 × 1 (all same) + 72 × 3 (two same, one different) + 84 × 6 (all different). Unordered count = 9 + 72 + 84 = 165.

(c) For {\*2, \*2, \*41}: 3! / 2! = 3 (the \*41 can sit on copy 1, 2 or 3).

**Final Answer:**
(a) **45** diplotypes at copy number 2 and **165** at copy number 3
(b) 9ⁿ counts ordered assignments and overcounts by the number of orderings of each multiset (up to n!)
(c) **3** ordered assignments

---

## Question 9 (Q3.9): "Population mean dosage" caller (r² = 0)

**Given:** the caller gives every sample dosage 2p, then hard-calls by rounding to the nearest integer. Hardy–Weinberg equilibrium. Alternate allele frequency p = 0.01 and p = 0.30.

**Formula:** HWE genotype frequencies: P(0) = (1 − p)², P(1) = 2p(1 − p), P(2) = p². Concordance = frequency of the genotype the caller always outputs.

**Derivation:**
- p = 0.01: 2p = 0.02 rounds to 0, so every sample is called 0/0. Concordance = (0.99)² = 0.9801.
- p = 0.30: 2p = 0.60 rounds to 1, so every sample is called 0/1. Concordance = 2 × 0.3 × 0.7 = 0.42.

**Final Answer:**
- **0.9801 at p = 0.01** and **0.42 at p = 0.30**
- Why it misleads: the caller has zero information (dosage r² = 0), yet it scores 98% at a rare variant simply because most people are homozygous reference. Overall concordance is dominated by the common genotype, so it can make a useless caller look excellent, and it moves with allele frequency rather than with skill. Validation studies should also report non-reference concordance, per-genotype sensitivity and precision, and dosage r², stratified by MAF.

---

## Question 10 (Q3.10): Epigenetic clock with missing probes

**Given:** 68 probes, 14 missing on the newer array, summed weights of the missing probes Σw = −18.6, cohort mean beta = 0.62, imputation beta = 0.50 for part (b). Assumption: weights and betas are on the scale of the clock output (e.g. years).

**Formula:**
- Prediction = intercept + Σ wᵢ·βᵢ
- Contribution of the missing probes ≈ Σw × mean β
- Bias = (prediction with the change) − (original prediction)

**Derivation:**

(a) The true contribution of the missing probes = −18.6 × 0.62 = −11.532. Excluding them sets it to 0.
Bias = 0 − (−11.532) = +11.532.

(b) Imputed contribution = −18.6 × 0.50 = −9.30.
Bias = −9.30 − (−11.532) = +2.232.

(c) Imputing with one constant still leaves a bias, and excluding probes leaves a large one. The better options are to refit or recalibrate the clock on the probes common to both arrays, or to impute each missing probe with its own training-cohort mean (or a model based on correlated probes), then validate against the original array.

**Final Answer:**
(a) **+11.53** (positive, so predictions are overestimated)
(b) **+2.23** (still an overestimate, about 5× smaller)
(c) **Preferred strategy:** refit or recalibrate on the shared probes, or impute per-probe from training data. **Report:** how many probes were missing (14 of 68, about 20.6%) and which ones; the imputation or refit method; the estimated bias or error (MAE, correlation against the original array on overlapping samples); and a warning that values are not directly comparable to the original-array clock.
