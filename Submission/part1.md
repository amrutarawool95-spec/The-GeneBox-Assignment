# Part 1: Technical short answers

**Selected questions: 1.1, 1.2, 1.5, 1.6, 1.7**

---

## 1.1

`GGB9712_SGB15212` is a MetaPhlAn 4 label. SGB means "species-level genome bin" (genomes grouped at about 95% ANI, roughly 5% genomic distance) and GGB means "genus-level genome bin" (a wider group, about 15% distance), so the label reads "species bin 15212, sitting inside genus bin 9712". The 1.2% is just the estimated relative abundance of that bin; it is not part of the name. The species is "unassigned" because this bin has no genome from a formally named species in it. It was built only from metagenome-assembled genomes (an uncharacterised SGB, or uSGB), so there is no species name to pass on. What we can honestly report is that a distinct species-level lineage is present at about 1.2% of the profiled community, identified by its bin ID, plus any higher ranks the report actually gives. We should not invent a Latin name, call it a "novel species", or match it to a named species in another da

## 1.2

MetaPhlAn counts reads that map to clade-specific marker genes, turns them into coverage (reads per marker length), and divides each clade's coverage by the total coverage of all detected clades. So its denominator is "all detected microbial clades" (unmapped and host reads are left out), every rank adds up to 100%, and it is closer to a share of genome copies or cells than a share of reads. [1][2] A k-mer classifier like Kraken 2 reports the reads in a clade divided by all input reads, unclassified ones included. If you renormalise to classified reads only, that has to be stated. [3][4] Because k-mer fractions follow the amount of DNA, big genomes and host DNA take a larger share, and ambiguous reads get pushed up to higher ranks (LCA), which shifts the kingdom numbers further. So two tables can both say "relative abundance" and still disagree; they only line up after both are put on the same denominator.

---

## 1.5

Both diplotypes have the same alleles, but the activity score counts every copy, so it matters which allele is duplicated. With the 2020 CPIC/DPWG consensus values (`*1` = 1.0, `*41` = 0.5), `*1/*41x2` = 1 + 2×0.5 = **2.0 (normal metabolizer)**, while `*1x2/*41` = 2×1 + 0.5 = **2.5 (ultrarapid)**. CPIC later updated its allele table (March 2023) and lowered `*41` (and `*9`) to 0.25. With that version the pair becomes 1 + 2×0.25 = **1.5** and 2×1 + 0.25 = **2.25**, both normal, but still 0.75 apart. For `*1/*2x2` vs `*1x2/*2` nothing changes: `*1` and `*2` are both normal function (1.0), so each totals 3×1.0 = **3.0 (ultrarapid)** in either version. The 2020 consensus is also where `*10` was cut from 0.5 to 0.25, so always say which table version you used; I used the 2020 values here because the assignment's `allele_definitions.tsv` lists `*41` = 0.5. [5][6][7]

---

## 1.6

Bracken reads the Kraken report by column position and expects the standard layout: % of reads, clade reads, taxon reads, rank, taxid, name. Running Kraken 2 with `--report-minimizer-data` adds two extra columns (total and distinct minimizer counts) after the taxon-reads column, so rank and taxid slide to columns 6 and 7. Bracken then reads the wrong fields, which can crash it or give meaningless abundances, and minimizer counts are not read counts anyway, so even a repaired file would give wrong numbers. The fix is to rerun Kraken 2 without that flag and run Bracken with a k-mer distribution file built for your database and read length. If you only have the minimizer report, drop the extra columns first. [3][4][8]

```
kraken2 --db DB --report sample.kreport --output sample.kraken sample.fastq.gz
bracken-build -d DB -k KMER_LEN -l READ_LEN          # once per database
bracken -d DB -i sample.kreport -o sample.bracken -r READ_LEN -l S
# only if you already have a minimizer report:
cut -f1-3,6-8 sample.minimizer.kreport > sample.kreport
```

---

## 1.7

If a locus has no record in an all-sites VCF, the file simply says nothing about it: the caller may have had too little coverage, filtered the site, or never looked, so you cannot treat it as reference. An explicit `0/0` record means the caller looked and saw the reference allele on both copies, so that position counts as reference with evidence behind it. The quality number that goes with it is **GQ** (genotype quality), the phred-scaled chance that the genotype call is wrong, usually read together with depth (DP, or MIN_DP in gVCF blocks). Standard QUAL thresholds do not work here because QUAL measures how sure we are that a variant exists, and at a reference call it is near zero or just ".", by design. A "QUAL ≥ 30" filter would therefore throw away every confident reference call, which is why these are filtered on GQ and depth instead. [9][10]

---

## References

1. Blanco-Míguez A, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nat Biotechnol* 2023. doi:10.1038/s41587-023-01688-w
2. Segata N, et al. Metagenomic microbial community profiling using unique clade-specific marker genes. *Nat Methods* 2012. doi:10.1038/nmeth.2066
3. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. *Genome Biol* 2019. doi:10.1186/s13059-019-1891-0
4. Kraken 2 manual (report format and `--report-minimizer-data` option). https://github.com/DerrickWood/kraken2/wiki/Manual
5. Caudle KE, et al. Standardizing CYP2D6 genotype to phenotype translation: consensus recommendations from CPIC and DPWG. *Clin Transl Sci* 2020. doi:10.1111/cts.12692
6. PharmGKB/CPIC. CYP2D6 allele function update (20 March 2023; `*9` and `*41` activity value lowered to 0.25). https://pharmgkb.blogspot.com/2023/03/cyp2d6-allele-function-update.html
7. Assignment file `problemB_cyp2d6/allele_definitions.tsv`, rows `*1` (1.0), `*2` (1.0), `*41` (0.5).
8. Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci* 2017. doi:10.7717/peerj-cs.104 (and the Bracken README, https://github.com/jenniferlu717/Bracken)
9. Danecek P, et al. The variant call format and VCFtools. *Bioinformatics* 2011. doi:10.1093/bioinformatics/btr330
10. VCF specification v4.3 (QUAL and GQ definitions). https://samtools.github.io/hts-specs/VCFv4.3.pdf
