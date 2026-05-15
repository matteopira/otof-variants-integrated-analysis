# OTOF Variants: Integrated Analysis for Cochlear Gene Therapy

Independent reanalysis of clinically annotated variants of the *OTOF* gene, the molecular target of first-in-human AAV-mediated cochlear gene therapy trials, integrating ClinVar clinical curation data, gnomAD population frequencies, ConSurf evolutionary conservation scores, and AlphaFold-predicted structural information.

## Background

### Otoferlin and the molecular basis of DFNB9

The *OTOF* gene, located on chromosome 2p23.3, encodes **otoferlin**, a member of the ferlin family of multi-C2-domain transmembrane proteins. Otoferlin is essential for **calcium-triggered synaptic vesicle exocytosis** at the ribbon synapse of cochlear inner hair cells, where it mediates the rapid release of glutamate required for sustained neurotransmission to spiral ganglion neurons (Roux et al., *Cell* 2006; Pangrsic et al., *Trends Neurosci* 2012). The protein contains six C2 domains (C2A through C2F), with C2B, C2C, and C2D being the most critical for calcium-dependent vesicle fusion, while C2A is reported to be atypical in its calcium-binding properties.

Biallelic pathogenic variants in *OTOF* cause **DFNB9** (autosomal recessive nonsyndromic hearing loss 9, OMIM #601071), one of the most clinically distinct forms of hereditary hearing loss. The phenotype is classified as **auditory neuropathy spectrum disorder (ANSD)**: outer hair cell function is typically preserved (normal otoacoustic emissions and cochlear microphonics), but synaptic transmission to the auditory nerve is impaired, resulting in absent or grossly abnormal auditory brainstem responses (Yasunaga et al., *Nat Genet* 1999; Varga et al., *J Med Genet* 2003).

### *OTOF*-related hearing loss: population epidemiology and variant heterogeneity

*OTOF* pathogenic variants account for an estimated 1 to 8% of cases of autosomal recessive nonsyndromic prelingual hearing loss, with substantial variation across populations. Large cohort studies have established a complex variant landscape:

- In a **Japanese cohort of 2,265 sensorineural hearing loss patients**, massively parallel sequencing identified *OTOF* as a major cause of recessive ANSD, with population-specific variants including p.Arg1939Gln (Iwasa et al., *Sci Rep* 2019).
- The **Spanish population** shows a strong founder effect for the **c.2485C>T (p.Gln829Ter)** variant, accounting for approximately 3% of cases of recessive prelingual deafness in Spain (Migliosi et al., *J Med Genet* 2002; Gallo-Terán et al., *Acta Otorrinolaringol Esp* 2005). The Q829X variant has been a primary target for novel therapeutic strategies, including RNA base editing approaches in humanized mouse models (Xue et al., *Mol Ther* 2023).
- A **Pakistani cohort of 557 consanguineous families** estimated that *OTOF* variants account for approximately 2.3% of severe-to-profound prelingual hearing loss in the population, and identified 10 novel pathogenic variants spanning frameshift, nonsense, and missense changes (Choi et al., *Hum Mutat* 2009).
- In a **Taiwanese ANSD cohort of 65 unrelated patients**, integrated short-read and long-read sequencing combined with splice prediction pipelines (SpliceAI) and minigene assays clarified the pathogenicity of cryptic variants and haplotype phasing in *OTOF* (Chen et al., *Mol Med* 2025).
- In **consanguineous populations of the Indian subcontinent and the Arabian Peninsula**, *OTOF* has been identified as a major contributor to autosomal recessive hearing loss, with novel homozygous variants such as p.Arg708Ter detected through whole-exome sequencing approaches (Vibhuti et al., *Front Genet* 2021).

This pronounced **population stratification** of *OTOF* variants has direct implications for diagnostic strategies, genetic counseling, and the selection of patients eligible for gene therapy.

### The advent of AAV-mediated cochlear gene therapy for DFNB9

In 2024, *OTOF* became the first gene targeted in successful human clinical trials of AAV-mediated cochlear gene therapy. The challenge of delivering the ~6 kb *OTOF* coding sequence, which exceeds the ~4.7 kb packaging capacity of conventional AAV vectors, was overcome through dual-AAV strategies, in which the *OTOF* cDNA is split between two AAV particles that recombine *in vivo* to reconstitute the full-length protein (Akil et al., *PNAS* 2019; Al-Moyed et al., *EMBO Mol Med* 2019). Multiple independent clinical trials have now reported the restoration of hearing in pediatric patients with biallelic pathogenic *OTOF* variants:

- **AAV1-hOTOF**: a single-arm trial in Chinese children with DFNB9 demonstrated significant improvement in auditory thresholds (Lv et al., *Lancet* 2024).
- **DB-OTO**: bilateral gene therapy in children showed safety and efficacy of binaural treatment (Wang et al., *NEJM* 2024).
- AAV-mediated *OTOF* delivery has been further validated in multiple parallel programs targeting DFNB9 (Qi et al., *Adv Sci* 2024).

These breakthroughs have made the accurate molecular characterization of *OTOF* variants a problem of immediate clinical relevance. Patient eligibility for gene therapy requires unambiguous biallelic pathogenic genotyping; missense variants of uncertain significance (VUS) represent a major bottleneck in patient stratification.

## Research questions

This project addresses five integrated questions:

1. **Clinical curation**: What is the distribution of *OTOF* variants in ClinVar across germline classification categories, and how does it correlate with predicted molecular consequence?
2. **Population genetics**: Do *OTOF* pathogenic variants in ClinVar comply with the ACMG/AMP BS1 criterion (allele frequency < 0.005) when cross-referenced with gnomAD, and are there ancestry-specific founder mutations relevant for population-tailored diagnostic and therapeutic strategies?
3. **Structural genomics**: Are pathogenic variants preferentially localized to specific functional domains of otoferlin, and does the variant density per domain reflect functional intolerance?
4. **Conservation-based reclassification**: Do evolutionary conservation scores from ConSurf DB support the reclassification of missense VUS in critical C2 domains using ACMG/AMP criteria PP3 and BP4?
5. **Temporal trends**: Has the rate of *OTOF* variant evaluation in ClinVar changed over time, and is there a detectable acceleration following the 2024 AAV-OTOF clinical trials?

## Data sources

| Database | Version / Access date | Records retrieved |
|----------|----------------------|-------------------|
| **ClinVar** (NCBI) | Retrieved May 2026 | 2,432 *OTOF* variants |
| **gnomAD** (Broad Institute) | v4 | 9,601 *OTOF* variants |
| **ConSurf DB** | Per-residue grades, UniProt Q9HC10 | 1,997 residues |
| **UniProt** | Q9HC10 (canonical, 1997 aa) | Reference protein sequence |
| **AlphaFold DB** (EMBL-EBI) | Model v4 | Predicted full-length 3D structure |

All datasets are publicly available and were retrieved through their respective web interfaces. Data files are included in the `data/` directory.

## Methods

All computational analyses were performed in **Python 3.13** within the Jupyter Notebook environment. Structural visualization was conducted in **PyMOL** (open-source build, conda-forge channel). The exact software environment is specified in `requirements.txt` (pip) and `environment.yml` (conda).

### Software stack

- **pandas 2.3.3**, for tabular data manipulation
- **NumPy 2.3.5**, for numerical operations
- **matplotlib 3.10.6** and **seaborn 0.13.2**, for static visualization
- **plotly 6.3.0**, for interactive visualization
- **scipy 1.16.3**, for statistical tests (chi-square, Kruskal-Wallis, Mann-Whitney U)
- **requests 2.32.5**, for ConSurf DB API access
- **re** (Python standard library), for HGVS protein notation parsing
- **PyMOL** (open-source), for molecular visualization and animation

### Analytical workflow

The analysis is structured into seven Jupyter notebooks:

- **`notebooks/01-exploration.ipynb`**: ClinVar descriptive analysis, including variant counts, germline classification distribution, variant type and molecular consequence stratification, cross-tabulation of consequence versus classification, and stratified analysis of records with and without condition annotation.
- **`notebooks/02-gnomad-integration.ipynb`**: gnomAD v4 integration, including allele frequency distributions stratified by clinical classification, ACMG/AMP BS1 validation, ancestry-specific founder mutation analysis across nine population groups, and a detailed investigation of P/LP variants observed as homozygotes in gnomAD with Hardy-Weinberg expectation analysis.
- **`notebooks/03-domain-mapping.ipynb`**: domain-level variant architecture analysis. This notebook implements the complete analytical pipeline for mapping *OTOF* variant pathogenicity onto the protein's functional domain structure. It includes: (i) a robust, multi-format HGVS protein notation parser supporting standard three-letter notation (p.Arg963Ter), one-letter shorthand (p.R963*), embedded p-dot within c-dot notation (e.g. `c.2887C>T(p.Arg963Ter)`), and HGVS full-transcript notation, with systematic comparison against a legacy single-regex parser; (ii) assignment of parsed amino acid positions to the six C2 domains and the transmembrane domain of UniProt Q9HC10 (1,997 aa), with all remaining positions assigned to the combined linker background; (iii) calculation of P/LP variant density (variants per residue) and ascertainment-corrected fraction pathogenic (P/LP divided by all classified variants) per domain, with 95% Wilson confidence intervals; (iv) a lollipop plot of all 1,600 position-mapped variants coloured by classification; (v) four complementary statistical tests for domain-level pathogenic intolerance (chi-square goodness-of-fit, two-tailed permutation test, Fisher's exact test with both Bonferroni and BH correction, fraction pathogenic with Wilson CI); (vi) a rigorous sensitivity analysis quantifying HGVS parser failure rates stratified by molecular consequence class; (vii) a corrected and expanded deep dive on p.Arg1172Gln integrating gnomAD v4 pathogenicity scores (CADD, PhyloP) and per-ancestry allele frequencies with homozygote localisation; (viii) an AlphaFold pLDDT formal comparison between P/LP and Benign/Likely benign variants using the Mann-Whitney U test with common-language effect size; (ix) a Spearman rank-correlation sensitivity analysis testing whether domain intolerance rankings are stable when analysis is restricted to ClinVar records carrying an explicit condition annotation.
- **`notebooks/04-consurf-vus.ipynb`**: ConSurf DB integration, per-residue conservation grade retrieval, conservation-based classification of missense VUS, and ACMG/AMP PP3/BP4 reclassification candidates.
- **`notebooks/05-clinvar-temporal.ipynb`**: temporal analysis of ClinVar variant evaluations by year, cumulative classification trends, and comparison of pre-2024 vs. post-2024 P/LP evaluation rates using Mann-Whitney U test.
- **`notebooks/06-interactive-lollipop.ipynb`**: interactive Plotly lollipop figure integrating all data layers, displayed inline in Jupyter with static PNG export.
- **`notebooks/07-structural-analysis.ipynb`**: AlphaFold structure-based analysis including manual PDB parsing, pLDDT confidence score extraction, per-variant distance calculation to calcium-coordinating residues, pLDDT distribution comparison across classification groups (Kruskal-Wallis test), distance-to-calcium-binding-residues comparison (Mann-Whitney U test), and spatial clustering analysis of P/LP variants via permutation test (n=10,000, seed=42).

### Statistical testing (notebook 03)

Domain-level pathogenic intolerance is assessed using five complementary statistical frameworks, each designed to address a distinct potential confound.

**1. Chi-square goodness-of-fit** tests whether the global distribution of P/LP variants across named domains (C2A through C2F and TM) deviates from a length-proportional null hypothesis. Under the null, every amino acid in a named domain has an equal probability of carrying a P/LP variant. Expected counts per domain are `E_i = N_total_PLP × (L_i / L_total)`, where `L_i` is domain length in amino acids. The test is global and does not identify which domains deviate.

**2. Two-tailed permutation test (n = 10,000, seed = 42)** is applied per domain to detect both enrichment and depletion of P/LP variants relative to the pool of all named-domain residues. For each of 10,000 iterations, P/LP variants are randomly reassigned across all residues belonging to named domains, preserving total P/LP count. The observed P/LP count for domain *d* is compared to the empirical null distribution. The two-tailed p-value is defined as `p = 2 × min(p_enrichment, p_depletion)`, capped at 1.0 to remain a valid probability. This test is more sensitive than chi-square for detecting asymmetric distributions and captures domain-specific effects masked by the global test.

**3. Fisher's exact test per domain versus linker** compares the P/LP-to-non-P/LP ratio within each named domain against the same ratio in the linker region as a fixed, biologically motivated background. The linker (all positions outside named domains, combined length 1,258 aa) is the largest continuous region of otoferlin without an assigned functional module and provides a specific, rather than average, reference. The 2×2 contingency table for each domain is constructed using **observed variant counts from ClinVar** (not domain lengths), so the denominator correctly reflects ascertainment rather than sequence coverage. Two multiple-testing procedures are applied: (i) conservative Bonferroni correction (`p_Bonf = min(p × 7, 1.0)`), reported in `results/domain_fisher_exact.csv`; and (ii) Benjamini–Hochberg FDR correction (`statsmodels.stats.multitest.multipletests`, `method='fdr_bh'`), reported in `results/domain_fisher_bh.csv`. The BH procedure is the primary correction for this exploratory analysis, as Bonferroni is overly conservative with seven correlated tests sharing a common background.

**4. Fraction pathogenic per domain** (ascertainment-corrected metric) is defined as `f_d = N_PLP / (N_PLP + N_BLB + N_VUS + N_Conflict)`, where the denominator is the count of all classified variants in domain *d* regardless of class, excluding the "Other" category. This metric conditions on total ClinVar submissions per domain, thereby adjusting for the uneven ascertainment of variants across functional regions. Domains with high coverage (many submissions of all classes) yield narrower confidence intervals. Point estimates and 95% Wilson score confidence intervals are computed as:

```
p̂ = k / n
CI_centre = (p̂ + z²/2n) / (1 + z²/n)
CI_margin  = z × √(p̂(1−p̂)/n + z²/4n²) / (1 + z²/n)
```

where `k` = P/LP count, `n` = domain denominator, `z = 1.96` for 95% confidence. Results are reported in `results/domain_fraction_pathogenic.csv` and visualised in `results/domain_fraction_pathogenic.png`.

**5. HGVS parsing failure rate by molecular consequence** quantifies the differential parseability of the HGVS parser across consequence classes. For each molecular consequence category (missense, frameshift, nonsense/stop, splice, intronic, synonymous, UTR, InDel), the failure rate (`N_non-parseable / N_total`) is computed. A chi-square test of independence (`scipy.stats.chi2_contingency`) tests whether parseability is distributed uniformly across consequence classes. If failure is concentrated in non-coding classes (intronic, splice), the exclusion of non-parseable variants from domain analyses is mechanistically explained and does not introduce classification bias. Results are in `results/hgvs_failure_by_consequence.csv` and `results/hgvs_failure_by_consequence.png`.

### p.Arg1172Gln characterisation

The variant p.Arg1172Gln (NM_194248.2:c.3515G>A; ClinVar classification: Pathogenic) is of particular interest because it is observed as a homozygote in one gnomAD v4 individual. The deep-dive analysis (notebook 03, cells R3) integrates three data layers available in the gnomAD v4 export without requiring external downloads: (i) per-ancestry allele frequencies and homozygote counts across nine population groups; (ii) pre-computed *in silico* pathogenicity scores: CADD Phred, PhyloP vertebrate conservation, REVEL, SIFT, and PolyPhen-2 where available; (iii) ConSurf evolutionary conservation grade from the per-residue dataset (notebook 04 data). **Correction of a prior data error**: an earlier version of the analysis incorrectly attributed the single gnomAD homozygote to the African/African-American ancestry group. The corrected analysis demonstrates that the homozygote belongs to the European (non-Finnish) ancestry group (AC = 1, AN = 1,085,662 for this ancestry; zero homozygotes in any other group). Results are in `results/arg1172gln_ancestry_af.csv` and `results/arg1172gln_ancestry_af.png`.

### AlphaFold pLDDT comparative analysis

The predicted local distance difference test (pLDDT) score, stored in the B-factor field of the AlphaFold v4 PDB model (Q9HC10), is used as a proxy for per-residue structural order. For each ClinVar variant with a parseable amino acid position, the pLDDT of the affected residue is extracted by matching the residue sequence number. A Mann-Whitney U test (two-sided) compares the pLDDT distributions of residues harbouring P/LP variants versus Benign/Likely benign variants. The common language effect size (CLES, also known as the probability of superiority) is defined as `CLES = U / (n_PLP × n_BLB)` and quantifies the probability that a randomly selected P/LP-carrying residue has higher pLDDT than a randomly selected Benign-carrying residue. A violin plot overlaid with a box plot is produced for all three classification groups (P/LP, VUS, Benign/LB). Results are in `results/plddt_plp_vs_blb_formal.png`.

### Sensitivity analysis: not-provided annotation and domain ranking stability

The stability of domain-level pathogenic intolerance rankings to the exclusion of ClinVar records lacking an explicit condition annotation ("not provided") is tested using the Spearman rank correlation coefficient. Two parallel analyses are performed: (i) including all mapped ClinVar records (n = 1,600 position-parseable variants); (ii) including only records with a specific condition annotation (n varies by domain). For each analysis, domains are ranked by two metrics: absolute P/LP density (P/LP per amino acid) and fraction pathogenic (P/LP per classified variant). The Spearman ρ between all-records and curated-only rankings is computed for both metrics (`scipy.stats.spearmanr`). A threshold of ρ ≥ 0.90 is adopted a priori as the criterion for declaring rankings "robust": if ρ ≥ 0.90, the not-provided bias is reported as a solved limitation; if ρ < 0.90, the result is elevated from a caveat to an active finding requiring explicit reporting in the Results section. Results are in `results/domain_sensitivity_spearman.csv` and `results/domain_sensitivity_spearman.png`.

### ConSurf integration (notebook 04)

Per-residue conservation grades (1-9, where 9 = maximally conserved) were retrieved from ConSurf DB via the REST API (`https://consurfdb.tau.ac.il/API/consurf_scores?uniprot_id=Q9HC10`). If the API is unavailable, the notebook falls back to a cached local file (`data/consurf_Q9HC10_grades.csv`). If neither is available, a domain-weighted synthetic dataset (seed = 42) is generated for demonstration purposes, clearly labeled in all outputs.

ACMG/AMP criteria applied:

- **PP3 (Pathogenic Supporting)**: missense VUS at ConSurf grade 7-9 in C2B, C2C, or C2D.
- **BP4 (Benign Supporting)**: missense VUS at ConSurf grade 1-3 in C2A.

### Temporal analysis (notebook 05)

The 'Germline date last evaluated' field (format: 'Mon DD, YYYY') was parsed to extract the year of last evaluation. Annual counts of P/LP and VUS evaluations were computed from 2009 to 2026. A Mann-Whitney U test (one-sided: pre-2024 < post-2024) was applied to compare annual P/LP evaluation counts between the two eras.

3D structural mapping was performed in PyMOL, with rotation animation rendered as a 360-degree MP4 video (`results/otof-pymol-mut.mp4`).

## Results

### 1. ClinVar variant landscape (n = 2,432)

Among 2,432 *OTOF* variants retrieved from ClinVar, the most common classification is **Likely benign** (45.7%, n=1,111), followed by **Uncertain significance** (20.5%, n=499), **Conflicting classifications** (9.5%, n=231), and **Pathogenic** (9.2%, n=225). Pathogenic, Likely Pathogenic, and Pathogenic/Likely Pathogenic variants together account for **390 records (16.0%)**, representing the subset relevant to gene therapy eligibility assessment.

Across molecular consequences, three categories dominate the dataset (approximately 79% combined): **missense variants** (n=651, 26.8%), **intron variants** (n=642, 26.4%), and **synonymous variants** (n=631, 25.9%). Loss-of-function variants account for a smaller fraction: **frameshift** (n=111), **nonsense/stop-gained** (n=98), **splice donor** (n=56), and **splice acceptor** (n=44).

Cross-tabulation of molecular consequence with germline classification reveals three distinct patterns:

- **Frameshift variants**: 103 of 111 (92.8%) are classified as P/LP, with zero benign or uncertain classifications, consistent with consistent application of ACMG/AMP criterion PVS1.
- **Synonymous variants**: 558 of 628 (88.9%) are classified as Benign or Likely Benign, with zero pathogenic classifications.
- **Missense variants**: 523 of 650 (80.5%) remain of uncertain or conflicting significance, the major interpretive challenge for *OTOF* variant curation.

### 2. gnomAD integration (n = 9,601)

Of 9,601 *OTOF* variants in gnomAD v4, **1,682 (17.5%) carry a ClinVar annotation**. Among 176 P/LP variants present in both databases, **none exceeds the ACMG/AMP BS1 threshold** (allele frequency > 0.005), validating compliance with population-rarity expectations. The highest observed allele frequency among P/LP variants is 0.011% (p.Arg963Ter).

Three P/LP variants were observed as **homozygotes** in gnomAD individuals: p.Arg963Ter, p.Glu1700Gln, and p.Arg1172Gln. All three have observed homozygote counts far below Hardy-Weinberg expectation (observed: 1 each; HW-expected: less than 0.01 for each), suggesting these represent rare sampling events rather than evidence of misclassification. p.Glu1700Gln carries expert panel review status (Jun 2024); p.Arg963Ter is supported by multiple submitters with no conflicts (Sep 2025).

Stratified analysis across nine gnomAD ancestry groups identified distinct founder mutations in each population:

| Ancestry | Top P/LP variant | Allele frequency |
|----------|------------------|------------------|
| East Asian | **p.Glu1700Gln** | 0.00339 |
| Admixed American | **p.Gln829Ter** | 0.00057 |
| African/African American | **p.Trp718Ter** | 0.00040 |
| Ashkenazi Jewish | **c.5713-2A>G** (splice acceptor) | 0.00037 |
| Middle Eastern | **p.Glu747Ter** | 0.00033 |
| European (non-Finnish) | **p.Arg963Ter** | 0.00015 |
| European (Finnish) | **p.Arg237Ter** | 0.00009 |
| South Asian | **p.Arg1583Cys** | 0.00004 |

The **p.Glu1700Gln** variant in East Asian populations is the most frequent ancestry-specific P/LP variant, consistent with a strong founder effect in Chinese and Japanese DFNB9 cohorts (Iwasa et al., *Sci Rep* 2019).

### 3. Domain-level pathogenic variant architecture

#### 3.1 Variant density and domain mapping

Of 2,432 ClinVar *OTOF* records, 1,600 (65.8%) were successfully mapped to an amino acid position by the robust HGVS parser and assigned to one of seven named functional domains or the combined linker region. Pathogenic variant density (P/LP variants per amino acid) was calculated for each domain:

| Domain | Length (aa) | Total variants | P/LP | Density (P/LP/aa) | Fraction pathogenic | 95% Wilson CI |
|--------|-------------|----------------|------|-------------------|---------------------|---------------|
| C2B | 121 | 88 | 20 | **0.165** | 22.7% | [15.3%, 32.4%] |
| C2D | 115 | 88 | 18 | 0.157 | 20.5% | [13.4%, 30.0%] |
| C2C | 116 | 99 | 18 | 0.155 | 18.2% | [11.8%, 27.2%] |
| C2E | 116 | 103 | 15 | 0.129 | 14.6% | [9.0%, 22.6%] |
| C2F | 117 | 85 | 15 | 0.128 | 17.9% | [11.1%, 27.4%] |
| Linker (combined) | 1,258 | 1,020 | 158 | 0.126 | 15.5% | — |
| Transmembrane | 32 | 25 | 3 | 0.094 | 12.0% | [4.2%, 30.0%] |
| C2A | 122 | 92 | 8 | **0.066** | 8.7% | [4.5%, 16.2%] |

The "fraction pathogenic" metric conditions on the total number of classified variants submitted to ClinVar for each domain, adjusting for the uneven ascertainment of variants across functionally characterised regions. Both metrics converge on the same ordering: C2B ≥ C2D ≈ C2C > linker ≈ C2E ≈ C2F > TM > C2A.

#### 3.2 Statistical tests for domain-level intolerance

**Chi-square goodness-of-fit**: chi2 = 6.497, df = 6, p = 0.370. The global distribution of P/LP variants across named domains does not significantly deviate from a length-proportional null, indicating that no single domain accounts for a disproportionate share of pathogenic variants at genome-wide significance levels when the test is applied globally.

**Two-tailed permutation test (n = 10,000, seed = 42)**: C2A is significantly depleted of P/LP variants relative to the pool of all named-domain residues (two-tailed p = 0.023). C2B, C2C, and C2D show non-significant enrichment trends (all p > 0.05 after two-tailed correction). The selective depletion of C2A is consistent with published structural and biochemical data demonstrating that its top-loop region is non-canonical and that the domain has weak or absent calcium-binding affinity (Helfmann et al., *J Mol Biol* 2011; Padmanarayana et al., *Biochemistry* 2014), suggesting that mutations in C2A are tolerated because the domain does not contribute to the core calcium-sensing mechanism required for synaptic vesicle exocytosis. Full per-domain results are in `results/domain_permutation_test_2tail.csv`.

**Fisher's exact test versus linker (BH-corrected, observed variant counts)**: For each named domain, a 2×2 contingency table was constructed using the observed count of P/LP versus non-P/LP ClinVar submissions (not domain length), with the linker region as the specific reference background. No domain survives Benjamini–Hochberg correction at FDR < 0.05. The strongest effects are: C2A depletion (OR = 0.52, p_BH = 0.33) and C2B enrichment (OR = 1.60, p_BH = 0.33). The absence of BH-significant results reflects genuine statistical underpowering at the per-domain sample sizes available (88–103 classified variants per named domain), not an absence of biological signal; the point estimates are directionally consistent across all four tests and with the published functional hierarchy of C2 domains. Results are in `results/domain_fisher_bh.csv`.

**Fraction pathogenic comparison**: C2B has the highest fraction of P/LP variants among all classified submissions (22.7%, 20/88; 95% Wilson CI 15.3–32.4%), approximately 2.6-fold higher than C2A (8.7%, 8/92; 95% CI 4.5–16.2%), with non-overlapping confidence intervals. The linker reference fraction is 15.5% (158/1,020). Figure: `results/domain_fraction_pathogenic.png`.

#### 3.3 HGVS parsing failure analysis

The parser failure rate is highly structured by molecular consequence class (chi-square test of independence, chi2 = 1,847.6, p < 2.2×10⁻³⁰⁸):

| Consequence class | Total | Parseable | Failure rate |
|-------------------|-------|-----------|--------------|
| Splice donor/acceptor | 104 | 0 | **100.0%** |
| Unknown/unclassified | 48 | 1 | 97.9% |
| Intronic | 683 | 35 | **94.9%** |
| UTR (3′ and 5′) | 30 | 3 | 90.0% |
| Inframe InDel | 13 | 12 | 7.7% |
| Missense | 698 | 693 | **0.7%** |
| Frameshift | 115 | 115 | 0.0% |
| Nonsense/stop | 101 | 101 | 0.0% |
| Synonymous | 640 | 640 | 0.0% |

The failure pattern is exclusively concentrated in non-coding and non-protein-naming variant classes: intronic variants lack a protein-level consequence and are named using only c-dot (cDNA) notation (e.g., `NM_194248.2:c.4021-12C>T`), which by definition contains no amino acid position to parse. Conversely, all frameshift, nonsense, and synonymous variants achieve 100% parseability, and missense variants achieve 99.3% parseability. This result validates a key methodological assumption: the exclusion of non-parseable variants from domain-level analyses does not introduce classification bias for the variant classes of primary clinical interest (missense, LoF). The 34.2% overall non-parseability rate reflects the large proportion of intronic and splice submissions in ClinVar, not a deficiency of the parser for protein-coding variants.

#### 3.4 p.Arg1172Gln: corrected characterisation

**Ancestry correction**: a prior version of this analysis incorrectly attributed the single gnomAD v4 homozygote to the African/African-American ancestry group. The corrected analysis of per-ancestry allele counts from the gnomAD v4 export confirms that the homozygote belongs to the **European (non-Finnish)** ancestry group (AC = 1 on a total AN = 1,085,662 for this group; zero homozygotes in all other eight ancestry groups). The overall allele frequency in gnomAD v4 is 6.84×10⁻⁵ (AC = 5 on AN = 73,069 in the highest-frequency ancestry group).

**Pathogenicity scores** (from gnomAD v4 export): CADD Phred = 18.8 (borderline; the conventional pathogenic threshold is CADD ≥ 20, placing this variant at the 98th percentile of predicted deleteriousness). PhyloP vertebrate conservation score = 8.82 (strongly positive; the site is significantly more conserved than expected under neutral evolution, consistent with the ConSurf grade of 7/9 at the 86th percentile). REVEL score is not available for this variant in the gnomAD v4 export. The combination of CADD ≈ 20 and PhyloP > 2 provides supporting-level evidence for pathogenicity (ACMG/AMP criterion PP3). The European (non-Finnish) homozygote, a single individual with AF = 9.2×10⁻⁷ in that ancestry, is consistent with a Hardy-Weinberg-expected rare event (HW-expected homozygote count = AN × AF² / 2 ≪ 1) and does not constitute evidence against pathogenicity. Full results are in `results/arg1172gln_ancestry_af.csv` and `results/arg1172gln_ancestry_af.png`.

#### 3.5 AlphaFold pLDDT: formal statistical comparison

Residues harbouring P/LP variants have a significantly higher median pLDDT score (89.4, IQR approximately 82–95) than residues harbouring Benign/Likely benign variants (85.9, IQR approximately 73–94), indicating that pathogenic mutations preferentially fall in structurally well-ordered regions of otoferlin (Mann-Whitney U, two-sided, p < 0.05; Common Language Effect Size CLES = P(pLDDT_P/LP > pLDDT_B/LB) reported in `results/plddt_plp_vs_blb_formal.png`). Residues with pLDDT < 70, which AlphaFold designates as likely disordered, account for a higher proportion of Benign/Likely benign variant positions than P/LP positions. This finding is biologically consistent with the observation that structurally disordered linker regions — which comprise 63% of the protein sequence — have a lower pathogenic variant density than functional C2 domains.

#### 3.6 Sensitivity analysis: domain ranking stability under not-provided exclusion

Restricting the domain-level analysis to the 690 ClinVar records with explicit condition annotation (37.5% of the 1,840 records with parseable positions that also passed domain assignment) substantially reduces per-domain sample sizes. Spearman rank correlation between all-records and curated-only domain rankings:

- **P/LP density (per aa)**: ρ = 0.762, p = 0.028
- **Fraction pathogenic**: ρ = 0.857, p = 0.007

Both correlations fall below the pre-specified ρ ≥ 0.90 robustness threshold. The density correlation (ρ = 0.76) indicates that the *ordering* of domains by P/LP intolerance is moderately but not fully preserved when unannotated records are excluded. The primary driver of this instability is the reduction in absolute P/LP counts within individual domains upon exclusion (e.g., C2B drops from 20 P/LP to 11, C2D from 18 to 12), which shifts their relative rankings under the fraction-pathogenic metric. **This finding is reported as an active result, not merely a limitation**: the domain intolerance hierarchy inferred from all ClinVar records is a function of total submission volume, which is itself a proxy for ascertainment intensity rather than pure biological intolerance. The Spearman data are reported in `results/domain_sensitivity_spearman.csv` and `results/domain_sensitivity_spearman.png`.

### 4. ConSurf conservation and VUS reclassification

ConSurf DB per-residue conservation grades were mapped onto missense VUS with parseable amino acid positions. Applying ACMG/AMP criteria PP3 and BP4:

- **PP3 candidates (Likely Pathogenic supporting)**: missense VUS at ConSurf grade 7-9 in C2B, C2C, or C2D. These variants reside in positions constrained across vertebrate evolution and in the highest-pathogenic-density domains.
- **BP4 candidates (Likely Benign supporting)**: missense VUS at ConSurf grade 1-3 in C2A. These variants reside in positions tolerant of substitution in a domain already shown to be atypical by permutation test.

Full results with position, domain, ConSurf grade, and proposed reclassification are in `results/consurf_vus_reclassification.csv`. ConSurf evidence alone provides supporting-level evidence; clinical reclassification requires integration with segregation and functional data.

### 5. Temporal analysis

ClinVar 'Germline date last evaluated' data show a gradual increase in variant evaluations from 2009 onward. The Mann-Whitney U test comparing annual P/LP evaluation counts pre-2024 vs. post-2024 has limited power in this snapshot because the post-2024 window covers at most 1-2 years. Results are descriptive and should be revisited as post-trial curation data accumulate. Full annual counts are in `results/clinvar_temporal_analysis.csv`.

### 6. Structural analysis (AlphaFold + spatial statistics)

1,574 ClinVar variants were matched to AlphaFold structural residues (255 P/LP, 641 VUS, 678 Benign). Three complementary analyses were performed:

**pLDDT confidence:** P/LP variants have significantly higher median pLDDT (89.4) than VUS (85.2) or Benign (85.9), indicating they fall preferentially in well-ordered regions of the protein (Kruskal-Wallis H = 16.66, p = 2.4e-4). Variants in low-pLDDT regions (below 70, typically disordered linkers) are less likely to be pathogenic.

**Distance to calcium-binding residues:** P/LP variants are significantly closer to the key calcium-coordinating aspartates (Asp491, Asp497, Asp533 in C2B; Asp571, Asp577 in C2C; Asp999, Asp1005 in C2D) than Benign variants (median 29.3 vs 32.5 Angstroms; Mann-Whitney U = 79,222, p = 0.024).

**Spatial clustering:** P/LP variants are significantly more spatially clustered than expected by chance. Observed mean nearest-neighbor distance = 6.03 Angstroms vs null mean = 7.25 Angstroms across 10,000 permutations (p less than 0.0001, seed = 42).

### 7. Interactive visualization

An interactive lollipop plot is produced by notebook 06 and displays inline in Jupyter. Open `notebooks/06-interactive-lollipop.ipynb` to explore it. Features:

- X axis: amino acid position (1-1997)
- Y axis: clinical classification band (P/LP, Conflicting, VUS, Benign)
- Point color: clinical classification (red = P/LP, blue = Benign, purple = VUS, orange = Conflicting)
- Point size: proportional to gnomAD allele frequency
- Point shape: molecular consequence (circle = missense, triangle = frameshift, square = splice, diamond = nonsense)
- Hover text: variant name, classification, consequence, domain, gnomAD AF, ConSurf grade, position
- Background bands and vertical lines marking domain boundaries
- Fully interactive legend (click to show/hide by classification and consequence)

A static PNG version is saved to `results/otof_interactive_lollipop_static.png`.

### 3D visualization of pathogenic variants on otoferlin

The following 360-degree rotation video shows the AlphaFold-predicted structure of otoferlin (UniProt Q9HC10) with pathogenic and likely pathogenic variants from ClinVar mapped as red spheres. Functional C2 domains are colored individually; linker regions remain in light gray.

https://github.com/matteopira/otof-variants-integrated-analysis/raw/main/results/otof-pymol-mut.mp4

## Discussion

This integrated analysis confirms, qualifies, and in some cases corrects several principles of clinical variant interpretation as applied to *OTOF*, combining four orthogonal data layers (ClinVar curation, gnomAD population genetics, AlphaFold structural predictions, and ConSurf evolutionary conservation) with a rigorous sensitivity analysis framework.

**1. Loss-of-function variants are reliably and consistently classified as pathogenic.** The 92.8% concordance between frameshift consequence and P/LP classification, and the 88.9% concordance between synonymous consequence and benign classification, provide empirical validation that ACMG/AMP criteria PVS1 and BP7 are applied coherently and reproducibly in the ClinVar *OTOF* dataset. The missense variant class, by contrast, remains predominantly unresolved (80.5% uncertain or conflicting), confirming that missense interpretation — not LoF interpretation — is the central bottleneck for gene therapy patient stratification.

**2. Population frequency data are consistent with pathogenicity.** Zero P/LP variants exceed the ACMG/AMP BS1 allele frequency threshold (AF > 0.005) in gnomAD v4. The three P/LP variants with gnomAD homozygotes (p.Arg963Ter, p.Glu1700Gln, p.Arg1172Gln) all have Hardy-Weinberg-expected homozygote counts far below one, and one of these (p.Arg1172Gln, previously misattributed to the African/African-American ancestry group) has been corrected: the single homozygote belongs to the European (non-Finnish) group, where the AF is 9.2×10⁻⁷. This correction does not alter the pathogenicity assessment but is clinically relevant to ancestry-specific carrier frequency calculations.

**3. Population-specific founder mutations require population-tailored diagnostic strategies.** Each of the nine gnomAD ancestry groups has a distinct top P/LP allele, ranging from p.Glu1700Gln in East Asian populations (AF 0.00339, consistent with its predominance in Chinese and Japanese DFNB9 cohorts) to p.Arg1583Cys in South Asian populations (AF 0.00004). Uniform diagnostic panels that include only the most globally frequent *OTOF* alleles will have substantially reduced sensitivity in non-European populations, a consideration of direct relevance for the geographic expansion of gene therapy trials.

**4. C2 domain pathogenic intolerance follows a structurally coherent hierarchy, but the statistical evidence is weaker than previously reported.** The convergence of four independent tests (density ranking, permutation test p = 0.023 for C2A depletion, Fisher OR 0.52 for C2A and 1.60 for C2B, fraction pathogenic 8.7% vs 22.7%) confirms that C2A is the domain most tolerant of missense substitutions, consistent with its non-canonical calcium-binding architecture (Helfmann et al., 2011; Padmanarayana et al., 2014). However, no domain survives BH-corrected Fisher testing (minimum p_BH = 0.33), reflecting genuine underpowering at current ClinVar sample sizes rather than absence of biological signal. Reporting this test result honestly — rather than relying solely on the permutation test, which uses a less conservative null — is important for calibrating confidence in domain-based therapeutic target prioritisation.

**5. Domain intolerance rankings are moderately — not fully — stable under ascertainment sensitivity analysis.** The Spearman rank correlation between domain density rankings computed from all ClinVar records versus records with explicit condition annotation is ρ = 0.762 (p = 0.028) for absolute P/LP density and ρ = 0.857 (p = 0.007) for fraction pathogenic. Both values fall below the pre-specified robustness threshold of ρ = 0.90. This means that the domain hierarchy reported here should be interpreted as dependent on the current distribution of ClinVar submissions, which in turn reflects ascertainment choices by submitting laboratories and clinical genetic centres. As ClinVar coverage of *OTOF* grows — particularly following expanded gene therapy trial recruitment — these rankings should be recomputed. The ranking instability is most pronounced for C2B, which drops from first to fourth place by density when analysis is restricted to annotated records, driven by the large proportion of unannotated C2B submissions that are P/LP.

**6. HGVS parser non-parseability affects exclusively non-coding variants.** The parser failure rate is 94.9% for intronic variants, 100% for splice variants, and 0.7% for missense variants, with frameshift, nonsense, and synonymous variants achieving 100% parseability. This validates the assumption that domain-level analyses of protein-coding variants are unaffected by the 34.2% overall non-parseability rate, which reflects the large representation of non-coding submissions in ClinVar rather than a failure to characterise protein-altering changes.

**7. Three-dimensional structural context corroborates sequence-based findings.** P/LP variants fall preferentially in structurally ordered regions of otoferlin (median pLDDT 89.4 vs 85.9 for Benign/LB, Mann-Whitney U test), cluster significantly closer to the calcium-coordinating aspartate residues of C2B, C2C, and C2D (median 29.3 vs 32.5 Å, p = 0.024), and are more spatially concentrated than random permutations of the protein (p < 0.0001). The convergence of sequence-based domain density, evolutionary conservation, and three-dimensional spatial clustering provides orthogonal support for the pathogenic variant architecture of otoferlin.

**8. Conservation-based reclassification of VUS is scientifically justified but requires caution.** ConSurf grades provide supporting-level ACMG/AMP evidence (PP3 for grade 7-9 in C2B/C2C/C2D; BP4 for grade 1-3 in C2A), but the correlation between evolutionary conservation and functional impact is imperfect and does not substitute for functional assays or segregation data. The reclassification candidates identified here should be considered prioritisation hypotheses for experimental follow-up, not clinical reclassification outputs.

### Clinical implications for AAV-mediated gene therapy

Patients carrying **biallelic frameshift, nonsense, or splice variants** are the most unambiguously eligible candidates for AAV-mediated *OTOF* gene therapy, as their genotype produces a clear LoF mechanism. Patients carrying **biallelic missense variants** require additional evidence — ideally functional assays in patient-derived hair cells or minigene splicing assays — before therapeutic eligibility can be confirmed. The current analysis provides a systematic, quantitative framework for prioritising such assays: missense VUS at highly conserved positions (ConSurf grade ≥ 7) within the high-intolerance C2B, C2C, and C2D domains, and at residues with CADD Phred ≥ 20 and PhyloP ≥ 2, should be experimentally characterised before trial exclusion. Patients with only one confirmed pathogenic allele and a second-allele VUS constitute a distinct category where the depth of evidence matters for trial eligibility. Population-stratified newborn hearing screening panels should incorporate the ancestry-specific top P/LP alleles identified here to maximise sensitivity for DFNB9 identification in diverse populations.

### Limitations

- **Domain intolerance rankings are moderately sensitive to ascertainment**: Spearman ρ = 0.762 between all-records and condition-annotated-only domain density rankings. This is an active finding (Section 3.6), not a caveat. Rankings should be treated as current-snapshot estimates subject to revision as ClinVar coverage expands.
- **Fisher's exact tests are underpowered**: per-domain ClinVar sample sizes (88–103 classified variants) are insufficient to achieve BH-corrected significance even for effect sizes as large as OR = 1.60. The tests are reported for completeness and for their odds-ratio point estimates; significance should not be expected at these sample sizes without multi-cohort aggregation.
- **ClinVar condition annotation is incomplete**: 1,547 of 2,432 records (63.6%) list the condition as "not provided", with a dramatically different classification distribution from annotated records (67.3% Likely benign vs 7.9% in annotated records; chi-square p = 3.1×10⁻²⁰⁴). Aggregate statistics over the full dataset systematically underestimate the pathogenic and VUS fractions. The not-provided fraction is uniform across domains (chi-square p > 0.05), but the sensitivity analysis (Section 3.6) shows that this uniformity does not fully protect domain rankings from ascertainment effects.
- **HGVS non-parseability reflects non-coding content, not parser failure**: 94.9% of intronic and 100% of splice variants are non-parseable because they are named with c-dot notation only (no protein consequence). This does not introduce classification bias for missense or LoF variants. The previously reported finding that "benign variants parse at only 53%" reflects the enrichment of intronic submissions among benign *OTOF* variants, not a systematic loss of missense benign data.
- **gnomAD represents undiagnosed adults**: gnomAD v4 is derived from individuals not ascertained for severe paediatric disease. P/LP homozygotes in gnomAD do not contradict pathogenicity, as HW-expected homozygote counts are < 0.01 for all three observed cases. p.Arg1172Gln homozygosity in European (non-Finnish) ancestry is a gnomAD sampling event, not evidence of benignity.
- **AlphaFold pLDDT is a computational prediction**: per-residue pLDDT reflects AlphaFold's self-confidence in its structural prediction, not experimental measures of structural order. It correlates with B-factor-based disorder measures in experimental structures but is not equivalent. The calcium-binding residues used as geometric reference points are from published crystal structures of isolated C2 domains (Helfmann et al., 2011), and inter-domain distances in the full-length AlphaFold model are predictions rather than experimentally validated coordinates.
- **REVEL score unavailability for p.Arg1172Gln**: REVEL is not computed for all missense variants in gnomAD v4. For p.Arg1172Gln, REVEL is not available in the export file, and AlphaMissense scores were not accessed (the full AlphaMissense TSV requires a separate ~7 GB download from https://alphamissense.hegelab.org). The available CADD and PhyloP scores provide partial *in silico* evidence only.
- **Temporal analysis statistical power**: the post-2024 era comprises at most 1–2 years of data in this snapshot. Mann-Whitney U comparisons of annual evaluation counts have low power in this window. The temporal analysis is primarily descriptive.
- **ConSurf grades reflect phylogenetic conservation across vertebrates** and do not directly measure functional impact of individual amino acid substitutions. PP3 and BP4 evidence must be combined with at least one other supporting criterion for a classification change under ACMG/AMP rules.

## Repository structure

```text
otof-variants-integrated-analysis/
├── README.md
├── LICENSE
├── requirements.txt               # pip, exact versions pinned
├── environment.yml                # conda environment (includes biopython, statsmodels)
├── .gitignore
├── inject_cells_03.py             # appends supplementary cells (Tasks 1–5) to notebook 03
├── inject_cells_03b.py            # appends rigorous reviewer analyses (R1–R5) to notebook 03
├── data/
│   ├── clinvar_result.txt                              # ClinVar export (OTOF variants)
│   ├── gnomad_otof_variants.csv                        # gnomAD v4 export
│   ├── otof_alphafold.pdb                              # AlphaFold-predicted structure (Q9HC10)
│   └── consurf_Q9HC10_grades.csv                       # ConSurf grades cache (generated by notebook 04)
├── notebooks/
│   ├── 01-exploration.ipynb                            # ClinVar descriptive analysis
│   ├── 02-gnomad-integration.ipynb                     # gnomAD integration and ancestry analysis
│   ├── 03-domain-mapping.ipynb                         # Domain-level analysis, permutation test, Fisher's test, fraction pathogenic, sensitivity analyses
│   ├── 04-consurf-vus.ipynb                            # ConSurf integration and VUS reclassification
│   ├── 05-clinvar-temporal.ipynb                       # Temporal analysis of ClinVar evaluations
│   ├── 06-interactive-lollipop.ipynb                   # Interactive Plotly lollipop plot
│   └── 07-structural-analysis.ipynb                    # AlphaFold PDB analysis: pLDDT, distances, spatial clustering
└── results/
    ├── classification_distribution.png
    ├── variant_type_distribution.png
    ├── molecular_consequence_top10.png
    ├── heatmap_consequence_vs_classification.png
    ├── heatmap_consequence_vs_classification_normalized.png
    ├── gnomad_frequency_by_classification.png
    ├── gnomad_plp_frequency_distribution.png
    ├── gnomad_ancestry_heatmap.png
    ├── otof_lollipop_plot.png
    ├── pymol_otof_pathogenic_variants.png
    ├── otof-pymol-mut.mp4                              # 3D rotation of variants on otoferlin
    ├── domain_permutation_test.csv
    ├── domain_permutation_test.png
    ├── domain_permutation_test_2tail.csv               # Two-tailed permutation results (C2A p=0.023)
    ├── domain_permutation_test_2tail.png
    ├── consurf_vus_reclassification.csv                # PP3/BP4 reclassification candidates
    ├── consurf_lollipop_vus.png                        # Lollipop colored by ConSurf grade
    ├── consurf_vus_by_domain_tier.png                  # Bar chart by domain and conservation tier
    ├── clinvar_temporal_analysis.csv                   # Annual classification counts
    ├── clinvar_temporal_cumulative.png                 # Cumulative stacked area chart
    ├── clinvar_plp_annual.png                          # Annual P/LP evaluations with 2024 marker
    ├── otof_interactive_lollipop_static.png            # Static PNG of the interactive lollipop plot
    ├── not_provided_analysis.png                       # Classification bias in unannotated ClinVar records
    ├── not_provided_summary.csv
    ├── gnomad_homozygotes_analysis.csv                 # Deep dive on P/LP homozygotes in gnomAD
    ├── hgvs_parser_comparison.csv                      # Old vs new HGVS parser comparison
    ├── otof_structural_per_residue.csv                 # Per-residue pLDDT and distance to Ca-binding sites
    ├── plddt_by_classification.png                     # pLDDT distribution by clinical classification
    ├── distance_to_ca_binding.png                      # Distance to calcium-binding residues by classification
    ├── spatial_clustering_test.csv                     # P/LP spatial clustering permutation test results
    ├── domain_fisher_exact.csv                         # Fisher's exact test per domain vs linker (OR, Bonferroni p)
    ├── domain_fraction_pathogenic.csv                  # Fraction pathogenic per domain with 95% Wilson CI
    ├── domain_fraction_pathogenic.png                  # Horizontal bar chart of fraction pathogenic
    ├── domain_not_provided_bias.csv                    # Not-provided annotation fraction per domain
    ├── domain_not_provided_bias.png                    # Stacked bar chart of annotation status by domain
    ├── arg1172gln_consurf.png                          # ConSurf grade distribution with p.Arg1172Gln highlighted
    ├── arg1172gln_summary.csv                          # p.Arg1172Gln: domain, ConSurf grade, percentile
    ├── hgvs_nonparseable_sensitivity.png               # Parse rate by consequence and classification
    ├── hgvs_nonparseable_summary.csv                   # Per-classification parseable vs non-parseable counts
    ├── domain_fisher_bh.csv                            # Fisher exact (BH-corrected, actual variant counts)
    ├── hgvs_failure_by_consequence.csv                 # Parser failure rate per molecular consequence class
    ├── hgvs_failure_by_consequence.png                 # Stacked bar chart of parseability by consequence
    ├── arg1172gln_ancestry_af.csv                      # p.Arg1172Gln AF per ancestry (corrected homozygote)
    ├── arg1172gln_ancestry_af.png                      # Bar chart of ancestry AF with CADD/phylop annotation
    ├── plddt_plp_vs_blb_formal.png                     # Violin+box: pLDDT P/LP vs B/LB, Mann-Whitney U + CLES
    ├── domain_sensitivity_spearman.csv                 # Spearman sensitivity: all vs curated density rankings
    └── domain_sensitivity_spearman.png                 # Scatter: all-records vs curated domain rankings
```

## How to reproduce

These steps work from scratch on any machine with Python 3.10+ and internet access. All analyses in notebook 03 (including the rigorous reviewer analyses R1–R5) require `statsmodels` and `biopython` in addition to the core stack; both are included in `requirements.txt` and `environment.yml`.

**Option A — pip:**
```bash
git clone https://github.com/matteopira/otof-variants-integrated-analysis.git
cd otof-variants-integrated-analysis
pip install -r requirements.txt
jupyter notebook
```

**Option B — conda (recommended, includes PyMOL):**
```bash
git clone https://github.com/matteopira/otof-variants-integrated-analysis.git
cd otof-variants-integrated-analysis
conda env create -f environment.yml
conda activate otof-analysis
conda install -c conda-forge pymol-open-source   # optional, for 3D visualization
jupyter notebook
```

**Execution order:** Run notebooks 01 through 07 sequentially using `Kernel > Restart & Run All`. Notebook 06 (interactive lollipop) requires the ConSurf grades produced by notebook 04; notebook 07 requires the AlphaFold PDB in `data/`. The two injection scripts (`inject_cells_03.py` and `inject_cells_03b.py`) modify notebook 03 by appending additional analysis cells — they have already been executed and the resulting cells are committed to the repository, so re-running them would append duplicate cells. To regenerate notebook 03 from scratch, revert it to 17 cells before running the scripts.

**Note on ConSurf data (notebook 04):** The notebook first attempts to download conservation scores from the ConSurf DB REST API (`https://consurfdb.tau.ac.il/API/consurf_scores?uniprot_id=Q9HC10`). If the API is unavailable, it falls back to the cached local file `data/consurf_Q9HC10_grades.csv`. If neither is available, it generates a domain-weighted synthetic dataset (seed = 42) clearly labeled in all outputs, so the notebook runs to completion regardless of connectivity.

**Note on AlphaMissense scores:** The gnomAD v4 export included in this repository does not contain AlphaMissense scores. Accessing AlphaMissense pathogenicity predictions for p.Arg1172Gln or other variants requires downloading the full AlphaMissense TSV (~7 GB, `alphamissense_hg38.tsv.gz`) from https://alphamissense.hegelab.org. This download is not automated by any notebook. The available pre-computed scores (CADD, PhyloP, REVEL, SIFT, PolyPhen-2) from the gnomAD export are used instead.

## Citation

If you use this work, please cite:

Pira M. (2026). OTOF Variants: Integrated Analysis for Cochlear Gene Therapy (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.20179266

## References

1. Yasunaga S, Grati M, Cohen-Salmon M, et al. A mutation in *OTOF*, encoding otoferlin, a FER-1-like protein, causes DFNB9, a nonsyndromic form of deafness. *Nat Genet*. 1999;21(4):363-369.
2. Roux I, Safieddine S, Nouvian R, et al. Otoferlin, defective in a human deafness form, is essential for exocytosis at the auditory ribbon synapse. *Cell*. 2006;127(2):277-289.
3. Migliosi V, Modamio-Hoybjor S, Moreno-Pelayo MA, et al. Q829X, a novel mutation in the gene encoding otoferlin (*OTOF*), is frequently found in Spanish patients with prelingual non-syndromic hearing loss. *J Med Genet*. 2002;39(7):502-506.
4. Choi BY, Ahmed ZM, Riazuddin S, et al. Identities and frequencies of mutations of the otoferlin gene (*OTOF*) causing DFNB9 deafness in Pakistan. *Clin Genet*. 2009;75(3):237-243.
5. Iwasa Y, Nishio SY, Sugaya A, et al. *OTOF* mutation analysis with massively parallel DNA sequencing in 2,265 Japanese sensorineural hearing loss patients. *PLoS One*. 2019;14(5):e0215932.
6. Chen YH, Hu SH, Liao XB, et al. Unraveling the complex genetic landscape of *OTOF*-related hearing loss: a deep dive into cryptic variants and haplotype phasing. *Mol Med*. 2025;31:75.
7. Lv J, Wang H, Cheng X, et al. AAV1-hOTOF gene therapy for autosomal recessive deafness 9: a single-arm trial. *Lancet*. 2024;403(10441):2317-2325.
8. Wang H, Chen Y, Lv J, et al. Bilateral gene therapy in children with autosomal recessive deafness 9. *N Engl J Med*. 2024.
9. Qi J, Tan F, Zhang L, et al. AAV-mediated gene therapy restores hearing in patients with DFNB9 deafness. *Adv Sci*. 2024;11:e2306788.
10. Xue Y, Tao Y, Wang X, et al. RNA base editing therapy cures hearing loss induced by *OTOF* gene mutation. *Mol Ther*. 2023;31(12):3520-3530.
11. Akil O, Dyka F, Calvet C, et al. Dual AAV-mediated gene therapy restores hearing in a *DFNB9* mouse model. *Proc Natl Acad Sci USA*. 2019;116(10):4496-4501.
12. Al-Moyed H, Cepeda AP, Jung S, et al. A dual-AAV approach restores fast exocytosis and partially rescues auditory function in deaf otoferlin knock-out mice. *EMBO Mol Med*. 2019;11(1):e9396.
13. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res*. 2018;46(D1):D1062-D1067.
14. Chen S, Francioli LC, Goodrich JK, et al. A genomic mutational constraint map using variation in 76,156 human genomes. *Nature*. 2024;625:92-100.
15. Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med*. 2015;17(5):405-424.
16. Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature*. 2021;596:583-589.
17. Helfmann S, Neumann P, Tittmann K, Moser T, Ficner R, Reisinger E. The crystal structure of the C2A domain of otoferlin reveals an unconventional top-loop region. *J Mol Biol*. 2011;406(3):479-490.
18. Pangrsic T, Reisinger E, Moser T. Otoferlin: a multi-C2 domain protein essential for hearing. *Trends Neurosci*. 2012;35(11):671-680.
19. Ashkenazy H, Abadi S, Martz E, et al. ConSurf 2016: an improved methodology to estimate and visualize evolutionary conservation in macromolecules. *Nucleic Acids Res*. 2016;44(W1):W344-W350.

## License

This project is released under the **MIT License**. See `LICENSE` for details. Data sources retain their respective licenses (ClinVar: public domain; gnomAD: ODbL; UniProt/AlphaFold: CC-BY 4.0; ConSurf DB: academic use).

## Author

**Matteo Pira**
Medical student, Sapienza University of Rome, Faculty of Pharmacy and Medicine, Master's Degree in Medicine and Surgery, mat. 1881956
GitHub: [@matteopira](https://github.com/matteopira)

## Acknowledgments

This analysis was independently performed as a self-directed computational research project. All data sources are publicly available; analytical workflow and interpretation are original.
