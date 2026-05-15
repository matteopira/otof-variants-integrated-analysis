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
- **`notebooks/03-domain-mapping.ipynb`**: structural mapping, including robust HGVS protein notation parsing with multi-format support and parser comparison logging, variant assignment to the six C2 domains and the transmembrane domain (UniProt Q9HC10), pathogenic variant density calculation per domain, lollipop plot visualization, and chi-square and permutation-based statistical testing.
- **`notebooks/04-consurf-vus.ipynb`**: ConSurf DB integration, per-residue conservation grade retrieval, conservation-based classification of missense VUS, and ACMG/AMP PP3/BP4 reclassification candidates.
- **`notebooks/05-clinvar-temporal.ipynb`**: temporal analysis of ClinVar variant evaluations by year, cumulative classification trends, and comparison of pre-2024 vs. post-2024 P/LP evaluation rates using Mann-Whitney U test.
- **`notebooks/06-interactive-lollipop.ipynb`**: interactive Plotly lollipop figure integrating all data layers, displayed inline in Jupyter with static PNG export.
- **`notebooks/07-structural-analysis.ipynb`**: AlphaFold structure-based analysis including manual PDB parsing, pLDDT confidence score extraction, per-variant distance calculation to calcium-coordinating residues, pLDDT distribution comparison across classification groups (Kruskal-Wallis test), distance-to-calcium-binding-residues comparison (Mann-Whitney U test), and spatial clustering analysis of P/LP variants via permutation test (n=10,000, seed=42).

### Statistical testing (notebook 03)

Two statistical tests were applied to assess whether the distribution of P/LP variants across C2 domains deviates from a length-proportional null:

1. **Chi-square goodness-of-fit**: tests whether the global distribution across all domains differs from length-proportional expectations.
2. **Permutation test (n = 10,000, two-tailed, seed = 42)**: for each domain, tests both enrichment and depletion by randomly reassigning P/LP variants across named-domain residues. The two-tailed p-value is `2 * min(p_enrichment, p_depletion)`, capped at 1.0.

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

### 3. Structural mapping and statistical testing

Variants with extractable amino acid positions (n = 1,600) were mapped onto the six C2 domains and the transmembrane domain (UniProt Q9HC10, 1997 aa). Pathogenic variant density (P/LP per residue) was calculated for each domain:

| Domain | Length (aa) | P/LP variants | Density (P/LP per aa) |
|--------|-------------|---------------|------------------------|
| C2B | 121 | 20 | **0.165** (highest) |
| C2D | 115 | 18 | 0.157 |
| C2C | 116 | 18 | 0.155 |
| C2E | 116 | 15 | 0.129 |
| C2F | 117 | 15 | 0.128 |
| Linker (combined) | 1258 | 158 | 0.126 (baseline) |
| Transmembrane | 32 | 3 | 0.094 |
| C2A | 122 | 8 | **0.066** (lowest) |

**Statistical results:**

- **Chi-square goodness-of-fit**: chi2 = 6.497, df = 6, p = 0.370 (not significant globally).
- **Permutation test (two-tailed, n = 10,000, seed = 42)**: C2A is significantly depleted (p = 0.023). C2B, C2C, and C2D show a trend toward enrichment but do not reach significance individually. Full results in `results/domain_permutation_test_2tail.csv`.

The significant depletion of C2A is consistent with published data describing it as atypical with weak or absent calcium-binding affinity (Helfmann et al., *J Mol Biol* 2011; Padmanarayana et al., *Biochemistry* 2014).

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

This integrated analysis confirms several principles of clinical variant interpretation as applied to *OTOF* and refines them with quantitative data:

1. **Loss-of-function variants are reliably classified as pathogenic**, supporting the consistent application of ACMG/AMP PVS1 for *OTOF*. The 92.8% concordance between frameshift consequence and pathogenic classification provides empirical validation of current curation practices.

2. **Population frequency data align with clinical classification**: zero P/LP variants exceed the BS1 threshold in gnomAD, demonstrating that the *OTOF* pathogenic variant pool in ClinVar is consistent with population rarity expectations for an autosomal recessive condition.

3. **Population-specific founder mutations dominate the pathogenic allele pool**: each ancestry group has distinct most-common pathogenic alleles, with strong implications for population-tailored diagnostic panels and gene therapy patient recruitment strategies.

4. **C2 functional domains are mutation-intolerant**, with C2B-D showing the highest pathogenic variant density. C2A, by contrast, is significantly depleted of pathogenic variants (permutation p = 0.023), likely reflecting its non-canonical biology rather than under-curation.

5. **Conservation-based evidence can prioritize VUS for reclassification**, but PP3 and BP4 provide only supporting-level evidence and require integration with functional and segregation data.

6. **Three-dimensional structure corroborates sequence-based findings**: P/LP variants cluster significantly closer to calcium-coordinating residues (p = 0.024) and are more spatially concentrated than random (p less than 0.0001), confirming that functional constraint operates at the structural level and not only at the domain-label level.

7. **ClinVar condition annotation status introduces a strong classification bias**: records without condition annotation (63.6% of the dataset) are disproportionately Likely benign (67.3% vs 7.9% in annotated records; chi-square p = 3.1e-204). Aggregate statistics that do not stratify by annotation status systematically underestimate the VUS fraction.

### Clinical implications for AAV-mediated gene therapy

Patients carrying **biallelic frameshift, nonsense, or splice variants** represent the most unambiguously eligible candidates for gene therapy. Patients carrying **biallelic missense variants** would benefit from additional functional characterization (minigene assays, cellular models) before therapeutic eligibility can be confirmed. Population-specific screening panels for newborn hearing screening should incorporate ancestry-relevant founder mutations.

### Limitations

- **ClinVar condition annotation is incomplete**: 1,547 of 2,432 records (63.6%) list the associated condition as "not provided", with a significantly different classification distribution from annotated records (chi-square p = 3.1e-204). Stratification by annotation status is required for unbiased aggregate statistics.
- **gnomAD represents adult individuals** and explicitly excludes those with diagnosed severe pediatric disease. The three P/LP homozygotes identified are consistent with HW-expected rare sampling events, not misclassification.
- The **AlphaFold-predicted structure** is computational; experimentally determined structures of full-length otoferlin are not yet available. The calcium-binding residues used as reference points in the structural analysis are from published crystallographic data of individual C2 domains (Helfmann et al. 2011), not full-length coordinates.
- **Short-read sequencing has limitations** for *OTOF* variant detection, including difficulties in phasing variants across the large gene and in detecting deep intronic and complex structural variants.
- The **temporal analysis** uses 'date last evaluated', a point-in-time snapshot, not historical submission records. True reclassification history requires the ClinVar variant history VCF.
- **ConSurf scores** reflect evolutionary conservation across vertebrates; they do not directly report functional impact and should be interpreted alongside experimental evidence.

## Repository structure

```text
otof-variants-integrated-analysis/
├── README.md
├── LICENSE
├── requirements.txt               # pip, exact versions pinned
├── environment.yml                # conda environment
├── .gitignore
├── data/
│   ├── clinvar_result.txt                              # ClinVar export (OTOF variants)
│   ├── gnomad_otof_variants.csv                        # gnomAD v4 export
│   ├── otof_alphafold.pdb                              # AlphaFold-predicted structure (Q9HC10)
│   └── consurf_Q9HC10_grades.csv                       # ConSurf grades cache (generated by notebook 04)
├── notebooks/
│   ├── 01-exploration.ipynb                            # ClinVar descriptive analysis
│   ├── 02-gnomad-integration.ipynb                     # gnomAD integration and ancestry analysis
│   ├── 03-domain-mapping.ipynb                         # Domain-level analysis, permutation test, lollipop plot
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
    └── spatial_clustering_test.csv                     # P/LP spatial clustering permutation test results
```

## How to reproduce

These steps work from scratch on any machine with Python 3.10+ and internet access.

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

For each notebook, use `Kernel > Restart & Run All`. Run in sequential order (01 through 07): notebook 06 reads ConSurf data produced by notebook 04, and notebook 07 reads the AlphaFold PDB from `data/`.

**Note on ConSurf data (notebook 04):** The notebook first attempts to download conservation scores from the ConSurf DB REST API. If the API is unavailable, it falls back to a cached local file (`data/consurf_Q9HC10_grades.csv`). If neither is available, it generates a synthetic dataset (seed = 42) clearly labeled in all outputs, so the notebook runs to completion regardless.

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
