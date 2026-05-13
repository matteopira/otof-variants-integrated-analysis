# OTOF Variants: Integrated Analysis for Cochlear Gene Therapy

Independent reanalysis of clinically annotated variants of the *OTOF* gene, the molecular target of first-in-human AAV-mediated cochlear gene therapy trials, through integration of ClinVar clinical curation data, gnomAD population frequencies, and AlphaFold-predicted structural information.

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

These breakthroughs have made the accurate molecular characterization of *OTOF* variants, particularly the distinction between pathogenic and benign variants, and the identification of population-specific founder mutations, a problem of immediate clinical relevance. Patient eligibility for gene therapy requires unambiguous biallelic pathogenic genotyping; missense variants of uncertain significance (VUS), in particular, represent a major bottleneck in patient stratification.

## Research questions

This project addresses four integrated questions:

1. **Clinical curation**: What is the distribution of *OTOF* variants in ClinVar across germline classification categories, and how does it correlate with predicted molecular consequence?
2. **Population genetics**: Do *OTOF* pathogenic variants in ClinVar comply with the ACMG/AMP BS1 criterion (allele frequency < 0.005) when cross-referenced with gnomAD, and are there ancestry-specific founder mutations relevant for population-tailored diagnostic and therapeutic strategies?
3. **Structural genomics**: Are pathogenic variants preferentially localized to specific functional domains of otoferlin, and does the variant density per domain reflect functional intolerance?
4. **Translational implications**: How can the integrated variant landscape inform patient stratification for AAV-mediated cochlear gene therapy?

## Data sources

| Database | Version / Access date | Records retrieved |
|----------|----------------------|-------------------|
| **ClinVar** (NCBI) | Retrieved May 2026 | 2,432 *OTOF* variants |
| **gnomAD** (Broad Institute) | v4 | 9,601 *OTOF* variants |
| **UniProt** | Q9HC10 (canonical, 1997 aa) | Reference protein sequence |
| **AlphaFold DB** (EMBL-EBI) | Model v4 | Predicted full-length 3D structure |

All datasets are publicly available and were retrieved through their respective web interfaces. Data files are included in the `data/` directory.

## Methods

All computational analyses were performed in **Python 3.10** within the Jupyter Notebook environment. Structural visualization was conducted in **PyMOL** (open-source build, conda-forge channel).

### Software stack
- **pandas** ≥ 2.0, for tabular data manipulation
- **NumPy** ≥ 1.24, for numerical operations
- **matplotlib** ≥ 3.7 and **seaborn** ≥ 0.12, for data visualization
- **re** (Python standard library), for regular expression parsing of HGVS protein annotations
- **PyMOL** (open-source), for molecular visualization and animation

### Analytical workflow

The analysis is structured into three Jupyter notebooks, each addressing one analytical dimension:

- **`notebooks/01-exploration.ipynb`**: ClinVar descriptive analysis, including variant counts, germline classification distribution, variant type and molecular consequence stratification, and cross-tabulation of consequence versus classification.
- **`notebooks/02-gnomad-integration.ipynb`**: gnomAD v4 integration, including allele frequency distributions stratified by clinical classification, ACMG/AMP BS1 validation, identification of suspicious variants, and ancestry-specific founder mutation analysis across nine population groups.
- **`notebooks/03-domain-mapping.ipynb`**: structural mapping, including parsing of amino acid positions from HGVS protein annotations using regular expressions, assignment of variants to the six C2 domains and the transmembrane domain (UniProt Q9HC10 annotation), pathogenic variant density calculation per domain, and lollipop plot visualization.

3D structural mapping of pathogenic variants onto the AlphaFold-predicted otoferlin structure was performed in PyMOL, with rotation animation rendered as a 360° MP4 video (`results/otof-pymol-mut.mp4`).

## Results

### 1. ClinVar variant landscape (n = 2,432)

Among 2,432 *OTOF* variants retrieved from ClinVar, the most common classification is **Likely benign** (45.7%, n=1,111), followed by **Uncertain significance** (20.5%, n=499), **Conflicting classifications** (9.5%, n=231), and **Pathogenic** (9.2%, n=225). Pathogenic, Likely Pathogenic, and Pathogenic/Likely Pathogenic variants together account for **390 records (16.0%)**, representing the subset relevant to gene therapy eligibility assessment.

Across molecular consequences, three categories dominate the dataset (approximately 79% combined): **missense variants** (n=651, 26.8%), **intron variants** (n=642, 26.4%), and **synonymous variants** (n=631, 25.9%). Predicted loss-of-function variants account for a smaller fraction: **frameshift** (n=111), **nonsense/stop-gained** (n=98), **splice donor** (n=56), and **splice acceptor** (n=44).

Cross-tabulation of molecular consequence with germline classification reveals three sharply distinct patterns:

- **Frameshift variants**: 103 of 111 (92.8%) are classified as P/LP/P-LP, with **zero benign or uncertain classifications**. This confirms the consistent application of ACMG/AMP criterion **PVS1** (predicted loss-of-function in a gene where loss-of-function is the established disease mechanism) for *OTOF*.
- **Synonymous variants**: 558 of 628 (88.9%) are classified as Benign or Likely Benign, with **zero pathogenic classifications**, consistent with the absence of expected functional impact on protein sequence.
- **Missense variants**: 523 of 650 (80.5%) remain of uncertain or conflicting significance, representing the major interpretive challenge for *OTOF* variant curation and the primary obstacle to patient stratification for gene therapy.

### 2. gnomAD integration (n = 9,601)

Of 9,601 *OTOF* variants observed in gnomAD v4, **1,682 (17.5%) carry a ClinVar annotation**, defining the integrated dataset for population-level analyses. Among 176 P/LP variants present in both databases, **none exceeds the ACMG/AMP BS1 threshold** (allele frequency > 0.005), validating compliance with population-rarity expectations for an autosomal recessive disorder. The highest observed allele frequency among P/LP variants is 0.011% (p.Arg963Ter in the global gnomAD population), well below the BS1 threshold.

Three P/LP variants were unexpectedly observed as **homozygotes** in gnomAD individuals (gnomAD specifically excludes individuals with diagnosed severe pediatric disease, making such observations notable):
- **p.Arg963Ter** (1 homozygote): a well-established European founder variant of moderate frequency.
- **p.Glu1700Gln** (1 homozygote): the East Asian founder variant.
- **p.Arg1172Gln** (1 homozygote): a missense variant warranting functional re-evaluation.

#### Ancestry-specific founder mutations

Stratified analysis across nine gnomAD ancestry groups identified distinct founder mutations in each population, consistent with literature reports:

| Ancestry | Top P/LP variant | Allele frequency |
|----------|-------------------|------------------|
| East Asian | **p.Glu1700Gln** | 0.00339 |
| Admixed American | **p.Gln829Ter** | 0.00057 |
| African/African American | **p.Trp718Ter** | 0.00040 |
| Ashkenazi Jewish | **c.5713-2A>G** (splice acceptor) | 0.00037 |
| Middle Eastern | **p.Glu747Ter** | 0.00033 |
| European (non-Finnish) | **p.Arg963Ter** | 0.00015 |
| European (Finnish) | **p.Arg237Ter** | 0.00009 |
| South Asian | **p.Arg1583Cys** | 0.00004 |

The **p.Glu1700Gln** variant in East Asian populations stands out as the most frequent ancestry-specific P/LP variant, with frequency nearly an order of magnitude higher than equivalent variants in other populations, consistent with a strong founder effect. This variant is reported in literature as one of the most common pathogenic alleles in Chinese and Japanese DFNB9 cohorts (Iwasa et al., *Sci Rep* 2019). The Spanish founder **p.Gln829Ter** appears here as the top variant in the Admixed American population, reflecting the partial Iberian genetic contribution to this gnomAD group; in Spain itself, Q829X accounts for approximately 3% of recessive prelingual deafness cases (Migliosi et al., *J Med Genet* 2002).

### 3. Structural mapping on the AlphaFold-predicted otoferlin structure

Variants with extractable amino acid positions (n = 1,600) were mapped onto the six C2 domains and the transmembrane domain of otoferlin (UniProt Q9HC10 canonical isoform, 1997 aa). Pathogenic variant density (P/LP per residue) was calculated for each domain:

| Domain | Length (aa) | P/LP variants | Density (P/LP per aa) |
|--------|-------------|---------------|------------------------|
| C2B | 121 | 20 | **0.165** (highest) |
| C2D | 115 | 18 | 0.157 |
| C2C | 116 | 18 | 0.155 |
| C2E | 116 | 15 | 0.129 |
| C2F | 117 | 15 | 0.128 |
| Linker (combined) | 1258 | 158 | 0.126 (baseline) |
| Transmembrane | 32 | 3 | 0.094 |
| C2A | 122 | 8 | 0.066 (lowest) |

### 3D visualization of pathogenic variants on otoferlin

The following 360° rotation video shows the AlphaFold-predicted structure of otoferlin (UniProt Q9HC10) with pathogenic and likely pathogenic variants from ClinVar mapped as red spheres. Functional C2 domains are colored individually; linker regions remain in light gray.

https://github.com/matteopira/otof-variants-integrated-analysis/raw/main/results/otof-pymol-mut.mp4

**Domains C2B, C2C, and C2D show approximately 30% higher pathogenic variant density than linker regions**, supporting their critical functional role in calcium-dependent vesicle fusion. Conversely, **C2A shows the lowest pathogenic variant density**, consistent with literature reports describing C2A of otoferlin as **atypical** with respect to canonical C2 domain function. Specifically, its calcium-binding affinity is reportedly weak or absent, and its functional contribution to synaptic transmission may be partially dispensable (Helfmann et al., *J Mol Biol* 2011; Padmanarayana et al., *Biochemistry* 2014).

## Discussion

This integrated analysis confirms several principles of clinical variant interpretation as applied to *OTOF* and refines them with quantitative data:

1. **Loss-of-function variants are reliably classified as pathogenic**, supporting the consistent application of ACMG/AMP PVS1 for *OTOF*. The 92.8% concordance between frameshift consequence and pathogenic classification provides empirical validation of current curation practices.

2. **Population frequency data align with clinical classification**: zero P/LP variants exceed the BS1 threshold in gnomAD, demonstrating that the *OTOF* pathogenic variant pool in ClinVar is consistent with population rarity expectations for an autosomal recessive condition.

3. **Population-specific founder mutations dominate the pathogenic allele pool**: each ancestry group has distinct most-common pathogenic alleles, with strong implications for population-tailored diagnostic panels and gene therapy patient recruitment strategies. The identification of p.Glu1700Gln as a high-frequency East Asian founder is consistent with the geographic origin of the first successful AAV1-hOTOF clinical trials in China (Lv et al., *Lancet* 2024); similarly, p.Gln829Ter has emerged as a target for RNA base-editing therapies aimed at the Spanish DFNB9 population (Xue et al., *Mol Ther* 2023; ClinicalTrials.gov NCT06025032).

4. **C2 functional domains are mutation-intolerant**, with C2B-D showing the highest pathogenic variant density. The atypical behavior of C2A, which shows lower pathogenic density than even the disordered linker regions, likely reflects its non-canonical biology rather than under-curation, and warrants further functional investigation.

### Clinical implications for AAV-mediated gene therapy

The findings have direct implications for ongoing and planned AAV-*OTOF* gene therapy programs. Patients carrying **biallelic frameshift, nonsense, or splice variants** represent the most unambiguously eligible candidates, with essentially universal pathogenic classification and population rarity. Patients carrying **biallelic missense variants**, accounting for a substantial fraction of clinically suspected DFNB9 cases, would benefit most from additional functional characterization (minigene assays, cellular models, humanized mouse studies) before therapeutic eligibility can be confirmed.

**Population-specific screening panels** for newborn hearing screening and pediatric audiology should incorporate ancestry-relevant founder mutations: p.Glu1700Gln for East Asian populations, p.Gln829Ter for Iberian and Latin American populations, p.Arg963Ter for non-Finnish European populations, c.5713-2A>G for Ashkenazi Jewish populations, and the corresponding variants identified for other ancestry groups.

### Limitations

- **ClinVar condition annotation is incomplete**: 1,547 of 2,432 records (63.6%) list the associated condition as "not provided", reflecting submitter-level inconsistencies that limit downstream phenotype-genotype association studies.
- **gnomAD represents adult individuals** and explicitly excludes those with diagnosed severe pediatric disease; the unexpected observation of homozygotes for P/LP variants in three cases warrants functional and clinical follow-up investigation.
- The **AlphaFold-predicted structure** is computational; experimentally determined cryo-EM or X-ray structures of full-length otoferlin are not yet available, although individual C2 domains have been structurally characterized.
- **Short-read sequencing has limitations** for *OTOF* variant detection, including difficulties in phasing variants across the large gene and in detecting deep intronic and complex structural variants. Recent studies advocate the use of long-read sequencing and integrated splice prediction pipelines (Chen et al., *Mol Med* 2025).
- **C2A's apparent tolerance to pathogenic variants** may reflect either authentic biology (atypical C2A function) or under-curation due to lack of functional studies; this distinction cannot be resolved from ClinVar data alone.

## Repository structure

```text
otof-variants-integrated-analysis/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── data/
│   ├── clinvar_result.txt                              # ClinVar export (OTOF variants)
│   ├── gnomad_otof_variants.csv                        # gnomAD v4 export
│   └── otof_alphafold.pdb                              # AlphaFold-predicted structure (Q9HC10)
├── notebooks/
│   ├── 01-exploration.ipynb                            # ClinVar descriptive analysis
│   ├── 02-gnomad-integration.ipynb                     # gnomAD integration & ancestry analysis
│   └── 03-domain-mapping.ipynb                         # Domain-level analysis & lollipop plot
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
    └── otof-pymol-mut.mp4                              # 3D rotation of variants on otoferlin
```

## How to reproduce

1. Clone this repository:
```bash
   git clone https://github.com/matteopira/otof-variants-integrated-analysis.git
   cd otof-variants-integrated-analysis
```

2. Install Python dependencies:
```bash
   pip install -r requirements.txt
```

3. Launch Jupyter and run notebooks in order:
```bash
   jupyter notebook
```
   Then open `notebooks/01-exploration.ipynb`, `02-gnomad-integration.ipynb`, and `03-domain-mapping.ipynb`. For each notebook, use `Kernel → Restart & Run All Cells`.

4. (Optional) For 3D structural visualization, install PyMOL via conda:
```bash
   conda create -n pymol-env -c conda-forge python=3.10 pymol-open-source
   conda activate pymol-env
```

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

## License

This project is released under the **MIT License**. See `LICENSE` for details. Data sources retain their respective licenses (ClinVar: public domain; gnomAD: ODbL; UniProt/AlphaFold: CC-BY 4.0).

## Author

**Matteo Pira**
Medical student, Sapienza University of Rome, Faculty of Pharmacy and Medicine, Master's Degree in Medicine and Surgery, mat. 1881956
GitHub: [@matteopira](https://github.com/matteopira)

## Acknowledgments

This analysis was independently performed as a self-directed computational research project. All data sources are publicly available; analytical workflow and interpretation are original.
