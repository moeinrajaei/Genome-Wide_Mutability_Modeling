# Genome-Wide Mutability Modeling

This repository contains computational workflows for investigating how genomic and sequence features influence mutation susceptibility across the *Caenorhabditis elegans* genome.

The project combines experimentally observed mutations, naturally occurring rare variants, and genomic annotations to identify features associated with mutation probability and to develop statistical models of genome-wide mutability.

## Project Overview

Mutation rates are not uniform across the genome. Local sequence context, nucleotide composition, repetitive sequences, recombination, replication, gene function, and other genomic properties may influence the probability that a genomic position acquires a mutation.

In this project, I developed computational and statistical workflows to quantify these relationships and evaluate how well genomic features predict mutation susceptibility.

The analyses include:

* Processing and quality control of mutation datasets
* Generation of non-mutated genomic sites for comparison
* Annotation of genomic and sequence features
* Sequence-context analysis at multiple motif lengths
* GC-content analysis across multiple genomic scales
* Analysis of mononucleotide repeats and local sequence complexity
* Integration of recombination and DNA replication features
* Functional and gene-level annotation
* Statistical modeling of mutation probability
* Cross-strain model validation
* Model evaluation using ROC curves and AUC
* Visualization and interpretation of genomic predictors of mutability

The primary analyses were performed using **R, Bash, and Linux/HPC computational environments**.

## Biological Questions

This project addresses several related questions:

1. Which genomic and sequence characteristics are associated with increased mutation probability?
2. At what genomic scale is local GC content associated with mutation susceptibility?
3. How strongly does local nucleotide context influence mutability?
4. Are repetitive sequences and local sequence complexity associated with increased mutation probability?
5. Do recombination rate, replication strand, or distance to origins of replication influence mutation susceptibility?
6. Does genomic functional context or germline gene expression contribute to variation in mutability?
7. Are models developed using one *C. elegans* genetic background predictive in another genetic background?
8. Which combination of genomic features provides the strongest model of genome-wide mutation susceptibility?

## Data

The analyses use genomic variation data from *C. elegans*, including:

* Experimentally accumulated mutations
* Naturally occurring rare variants
* Randomly sampled non-mutated genomic positions used as controls
* Reference genome sequence
* Gene and functional annotations
* Recombination and DNA replication annotations
* Germline gene-expression annotations

The mutation datasets include approximately:

* **47,000 experimentally accumulated SNVs/INDELs**
* **1,000,000+ non-mutant sites**

Large raw and intermediate datasets are not necessarily included directly in the repository. Scripts used to generate, annotate, and process these datasets are provided where possible.

## Genomic Features

The modeling framework incorporates features representing local sequence composition, broader genomic context, DNA replication and recombination, functional annotation, and gene expression.

### Mutation and Genomic Coordinates

| Feature    | Description                                                                        |
| ---------- | ---------------------------------------------------------------------------------- |
| `Strain`   | *C. elegans* strain in which the mutation was observed.                            |
| `Mut_type` | Mutation class: single nucleotide variant (`SNV`) or insertion/deletion (`indel`). |
| `CHROM`    | Chromosome containing the mutation.                                                |
| `POS`      | Genomic position of the mutation using a 1-based coordinate system.                |
| `REF`      | Reference allele at the mutation position.                                         |
| `ALT`      | Alternate (mutant) allele observed at the position.                                |

### GC Content

GC content was calculated using five window sizes centered on each genomic position, allowing the relationship between mutation susceptibility and nucleotide composition to be evaluated across local and broader genomic scales.

| Feature       | Description                                                          |
| ------------- | -------------------------------------------------------------------- |
| `GC_11bp`     | GC content within an 11 bp window centered on the mutation site.     |
| `GC_101bp`    | GC content within a 101 bp window centered on the mutation site.     |
| `GC_1001bp`   | GC content within a 1,001 bp window centered on the mutation site.   |
| `GC_10001bp`  | GC content within a 10,001 bp window centered on the mutation site.  |
| `GC_100001bp` | GC content within a 100,001 bp window centered on the mutation site. |

Using multiple window sizes allows the models to distinguish highly localized sequence effects from broader regional genomic effects.

### Local Sequence Context

Multiple representations of the nucleotide sequence surrounding each genomic position were evaluated.

| Feature               | Description                                                                  |
| --------------------- | ---------------------------------------------------------------------------- |
| `motif_3bp_mutfirst`  | 3 bp sequence motif with the mutation site at the first nucleotide position. |
| `motif_3bp_mutlast`   | 3 bp sequence motif with the mutation site at the last nucleotide position.  |
| `motif_3bp_mutcenter` | 3 bp sequence motif centered on the mutation site.                           |
| `motif_5bp_center`    | 5 bp sequence motif centered on the mutation site.                           |
| `motif_7bp_center`    | 7 bp sequence motif centered on the mutation site.                           |
| `motif_9bp_center`    | 9 bp sequence motif centered on the mutation site.                           |

These features allow the models to capture context-dependent differences in mutation probability at multiple sequence scales.

### Repetitive Sequence and Sequence Complexity

| Feature          | Description                                                                                                  |
| ---------------- | ------------------------------------------------------------------------------------------------------------ |
| `seq_complexity` | Sequence composition/complexity metric calculated within a 99 bp window centered on the mutation site.       |
| `mononuc`        | Indicator of whether the genomic position occurs within a mononucleotide repeat tract (`1 = yes`, `0 = no`). |

These features were included to investigate whether repetitive or low-complexity genomic regions have distinct mutation susceptibilities.

### Recombination and DNA Replication

| Feature           | Description                                                                                                                                                               |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `recomb_rate`     | Local recombination rate at the mutation site, measured in centimorgans per megabase (cM/Mb).                                                                             |
| `lead_lag`        | DNA replication strand designation indicating leading- or lagging-strand replication context. Missing values indicate positions for which the annotation was unavailable. |
| `distance_to_ORI` | Distance from the genomic position to the nearest origin of replication (ORI), measured in base pairs.                                                                    |

These variables allow the analysis to test whether genomic processes associated with recombination and DNA replication contribute to mutation-rate variation.

### Functional Genomic Context

| Feature                   | Description                                                                                                                |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `Coding_context`          | Functional genomic region containing the mutation, such as exon, intron, UTR, intergenic, upstream, or downstream regions. |
| `putative_effect`         | Highest predicted functional impact based on SnpEff annotation, including `MODIFIER`, `LOW`, `MODERATE`, and `HIGH`.       |
| `Gene_annotation`         | Gene associated with the genomic position.                                                                                 |
| `gene_expressed_germline` | Indicator of whether the associated gene is expressed in the germline (`1 = yes`).                                         |

These features allow relationships between mutation susceptibility and functional genomic context to be evaluated.

## Analytical Workflow

The overall analytical framework is:

```text
Mutation and rare-variant datasets
                │
                ▼
        Data cleaning and QC
                │
                ▼
Generation of non-mutated control sites
                │
                ▼
       Genomic feature annotation
                │
       ┌────────┼─────────┐
       ▼        ▼         ▼
   Sequence     GC     Functional
    context   content   annotation
       │        │         │
       ├────────┼─────────┤
       │        │         │
       ▼        ▼         ▼
 Recombination  Replication  Gene expression
                │
                ▼
        Feature harmonization
                │
                ▼
        Statistical modeling
                │
        ┌───────┴────────┐
        ▼                ▼
 Logistic regression   Cross-strain
                        validation
        │                │
        └───────┬────────┘
                ▼
          Model evaluation
                │
         ┌──────┼──────┐
         ▼      ▼      ▼
        ROC    AUC   Odds ratios
                │
                ▼
 Biological interpretation
      and visualization
```

## Statistical Modeling

Mutation susceptibility was modeled primarily using **logistic regression**.

Genomic positions were classified as:

```text
1 = mutated site
0 = non-mutated control site
```

Models incorporated combinations of sequence and genomic features to estimate their associations with mutation probability.

Predictor effects were evaluated using regression coefficients and odds ratios, while predictive performance was assessed using receiver operating characteristic (ROC) curves and area under the ROC curve (AUC).

## Cross-Strain Validation

A major objective of the project was to determine whether genomic determinants of mutation susceptibility generalize across genetic backgrounds.

Models were therefore evaluated using reciprocal train/test designs:

```text
Train: PB306
Test:  N2
```

and:

```text
Train: N2
Test:  PB306
```

Models combining information from both genetic backgrounds were also evaluated.

This approach provides an independent assessment of model generalizability rather than relying exclusively on performance within the data used for model development.

## Model Evaluation

Model performance and predictor effects were evaluated using:

* Receiver operating characteristic (ROC) curves
* Area under the ROC curve (AUC)
* Regression coefficients
* Odds ratios
* Cross-strain validation
* Comparison of models incorporating different genomic features

These analyses were used to identify genomic characteristics associated with mutation susceptibility and determine which combinations of features provide the strongest predictive performance.

## Repository Structure

The repository contains separate directories corresponding to major stages and components of the analysis.

```text
Genome-Wide-Mutability-Modeling/
│
├── rawdata/
│   └── Input and intermediate genomic datasets
│
├── data_processing/
│   └── Data cleaning, filtering, and harmonization
│
├── control_sites/
│   └── Generation and processing of non-mutated genomic sites
│
├── feature_annotation/
│   └── Annotation of genomic and sequence features
│
├── sequence_context/
│   └── Sequence motif and local-context analyses
│
├── gc_content/
│   └── GC-content calculations across multiple window sizes
│
├── statistical_models/
│   └── Logistic regression and related statistical analyses
│
├── validation/
│   └── Train/test and cross-strain model validation
│
├── figures/
│   └── Visualization scripts and generated figures
│
└── README.md
```

The exact directory structure may differ from this overview as the repository is being organized to improve clarity and reproducibility.

## Reproducibility

Analyses were primarily performed using:

* **R**
* **Bash**
* **Linux**
* **High-performance computing (HPC)**

R was used for statistical modeling, model evaluation, data visualization, and downstream interpretation.

Command-line workflows were used for genomic data processing, sequence extraction, annotation, and manipulation of large genomic datasets.

## Running the Analysis

The overall analysis can be reproduced approximately in the following order:

```text
1. Prepare and quality-control mutation datasets
2. Generate non-mutated control positions
3. Annotate genomic and sequence features
4. Calculate sequence-context and GC-content features
5. Add recombination, replication, and functional annotations
6. Harmonize mutation and control datasets
7. Fit statistical models
8. Perform cross-strain validation
9. Evaluate predictive performance
10. Generate figures and summary results
```

Individual analysis directories contain scripts corresponding to each stage of the workflow.

Some source and intermediate datasets are too large to distribute directly through GitHub and must be generated locally before running downstream analyses.

## Key Outputs

The workflow generates:

* Harmonized mutation and non-mutated control datasets
* Genomic feature matrices
* Statistical models of mutation susceptibility
* Estimated effects of individual genomic predictors
* Regression coefficients and odds ratios
* ROC curves and AUC estimates
* Cross-strain validation results
* Summary visualizations and publication-quality figures

## Related Publications

This repository supports research investigating mutation rates, mutational spectra, genomic determinants of mutability, and the evolutionary consequences of spontaneous mutations in *Caenorhabditis elegans*. Related publications include:

1. **Rajaei M**, Saxena AS, Johnson LM, Snyder MC, Crombie TA, Tanny RE, Andersen EC, Joyner-Matos J, Baer CF. Mutability of mononucleotide repeats, not oxidative stress, explains the discrepancy between laboratory-accumulated mutations and the natural allele-frequency spectrum in *C. elegans*. **.*Genome Research*.** 2021;31(9):1602–1613.

2. Saber S, Rifat MM, Rahimi F, Snyder M, Singh A, Eickwort B, Newhall Y, **Rajaei M**, Saxena AS, Tanny R, Katju V, Andersen EC, Zhou J, Baer CF. Evolution of the rate, molecular spectrum, and fitness effects of mutation under minimal selection in *Caenorhabditis elegans*. **.*Genetics*.** 2026;233(1).

3. Crombie TA, **Rajaei M**, Saxena AS, Johnson LM, Saber S, Tanny RE, Ponciano JM, Andersen EC, Zhou J, Baer CF. Direct inference of the distribution of fitness effects of spontaneous mutations from recombinant inbred *Caenorhabditis elegans* mutation accumulation lines. **.*Genetics*.** 2024;228(2).

4. Saber S, Snyder M, **Rajaei M**, Baer CF. Mutation, selection, and the prevalence of the *Caenorhabditis elegans* heat-sensitive mortal germline phenotype. **.*G3: Genes|Genomes|Genetics*.** 2022;12(5):jkac063.


## Author

**Moein Rajaei, Ph.D.**
Computational Biologist | Genomics & Bioinformatics Scientist

* GitHub: https://github.com/moeinrajaei
* ORCID: https://orcid.org/my-orcid?orcid=0009-0006-7079-2033
* Google Scholar: https://scholar.google.com/citations?user=LbgtrokAAAAJ&hl=en&oi=ao
* LinkedIn: https://www.linkedin.com/in/moein-rajaei-ph-d-3a5a136b/ 

