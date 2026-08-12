# Genome-Wide Mutability Modeling

This repository contains computational workflows used to investigate how genomic features influence mutation susceptibility across the genome.

The project evaluates whether local sequence and genomic characteristics can predict the probability that a genomic position accumulates mutations. The analyses combine experimentally observed mutations, naturally occurring rare variants, and genomic annotations to build and validate statistical models of mutability.

## Project Overview

Mutation rates are not uniform across the genome. Local sequence context, nucleotide composition, genomic position, and other features can influence the probability that a site acquires a mutation.

In this project, I developed statistical models to quantify these effects and evaluate how well genomic features predict mutation susceptibility.

The main analyses include:

       - Construction and quality control of mutation datasets
       - Generation of matched non-mutated genomic sites
       - Genomic feature annotation
       - Sequence-context analysis
       - GC-content analysis across multiple genomic window sizes
       - Statistical modeling of mutation probability
       - Independent train/test validation across genetic backgrounds
       - Model evaluation using ROC curves and AUC
       - Visualization and interpretation of genomic predictors of mutability

The primary analyses were performed in R and Bash/Linux using reproducible computational workflows.

## Biological Questions

This project addresses several related questions:

       1- Which genomic characteristics are associated with increased mutation probability?
       2- How strongly does local sequence context affect genome-wide mutability?
       3- Does GC content influence mutation susceptibility, and at what genomic scale?
       4- Are mutability models learned in one genetic background predictive in another?
       5- Can experimentally accumulated mutations and naturally occurring rare variants reveal common genomic determinants of mutation?
## Data

The analyses use genomic variation data from Caenorhabditis elegans, including:

Experimentally accumulated single-nucleotide variants
Rare naturally occurring variants
Randomly sampled non-mutated genomic positions used as controls
Reference genome sequence and genomic annotations

The main mutation datasets include approximately:

7,000 experimentally accumulated SNVs
11,000+ rare variants

Analyses include the N2 and PB306 genetic backgrounds, allowing models trained on one background to be evaluated in an independent genomic background.

Raw or intermediate datasets that are too large for GitHub are not stored directly in the repository. Scripts used to generate or process these files are provided where possible.

Analytical Workflow

The overall workflow is:

Mutation datasets
       │
       ▼
Data cleaning and filtering
       │
       ▼
Generation of matched non-mutated sites
       │
       ▼
Genomic feature annotation
       │
       ├── Sequence context
       ├── GC content
       ├── Strand
       └── Additional genomic features
       │
       ▼
Feature harmonization
       │
       ▼
Statistical modeling
       │
       ├── Logistic regression
       ├── Train/test analyses
       └── Cross-background validation
       │
       ▼
Model evaluation
       │
       ├── ROC curves
       ├── AUC
       └── Odds ratios
       │
       ▼
Biological interpretation and visualization
Genomic Features

Features evaluated in the modeling framework include:

### Sequence Context

Local nucleotide context was incorporated to capture sequence-dependent differences in mutation probability.

Examples include:

Central nucleotide identity
3-bp sequence context
Flanking nucleotide composition

### GC Content

GC content was calculated across multiple genomic scales, including approximately:

10 bp
100 bp
1 kb
10 kb
100 kb

This allows comparison of local versus broader genomic effects on mutation probability.

Additional Features:
Strand information
Mononucleotide sequence features
Genomic position
Reference nucleotide
Genomic annotations
Statistical Modeling

Mutation probability was modeled primarily using logistic regression.

Genomic positions were classified as:

1 = mutated site
0 = non-mutated control site

Models were trained using combinations of genomic features to estimate their association with mutation probability.

Model coefficients were interpreted using odds ratios where appropriate.

Cross-Background Validation

To evaluate whether genomic determinants of mutability generalize across genetic backgrounds, models were trained and tested independently between strains.

Example framework:

Train: PB306
Test: N2

and

Train: N2
Test: PB306

A combined model was also evaluated using data from both genetic backgrounds.

This design provides a stronger test of generalizability than evaluating model performance only on the dataset used for training.

Model Evaluation

Predictive performance was evaluated using:

Receiver operating characteristic curves
Area under the ROC curve (AUC)
Regression coefficients
Odds ratios
Cross-dataset validation

These analyses were used to determine which genomic features contribute most strongly to variation in mutation susceptibility.

Repository Structure

The repository is organized by major stages of the analysis.

Genome-Wide-Mutability-Modeling/
│
├── rawdata/
│   └── Input and intermediate genomic datasets
│
├── data_processing/
│   └── Cleaning, filtering, and harmonization scripts
│
├── control_sites/
│   └── Generation and processing of non-mutated genomic sites
│
├── feature_annotation/
│   └── Annotation of genomic and sequence features
│
├── sequence_context/
│   └── Sequence-context analyses
│
├── gc_content/
│   └── GC-content calculations and analyses
│
├── statistical_models/
│   └── Logistic regression and related modeling workflows
│
├── validation/
│   └── Train/test and cross-background validation
│
├── figures/
│   └── Scripts and outputs used for visualization
│
└── README.md

The exact folder names may differ from this overview as the repository is being reorganized for clarity and reproducibility.

Reproducibility

Analyses were primarily performed using:

R
Bash
Linux
High-performance computing environments

Common R analyses include statistical modeling, model evaluation, data visualization, and downstream interpretation.

Command-line workflows were used for genomic file processing, sequence extraction, annotation, and large-scale data manipulation.

Running the Analysis

The repository is organized so that analyses can be reproduced approximately in the following order:

1. Prepare mutation datasets
2. Generate non-mutated control positions
3. Annotate genomic features
4. Combine mutation and control datasets
5. Fit statistical models
6. Perform cross-background validation
7. Evaluate predictive performance
8. Generate figures and summary results

Individual analysis folders contain scripts for each stage.

Because some source datasets and intermediate files are large, users may need to generate these files locally before running downstream analyses.

Key Outputs

The project generates:

Annotated mutation and control datasets
Genomic feature matrices
Fitted mutability models
Estimated genomic predictor effects
Odds ratios and model coefficients
ROC curves
AUC estimates
Cross-background validation results
Publication-quality figures
Related Publication

This work builds on research into genomic determinants of mutation-rate variation in C. elegans.

Rajaei M., et al.
Mutability of mononucleotide repeats, not oxidative stress, explains the discrepancy between laboratory-accumulated mutations and the natural allele-frequency spectrum in C. elegans.
Genome Research

Author

Moein Rajaei, Ph.D.
Computational Biologist | Cancer Genomics & Bioinformatics

GitHub: [add profile link]
ORCID: [add ORCID link]
Google Scholar: [add Google Scholar link]
LinkedIn: [add LinkedIn link]
License

Add an appropriate open-source license before public release. For scientific code, common choices include MIT, BSD-3-Clause, or GPL-3.0.
