# AlphaSimR simulation using real maize genomic data

This folder contains scripts for an exploratory recurrent selection simulation using AlphaSimR with empirical maize genotype data.

The simulation uses:

- real genotype matrix from 350 tropical maize inbred lines;
- 13,826 SNP markers;
- a physical map rescaled to approximate genetic positions;
- marker effects previously estimated using BayesB;
- AlphaSimR for crossing, recombination and doubled-haploid generation;
- custom R routines for testcross evaluation, phenotypic error simulation, selection and elite recycling.

## Workflow

### 00_validate_import_AlphaSimR_real_data.R

Validates the import of real inbred genotype data into AlphaSimR using `importInbredGeno()`.

This script checks:

- genotype-map-marker effect alignment;
- import of real inbred genotypes;
- SNP chip registration;
- crossing with `makeCross()`;
- doubled-haploid generation with `makeDH()`;
- SNP genotype extraction with `pullSnpGeno()`;
- genomic value calculation using real BayesB marker effects.

### 01_run_recurrent_simulation_AlphaSimR_real_data.R

Runs the recurrent selection simulation using real genomic data and AlphaSimR.

The simulation includes six scenarios:

- AP;
- AE;
- FF;
- FM;
- MT;
- Culling.

Each scenario is simulated for six cycles. In each cycle, the pipeline includes:

- selection of crosses;
- F1 generation;
- doubled-haploid generation;
- TC1 selection;
- TC2 selection;
- TC3 selection;
- elite selection;
- elite recycling.

### 02_plot_AlphaSimR_results.R

Generates figures and summary tables from the AlphaSimR recurrent simulation outputs.

Main figures include:

- accumulated gain across cycles;
- mean elite performance;
- testcross selection funnel;
- final accumulated gain;
- distribution of selected elites;
- heatmap of elite performance;
- recycled elite dynamics.

## Methodological note

This simulation does not use `runMacs()` or artificial QTLs.

The founder population is empirical and corresponds to the real maize inbred panel. AlphaSimR is used only as the genetic simulation engine for crossing, recombination and doubled-haploid generation.

The genetic values of simulated progenies are calculated using real marker effects previously estimated with BayesB.
