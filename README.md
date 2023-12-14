# CNAKEPIT

## Definition

> A snake pit is, in a literal sense, a hole filled with snakes. In idiomatic speech, "snake pits" are places of horror, torture and death in European legends and fairy tales. [(from Wikipedia)](https://en.wikipedia.org/wiki/Snake_pit)

In the field of bioinformatics, a cnakepit is now a **C**opy **N**umber **A**lterations detection sna**KE**make **PI**peline for tumor-only **T**argeted sequencing data. In the truest sense of the word, this pit is filled with snakes such as pythons, mambas and anacondas.

## What is it?

This is a pipeline to (attempt to) call CNAs / CNVs by addressing the following challenges:
- requires **no matched** normal sample
- *ideally, requires no reference sample or panel of normals*
- suited for **short reads** *(Illumina)*
- detects (and filters for) **somatic** CNVs
- suitable for **panels** (/targeted sequencing)
- maintained
- publicly used *(>50 citations if not too new)*
- *ideally estimates infer tumor structure and estimate tumor characteristics such as tumor ploidy & purity*

This pipeline was originally forked (and later *unforked*) from [this private repo named FFPE-panel-pipeline](https://github.com/Pregelnuss/FFPE-panel-pipeline) which pre-processes raw reads and later maps them.

## General workflow

- Read preprocessing: QC of the raw data, trimming, QC of the trimmed data
- Mapping and QC
- Variant detection and filtering with [Mutect2](https://www.biorxiv.org/content/10.1101/861054v1.full.pdf) / [GATK](https://gatk.broadinstitute.org/hc/en-us)
- CNV calling by [CNVkit](https://doi.org/10.1371/journal.pcbi.1004873) with different segmentation methods
- Tumor analysis and CNV correction by [PureCN](https://doi.org/10.1186/s13029-016-0060-z):
    - either hierarchical clustering of CNVs or
    - minimal re-segmentation with [PSCBS](https://doi.org/10.1093%2Fbioinformatics%2Fbtr329)

```mermaid
flowchart TD
    A[raw reads] -->|QC & trimming| B[trimmed reads];
    B -->|QC & mapping| C[mapped reads];
    C -->|variant calling| D[variants];
    D -->|variant filtering| E[filtered variants];
    C & E -->|CNV calling by CNVkit| F[provisional CNVs];
    F & E -->|Tumor structure inference and CNV calling by PureCN| G[final CNVs];
```

## Full SnakeMake DAG of pipeline

...in its current configuration (for amplicon seq. data)

![workflow/DAG/DAG.png](https://github.com/pedricolino/map-and-CNV-calling-benchmark/blob/main/workflow/DAG/DAG.png?raw=true)

