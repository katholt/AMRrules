# Organism-specific interpretation of AMR genotypes
[![DOI](https://zenodo.org/badge/651623725.svg)](https://zenodo.org/badge/latestdoi/651623725)

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of [CLSI](https://clsi.org/) and [EUCAST](https://www.eucast.org/). EUCAST also maintains lists of [expert rules](https://www.eucast.org/expert_rules_and_expected_phenotypes) for some species, including [expected (intrinsic) resistance](https://www.eucast.org/expert_rules_and_expected_phenotypes/expected_phenotypes) and expected susceptibility phenotypes, to guide clinical labs in deciding which drugs to test and whether/how to report them.

We believe there is a similar need for systematic rules for the organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data. 

Current solutions focus on bespoke solutions for specific organisms (e.g. our [Kleborate](https://github.com/klebgenomics/Kleborate) tool for _Klebsiella pneumoniae_; Pathogenwatch [AMR libraries](https://gitlab.com/cgps/pathogenwatch/amr-libraries) for [_Salmonella_ Typhi](https://doi.org/10.1038/s41467-021-23091-2), [_Neisseria gonorrhoeae_](https://doi.org/10.1186/s13073-021-00858-2) and others; [Resfinder 4.0](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) for _E. coli_ and others; [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) and [TBProfiler](https://github.com/jodyphelan/TBProfiler) tools for _Mycobacterium tuberculosis_), but this complicates bioinformatics analyses and promotes fragmentation rather than consolidation of expertise. [AbritAMR](https://github.com/MDU-PHL/abritamr) offers a potential solution for multiple organisms, but organism-specific interpretation rules are hard-coded in Python and separated from supporting evidence, making the logic difficult for others to curate and update.

This repo outlines a proposal for a simple data structure to store organism-specific rules for the interpretation of AMR genotype data, that could be used to enrich the outputs of standard AMR genotyping tools (such as AMRfinderplus and other tools, with or without [hAMRonization](https://github.com/pha4ge/hAMRonization)) and generate informative genome reports that capture expert knowledge about how core genes contribute to antimicrobial susceptibility.


## Data analysis pipeline
![pipeline_image](https://github.com/katholt/orgspecAMR/blob/main/pipeline.png?raw=true)

## Organism-specific rules
![rules_table](https://github.com/katholt/orgspecAMR/blob/main/organism_specific_rules.png?raw=true)

Example file: [organism_specific_rules.txt](https://github.com/katholt/orgspecAMR/blob/main/organism_specific_rules.txt)

## Annotated gene table
![rules_table](https://github.com/katholt/orgspecAMR/blob/main/annotated_gene_report.png?raw=true)

Example file: [annotated_gene_report.txt](https://github.com/katholt/orgspecAMR/blob/main/annotated_gene_report.txt)

## Genome report
![genome_report](https://github.com/katholt/orgspecAMR/blob/main/genome_report.png?raw=true)

Example file (PDF): [genome_report.pdf](https://github.com/katholt/orgspecAMR/blob/main/genome_report.pdf)

Example file (RTF): [genome_report.pdf](https://github.com/katholt/orgspecAMR/blob/main/genome_report.rtf)

## Contributors
This concept was workshopped by members of the [Holt lab](https://holtlab.net) at [London School of Hygiene and Tropical Medicine](https://www.lshtm.ac.uk)
