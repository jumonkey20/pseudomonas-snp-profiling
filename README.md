# In Silico Profiling of Resistance and Virulence Markers in Pseudomonas aeruginosa

## Overview
This project investigates the genetic divergence of key resistance and outer-membrane porin genes (such as oprD) across different strains of *Pseudomonas aeruginosa*. By leveraging programmatic sequence alignment tools in Python, this workflow automates the detection of Single Nucleotide Polymorphisms (SNPs) and calculates percentage identity matrices between reference and clinical isolates.

## Objectives
* **Automate Data Retrieval:** Curate and manage FASTA-formatted nucleotide sequences from NCBI GenBank.
* **Sequence Alignment:** Implement global pairwise sequence alignments using Biopython.
* **Variant Analysis:** Quantify sequence identity scores and map divergence to understand potential mechanisms of antimicrobial resistance.

## Methodology
1. **Data Collection:** Retrieved full-length coding sequences for target resistance-associated loci from *Pseudomonas aeruginosa* PAO1 (reference) and select variant strains via NCBI.
2. **Computational Pipeline:** Developed a Python-based analysis script utilizing the `Biopython` library (`Bio.Align`).
3. **Scoring Matrix:** Executed global Needleman-Wunsch alignments with optimized gap penalties to evaluate sequence conservation and highlight structural mutation hotspots.

## Tech Stack & Tools
* **Language:** Python 3.x
* **Libraries:** Biopython, Pandas, Matplotlib
* **Databases:** NCBI Nucleotide / GenBank

## Results & Findings
* Successfully established an automated pipeline capable of parsing multi-FASTA files and outputting comparative alignment identity metrics.
* Identified regional sequence variations that correlate with known porin-mediated resistance profiles in clinical variants.

## Future Extensions
* Integrating structural modeling (AlphaFold/PyMOL) to visualize how identified SNP substitutions impact extracellular loops of the target membrane proteins.
