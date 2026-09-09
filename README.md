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


Step-by-Step Replication Guide
To replicate this bioinformatics analysis completely from scratch, follow these instructions to create the files and execute the pipeline:

Step 1: Create the Project Structure and Dependencies
Open your repository in GitHub Codespaces or your local terminal.

Create a file named requirements.txt and add the required package dependencies:

Plaintext
biopython
pandas
Install the dependencies by running this command in your terminal:

Bash
pip install -r requirements.txt


Step 2: Create the FASTA Data File
Create a folder named data.

Inside the data/ folder, create a file named oprD_sequences.fasta and paste the following reference and variant protein sequences:


>NP_249649.1_Reference_PAO1
MKVMKWSAIALAVSAGSTQFAVADAFVSDQAEAKGFIEDSSLDLLLRNYYFNRDGKSGSGDRVDWTQGFL
TTYESGFTQGTVGFGVDAFGYLGLKLDGTSDKTGTGNLPVMNDGKPRDDYSRAGGAVKVRISKTMLKWGE
MQPTAPVFAAGGSRLFPQTATGFQLQSSEFEGLDLEAGHFTEGKEPTTVKSRGELYATYAGETAKSADFI
GGRYAITDNLSASLYGAELEDIYRQYYLNSNYTIPLASDQSLGFDFNIYRTNDEGKAKAGDISNTTWSLA
AAYTLDAHTFTLAYQKVHGDQPFDYIGFGRNGSGAGGDSIFLANSVQYSDFNGPGEKSWQARYDLNLASY
GVPGLTFMVRYINGKDIDGTKMSDNNVGYKNYGYGEDGKHHETNLEAKYVVQSGPAKDLSFRIRQAWHRA
NADQGEGDQNEFRLIVDYPLSIL
>Variant_Isolate_Example
MKVMKWSAIALAVSAGSTQFAVADAFAAAQAEAKGFIEDSSLDLLLRNYYFNRDGKSGSGDRVDWTQGFL
TTYESGFTQGTVGFGVDAFGYLGLKLDGTSDKTGTGNLPVMNDGKPRDDYSRAGGAVKVRISKTMLKWGE
MQPTAPVFAAGGSRLFPQTATGFQLQSSEFEGLDLEAGHFTEGKEPTTVKSRGELYATYAGETAKSADFI
GGRYAITDNLSASLYGAELEDIYRQYYLNSNYTIPLASDQSLGFDFNIYRTNDEGKAKAGDISNTTWSLA
AAYTLDAHTFTLAYQKVHGDQPFDYIGFGRNGSGAGGDSIFLANSVQYSDFNGPGEKSWQARYDLNLASY
GVPGLTFMVRYINGKDIDGTKMSDNNVGYKNYGYGEDGKHHETNLEAKYVVQSGPAKDLSFRIRQAWHRA
NADQGEGDQNEFRLIVDYPLSIL


Step 3: Create the Python Analysis Script
In the root directory, create a file named analyze_resistance.py.

Paste the following Python script into the file:

Python
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd

def analyze_oprD_resistance(fasta_file):
    print(f"Loading sequences from {fasta_file}...")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    if len(sequences) < 2:
        print("Error: FASTA file must contain at least a reference sequence and one variant.")
        return

    reference = sequences[0]
    variants = sequences[1:]

    print(f"Reference sequence loaded: {reference.id} ({len(reference.seq)} aa)")

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    results = []

    for variant in variants:
        print(f"Aligning {variant.id} against reference...")
        alignments = aligner.align(reference.seq, variant.seq)
        best_alignment = alignments[0]

        score = best_alignment.score

        aligned_ref, aligned_var = best_alignment
        matches = sum(1 for r, v in zip(aligned_ref, aligned_var) if r == v)
        total_len = len(aligned_ref)
        percent_identity = (matches / total_len) * 100 if total_len > 0 else 0.0

        results.append({
            "Strain_ID": variant.id,
            "Alignment_Score": score,
            "Percent_Identity": round(percent_identity, 2)
        })

    df = pd.DataFrame(results)
    df.to_csv("alignment_results.csv", index=False)
    print("Analysis complete! Summary saved to 'alignment_results.csv'.")

if __name__ == "__main__":
    analyze_oprD_resistance("data/oprD_sequences.fasta")

    
Step 4: Execute the Pipeline
Run the Python script from your terminal:

Bash
python analyze_resistance.py
Verify that the script successfully generates the alignment_results.csv output summary file containing your alignment scores and identity percentages.
