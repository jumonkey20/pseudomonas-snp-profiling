from Bio import Align, SeqIO
import pandas as pd

def analyze_sequences(fasta_file):
    """
    Parses a multi-FASTA file, aligns variant sequences against a reference strain,
    calculates identity metrics, and exports a summary CSV table.
    """
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"Error: Could not find the file at {fasta_file}. Check your folder path.")
        return
    
    if len(records) < 2:
        print("Error: FASTA file must contain at least a reference sequence and one variant sequence.")
        return

    ref_record = records[0]  # First sequence is treated as the reference (e.g., PAO1)
    variants = records[1:]   # Remaining sequences are variants/clinical isolates

    print(f"=== Reference Strain Loaded: {ref_record.id} ===")
    print(f"Reference Length: {len(ref_record.seq)} bp\n")

    # Set up global aligner (Needleman-Wunsch algorithm)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    results = []

    # Align each variant against the reference
    for var in variants:
        alignments = aligner.align(ref_record.seq, var.seq)
        best_alignment = alignments[0]
        
        # Calculate percentage identity
        aligned_ref, aligned_var = best_alignment[0], best_alignment[1]
        matches = sum(1 for a, b in zip(aligned_ref, aligned_var) if a == b)
        identity = (matches / len(aligned_ref)) * 100
        
        print(f"Comparing with: {var.id}")
        print(f"  -> Alignment Score: {best_alignment.score}")
        print(f"  -> Sequence Identity: {identity:.2f}%\n")

        results.append({
            "Strain_ID": var.id,
            "Alignment_Score": best_alignment.score,
            "Percent_Identity": round(identity, 2)
        })

    # Save output into a clean CSV table
    df = pd.DataFrame(results)
    df.to_csv("alignment_results.csv", index=False)
    print("Analysis complete! Summary saved to 'alignment_results.csv'.")

if __name__ == "__main__":
    # Point this to your FASTA file location inside the repository
    analyze_sequences("data/oprD_sequences.fasta")
