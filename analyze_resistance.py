from Bio import Align, SeqIO

# 1. Load reference and variant sequences
records = list(SeqIO.parse("data/oprD_sequences.fasta", "fasta"))
ref_record = records[0]  # Reference strain (e.g., PAO1)
variants = records[1:]   # Clinical/variant strains

print(f"Reference Strain: {ref_record.id}")

# 2. Set up global aligner
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

# 3. Align each variant against the reference
for var in variants:
    alignments = aligner.align(ref_record.seq, var.seq)
    best_alignment = alignments[0]
    
    # Calculate identity percentage
    matches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b)
    identity = (matches / len(best_alignment[0])) * 100
    
    print(f"\n--- Alignment: {ref_record.id} vs {var.id} ---")
    print(f"Sequence Identity: {identity:.2f}%")
    print(f"Alignment Score: {best_alignment.score}")
