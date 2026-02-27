

#!/usr/bin/env python3
"""
Subset a FASTA alignment file based on metabolic potential from annotated CSV.
Creates two output FASTA files:
  1. DNRA sequences (nrfA present, NO denitrification genes)
  2. Denitrifier sequences (nirK OR nosZ present, NO nrfA)
  
Sequences with BOTH pathways are EXCLUDED from both files.
"""

import argparse
import csv
from pathlib import Path

# ================= CONFIGURATION =================
INPUT_CSV = "../data/sca_df_annotated.csv"
INPUT_FASTA = "../data/MSA_TIGR01580_noX_minlength1100_maxlength1350_with_soil_seqs_57.aln-fasta"  
OUTPUT_DNRA_FASTA = "../data/MSA_dnra.aln-fasta"
OUTPUT_DENIT_FASTA = "../data/MSA_denitrifier.aln-fasta"


def parse_fasta_headers(fasta_path):
    """
    Parse FASTA file and extract sequence headers.
    Returns dict: {base_seqid: full_header}
    where base_seqid is the part before the first '|'
    """
    seqid_to_header = {}
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                full_header = line[1:]  # Remove '>'
                # Extract base seqid (before first '|')
                base_seqid = full_header.split('|')[0]
                seqid_to_header[base_seqid] = full_header
    
    print(f"Parsed {len(seqid_to_header)} sequences from FASTA")
    return seqid_to_header

def load_annotated_csv(csv_path):
    """
    Load annotated CSV and categorize sequences by metabolic potential.
    Returns three sets:
      - dnra_only_seqids: seqids with nrfA but NO denitrification genes
      - denit_only_seqids: seqids with nirK/nosZ but NO nrfA
      - both_seqids: seqids with BOTH pathways (excluded from outputs)
    """
    dnra_only_seqids = set()
    denit_only_seqids = set()
    both_seqids = set()
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seqid = row.get('seqid', '')
            if not seqid:
                continue
            
            # Extract base seqid (before first '|')
            base_seqid = seqid.split('|')[0]
            
            # Check DNRA: must have nrfA (strict, not just nirB)
            genes_dnra = row.get('genes_found_dnra', '')
            has_dnra = 'nrfA' in genes_dnra
            
            # Check Denitrification: has nirK OR nosZ
            genes_denit = row.get('genes_found_denitrifier', '')
            has_denit = genes_denit and genes_denit.strip()  # Not empty
            
            # Categorize
            if has_dnra and has_denit:
                both_seqids.add(base_seqid)  # EXCLUDED from both outputs
            elif has_dnra:
                dnra_only_seqids.add(base_seqid)  # Pure DNRA
            elif has_denit:
                denit_only_seqids.add(base_seqid)  # Pure denitrifier
            # else: neither - ignore
    
    print(f"DNRA-only sequences (nrfA, no denit genes): {len(dnra_only_seqids)}")
    print(f"Denitrifier-only sequences (nirK/nosZ, no nrfA): {len(denit_only_seqids)}")
    print(f"Sequences with BOTH pathways (excluded): {len(both_seqids)}")
    
    return dnra_only_seqids, denit_only_seqids, both_seqids

def extract_subset_fasta(input_fasta, output_fasta, target_seqids, seqid_to_header):
    """
    Extract sequences from FASTA file that match target_seqids.
    Writes matching sequences to output_fasta.
    Returns count of sequences written.
    """
    written_count = 0
    current_header = None
    current_sequence = []
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            line = line.rstrip('\n')
            
            if line.startswith('>'):
                # Write previous sequence if it matched
                if current_header and current_sequence:
                    base_seqid = current_header.split('|')[0]
                    if base_seqid in target_seqids:
                        outfile.write(f">{current_header}\n")
                        outfile.write('\n'.join(current_sequence) + '\n')
                        written_count += 1
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header and current_sequence:
            base_seqid = current_header.split('|')[0]
            if base_seqid in target_seqids:
                outfile.write(f">{current_header}\n")
                outfile.write('\n'.join(current_sequence) + '\n')
                written_count += 1
    
    return written_count

def main():
    parser = argparse.ArgumentParser(
        description="Subset FASTA alignment by metabolic potential (DNRA vs Denitrification, exclusive)"
    )
    parser.add_argument("--csv", type=str, default=INPUT_CSV, help="Path to annotated CSV")
    parser.add_argument("--fasta", type=str, default=INPUT_FASTA, help="Path to input FASTA alignment")
    parser.add_argument("--out-dnra", type=str, default=OUTPUT_DNRA_FASTA, help="Output FASTA for DNRA")
    parser.add_argument("--out-denit", type=str, default=OUTPUT_DENIT_FASTA, help="Output FASTA for Denitrifiers")
    args = parser.parse_args()
    
    # Check input files exist
    if not Path(args.csv).exists():
        print(f"Error: CSV file not found: {args.csv}")
        return
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}")
        return
    
    print("=" * 60)
    print("SUBSETTING ALIGNMENT BY METABOLIC POTENTIAL (EXCLUSIVE)")
    print("Sequences with BOTH pathways will be EXCLUDED")
    print("=" * 60)
    
    # Step 1: Load CSV and identify target seqids
    print(f"\n[1/4] Loading annotated CSV: {args.csv}")
    dnra_only_seqids, denit_only_seqids, both_seqids = load_annotated_csv(args.csv)
    
    # Step 2: Parse FASTA headers
    print(f"\n[2/4] Parsing FASTA file: {args.fasta}")
    seqid_to_header = parse_fasta_headers(args.fasta)
    
    # Step 3: Find matches between CSV and FASTA
    dnra_matches = dnra_only_seqids & set(seqid_to_header.keys())
    denit_matches = denit_only_seqids & set(seqid_to_header.keys())
    both_matches = both_seqids & set(seqid_to_header.keys())
    
    print(f"\n[3/4] Matching sequences:")
    print(f"  DNRA-only: {len(dnra_matches)} of {len(dnra_only_seqids)} CSV entries found in FASTA")
    print(f"  Denitrifier-only: {len(denit_matches)} of {len(denit_only_seqids)} CSV entries found in FASTA")
    print(f"  BOTH pathways (excluded): {len(both_matches)} of {len(both_seqids)} found in FASTA")
    
    # Check for missing sequences
    dnra_missing = dnra_only_seqids - set(seqid_to_header.keys())
    denit_missing = denit_only_seqids - set(seqid_to_header.keys())
    both_missing = both_seqids - set(seqid_to_header.keys())
    
    if dnra_missing:
        print(f"\n  ⚠ {len(dnra_missing)} DNRA-only seqids NOT found in FASTA (e.g., {list(dnra_missing)[:3]})")
    if denit_missing:
        print(f"  ⚠ {len(denit_missing)} Denitrifier-only seqids NOT found in FASTA (e.g., {list(denit_missing)[:3]})")
    if both_missing:
        print(f"  ℹ {len(both_missing)} 'both' seqids NOT found in FASTA (already excluded)")
    
    # Report excluded sequences
    if both_matches:
        print(f"\n  🚫 EXCLUDED from both files ({len(both_matches)} sequences with both pathways):")
        print(f"     Examples: {list(both_matches)[:5]}")
    
    # Step 4: Extract subset FASTA files
    print(f"\n[4/4] Writing subset FASTA files...")
    
    dnra_written = extract_subset_fasta(args.fasta, args.out_dnra, dnra_matches, seqid_to_header)
    print(f"  ✓ Wrote {dnra_written} sequences to {args.out_dnra}")
    
    denit_written = extract_subset_fasta(args.fasta, args.out_denit, denit_matches, seqid_to_header)
    print(f"  ✓ Wrote {denit_written} sequences to {args.out_denit}")
    
    # Verify no overlap
    overlap_check = dnra_matches & denit_matches
    if overlap_check:
        print(f"\n  ⚠ WARNING: {len(overlap_check)} sequences appear in BOTH output files!")
        print(f"     This should not happen. Check: {list(overlap_check)[:5]}")
    else:
        print(f"\n  ✓ Verified: No overlap between output files")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Input CSV: {args.csv}")
    print(f"Input FASTA: {args.fasta}")
    print(f"DNRA-only output: {args.out_dnra} ({dnra_written} sequences)")
    print(f"Denitrifier-only output: {args.out_denit} ({denit_written} sequences)")
    print(f"Excluded (both pathways): {len(both_matches)} sequences")
    print(f"Total sequences in FASTA: {len(seqid_to_header)}")
    print(f"Total sequences in outputs: {dnra_written + denit_written}")
    print("\nDone!")

if __name__ == "__main__":
    main()