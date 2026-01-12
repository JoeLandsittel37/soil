import os
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_16S_sequences(annotations_folder, output_fasta):
    """
    Extract ONLY the actual 16S rRNA sequences (not modification enzymes) 
    from annotation files and create a FASTA file.
    
    Args:
        annotations_folder: Path to folder containing annotation .txt files
        output_fasta: Path for output FASTA file
    """
    
    records = []
    
    # Get all text files in the annotations folder
    annotation_files = glob.glob(os.path.join(annotations_folder, "*.txt"))
    
    print(f"Found {len(annotation_files)} annotation files")
    
    IDs = ['ENT01', 'ENT02', 'RLT01', 'PNT01', 'ENT03', 'ENT02', 'ENT01', 'RLT01', 'PDM33', 'PDM17', 'PDM23', 'PDM21', 'PDM22', 'PDM19', 'PDM20', 'ACM05', 'ACM04', 'ACM03', 'ACM02', 'ACM01', 'BOR01', 'BKH01', 'PDM02', 'PDM01', 'PDM03', 'PDM07', 'PDM08', 'PDM09', 'PDM24', 'PDM10', 'PDM25', 'PDM06', 'PDM04', 'PDM05', 'ARZ01', 'DLF01', 'ACV01', 'ACV02', 'VRV01', 'STM01', 'PAR01', 'PAR19367', 'BOR01', 'PDM33', 'PDM17', 'PDM23', 'PDM21', 'PDM19', 'PDM22', 'PDM20', 'PDM18']
    
    for file_path in annotation_files:
        # Get the base filename (e.g., "ACM01" from "ACM01.txt")
        filename = os.path.basename(file_path)
        sample_id = filename.replace('.txt', '')
        
        # if sample_id not in set(IDs):
        #     continue
        
        try:
            with open(file_path, 'r') as f:
                headers = f.readline().strip().split('\t')
                
                # Find the column indices we need
                function_idx = headers.index('function')
                nucleotide_sequence_idx = headers.index('nucleotide_sequence')
                
                found_16s_in_file = False
                
                for line_num, line in enumerate(f, start=2):
                    fields = line.strip().split('\t')
                    
                    if len(fields) <= max(function_idx, nucleotide_sequence_idx): 
                        continue  # Skip incomplete lines
                    
                    function = fields[function_idx]
                    nucleotide_sequence = fields[nucleotide_sequence_idx]
                    
                    # Check if this is specifically the 16S rRNA RNA molecule itself
                    # Looking for: "SSU rRNA ## 16S rRNA, small subunit ribosomal RNA"
                    # or similar patterns that indicate the actual rRNA, not modification enzymes
                    
                    is_actual_16s = (
                        ('SSU rRNA' in function and '16S rRNA' in function) or
                        ('16S ribosomal RNA' in function) or
                        ('small subunit ribosomal RNA' in function and '16S' in function) or
                        (function.startswith('16S rRNA') and 'methyltransferase' not in function and 'processing' not in function and 'nuclease' not in function)
                    )
                    
                    # Also check for negative patterns - exclude modification enzymes
                    is_modification_enzyme = any(
                        term in function.lower() for term in [
                            'methyltransferase', 
                            'processing', 
                            'nuclease',
                            'modification',
                            'modifying',
                            'modifies',
                            'yqgf',
                            'rimm',
                            'transferase',
                            'enzyme'
                        ]
                    )
                    
                    # We want the actual 16S rRNA, not modification enzymes
                    if is_actual_16s and not is_modification_enzyme:
                        
                        # Skip if nucleotide sequence is empty or too short
                        if nucleotide_sequence and len(nucleotide_sequence) > 1300:  # 16S rRNA is typically ~1500 bp
                            # Create sequence record
                            record = SeqRecord(
                                Seq(nucleotide_sequence),
                                id=f"{sample_id}_16S",
                                description=f"16S_rRNA_{sample_id} | {function}"
                            )
                            records.append(record)
                            print(f"Found actual 16S rRNA in {filename}: {function}")
                            found_16s_in_file = True
                            break  # Only take one per file
                        
                        elif nucleotide_sequence and len(nucleotide_sequence) > 100:
                            print(f"Warning: Short 16S sequence in {filename} ({len(nucleotide_sequence)} bp): {function}")
                    
                if not found_16s_in_file:
                    print(f"Note: No actual 16S rRNA found in {filename}")
                            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    # Write all 16S sequences to FASTA file
    if records:
        from Bio import SeqIO
        SeqIO.write(records, output_fasta, "fasta")
        print(f"\nSuccessfully created {output_fasta} with {len(records)} 16S rRNA sequences")
        
        # Print some statistics
        seq_lengths = [len(r.seq) for r in records]
        print(f"Sequence length range: {min(seq_lengths)} - {max(seq_lengths)} bp")
        print(f"Average length: {sum(seq_lengths)/len(seq_lengths):.0f} bp")
        
        # Show some examples
        print("\nSample of extracted sequences:")
        for r in records[:5]:
            print(f"  {r.id}: {len(r.seq)} bp - {r.description}")
    else:
        print("No 16S rRNA sequences found!")
    
    return records


if __name__ == "__main__":
    annotations_folder = "../data/isolate_seqs/annotations"  
    output_fasta = "../data/isolate_seqs/16S_sequences.fasta"
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    
    ssu_sequences = extract_16S_sequences(annotations_folder, output_fasta)