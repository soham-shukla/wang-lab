import pandas as pd
import subprocess
import pysam
import os

# ==========================================
# 1. CONFIGURATION & PATHS
# ==========================================
# Pointing to the Bowtie2 index prefix (No .fa extension)
RRNA_INDEX = "~~/wang lab/mouse_combined_with_rRNA_clean_dedup"

# Your three datasets
TARGET_FILES = [
    "BMDM_longRNA.Top100_rRNA_sequences_RPM.csv",
    "M1_longRNA.Top100_rRNA_sequences_RPM.csv",
    "BMDC_trimmed.Top100_rRNA_sequences_RPM.csv"
]

def run_cmd(cmd):
    """Utility to run shell commands safely via subprocess."""
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

# ==========================================
# 2. STEP 1: PREPARE FASTA FOR ALIGNMENT
# ==========================================
def prepare_queries(csv_path, fasta_out):
    print(f"Loading {csv_path} and generating FASTA queries...")
    df = pd.read_csv(csv_path)
    
    # Dynamically catch column names
    seq_col = 'Sequence' if 'Sequence' in df.columns else 'sequence'
    rank_col = 'Rank' if 'Rank' in df.columns else 'rank' if 'rank' in df.columns else None
    
    if not rank_col:
        df['rank'] = df.index + 1
        rank_col = 'rank'
        
    # Clean sequences: Ensure DNA alphabet for Bowtie2 mapping
    df["sequence_clean"] = df[seq_col].str.upper().str.replace("U", "T", regex=False)
    
    with open(fasta_out, "w") as f:
        for _, row in df.iterrows():
            rank = int(row[rank_col])
            seq = row['sequence_clean']
            f.write(f">{rank}\n{seq}\n")
            
    return df

# ==========================================
# 3. STEP 2: TARGETED BOWTIE2 MAPPING
# ==========================================
def run_targeted_mapping(fasta_in, sam_out):
    print("Running targeted Bowtie2 alignment against parent rRNAs...")
    
    bowtie_cmd = [
        "bowtie2", 
        "-x", RRNA_INDEX, 
        "-f", 
        "-U", fasta_in, 
        "-S", sam_out, 
        "--end-to-end", 
        "--very-sensitive", 
        "-p", "8" 
    ]
    run_cmd(bowtie_cmd)
    print("Bowtie2 mapping complete. Bypassing samtools conversion...")

# ==========================================
# 4. STEP 3: COORDINATE EXTRACTION & FOLDING
# ==========================================
def extract_and_fold(df, sam_in, final_csv_out, folding_temp=37):
    print(f"Extracting mapping coordinates directly from SAM and predicting 2D structures at {folding_temp}°C...")
    
    mapping_data = []
    
    # Parse the raw SAM file directly ("r" mode)
    with pysam.AlignmentFile(sam_in, "r") as sam:
        for read in sam.fetch(until_eof=True):
            rank = int(read.query_name)
            
            if read.is_unmapped:
                mapping_data.append({"rank": rank, "parent_rRNA": "Unmapped", "start": None, "end": None})
                continue
                
            mapping_data.append({
                "rank": rank,
                "parent_rRNA": read.reference_name,
                "start": read.reference_start,
                "end": read.reference_end
            })
            
    # Merge mapping coordinates back into the main dataframe
    map_df = pd.DataFrame(mapping_data).drop_duplicates(subset=['rank'])
    
    rank_col = 'Rank' if 'Rank' in df.columns else 'rank'
    if 'rank' in map_df.columns and rank_col != 'rank':
        map_df = map_df.rename(columns={'rank': rank_col})
        
    merged_df = pd.merge(df, map_df, on=rank_col, how="left")
    
    # Run RNAfold on each sequence
    structures = []
    mfes = []
    seq_col = 'Sequence' if 'Sequence' in df.columns else 'sequence'
    
    for _, row in merged_df.iterrows():
        # Fold the ORIGINAL RNA sequence, not the clean DNA sequence
        seq_rna = row[seq_col] 
        
        process = subprocess.Popen(
            ['RNAfold', '-T', str(folding_temp)], 
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, _ = process.communicate(input=seq_rna)
        lines = stdout.strip().split('\n')
        
        if len(lines) >= 2:
            struct_data = lines[1].rsplit(' ', 1)
            structures.append(struct_data[0])
            mfes.append(float(struct_data[1].strip('() ')))
        else:
            structures.append("Error")
            mfes.append(None)
            
    merged_df["2D_structure"] = structures
    merged_df["minimum_free_energy"] = mfes
    
    merged_df.to_csv(final_csv_out, index=False)
    print(f"Pipeline complete! Final dataset saved to {final_csv_out}")

# ==========================================
# EXECUTION BLOCK
# ==========================================
if __name__ == "__main__":
    for input_csv in TARGET_FILES:
        if not os.path.exists(input_csv):
            print(f"Skipping {input_csv}: File not found.")
            continue
            
        base_name = os.path.basename(input_csv).split('.')[0]
        prefix = f"{base_name}_Structural"
        
        print(f"\n{'='*50}")
        print(f"🚀 Starting pipeline for: {input_csv}")
        print(f"{'='*50}\n")
        
        fasta_file = f"{prefix}_queries.fa"
        sam_file = f"{prefix}.sam"
        final_csv = f"{prefix}_Final_Results.csv"
        
        master_df = prepare_queries(input_csv, fasta_file)
        run_targeted_mapping(fasta_file, sam_file)
        
        # NOTE: Verify the 37 degree parameter with collaborators
        extract_and_fold(master_df, sam_file, final_csv, folding_temp=37) 
        
        # Cleanup temporary text files safely
        if os.path.exists(fasta_file): os.remove(fasta_file)
        if os.path.exists(sam_file): os.remove(sam_file)