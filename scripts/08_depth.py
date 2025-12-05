import pysam
import pandas as pd
import glob
import re

def get_rdaf_from_bam(bam_file, ref_fasta, output_csv):
    """Extract read depth allele frequency at all sites"""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        ref = pysam.FastaFile(ref_fasta)
    except Exception as e:
        print(f"Error opening files for {bam_file}: {e}")
        return None
    
    results = []
    
    try:
        for pileupcolumn in bam.pileup():
            pos = pileupcolumn.pos
            ref_name = pileupcolumn.reference_name
            ref_base = ref.fetch(ref_name, pos, pos+1).upper()
            
            # Count bases and total depth (including all mapped reads)
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            total_depth = 0
            for pileupread in pileupcolumn.pileups:
                # Count all mapped reads (skip only soft-clipped or unmapped)
                if not pileupread.is_del and not pileupread.is_refskip:
                    total_depth += 1
                    if pileupread.query_position is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        if base in base_counts:
                            base_counts[base] += 1
            if total_depth == 0:
                continue
            
            # Calculate non-ref RDAF
            non_ref_count = total_depth - base_counts.get(ref_base, 0)
            rdaf = non_ref_count / total_depth
            
            results.append({
                'pos': pos + 1,  # 1-based
                'ref': ref_base,
                'depth': total_depth,
                'rdaf': rdaf,
                'A': base_counts['A'],
                'C': base_counts['C'],
                'G': base_counts['G'],
                'T': base_counts['T']
            })
    finally:
        bam.close()
        ref.close()
    
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Processed {bam_file} -> {output_csv} ({len(df)} sites)")
    return df

def find_reference_for_bam(bam_file):
    """Determine which reference was used for this BAM"""
    # Extract the reference version from BAM filename
    if '68U201_7141_v' in bam_file:
        ref_file = '68U201_7141_v.fa'
    elif '68U201_v' in bam_file:
        ref_file = '68U201_v.fa'
    else:
        raise ValueError(f"Cannot determine reference for {bam_file}")
    
    return ref_file

# Process all BAM files
bam_files = sorted(glob.glob('*.bam'))
print(f"Found {len(bam_files)} BAM files")

for bam in bam_files:
    try:
        ref = find_reference_for_bam(bam)
        output = bam.replace('.bam', '_rdaf.csv')
        get_rdaf_from_bam(bam, ref, output)
    except Exception as e:
        print(f"Error processing {bam}: {e}")

print("Done!")
