import os
import sys
import gzip
import pandas as pd

def process_vcf_files(file_paths, output_file, output_dir='./output'):
    """
    Process VCF files, split multiallelic records, and create a consolidated DataFrame
    """
    dataframes = []
    all_sample_names = set()

    for file in file_paths:
        if os.path.isdir(file):
            gz_files = [os.path.join(file, f) for f in os.listdir(file) if f.endswith('.gz')]
        else:
            gz_files = [file] if file.endswith('.gz') else []

        for gz_file in gz_files:
            print(f"Processing file: {gz_file}")
            df, sample_names = _process_single_vcf(gz_file)
            dataframes.append(df)
            all_sample_names.update(sample_names)

    if not dataframes:
        print("No valid VCF files processed.")
        return

    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = pd.merge(merged_df, df, how='outer', 
                              on=['CHROM', 'POS', 'REF', 'ALT', 'END', 'REP_UNIT', 'VAR_ID'])

    merged_df = merged_df.fillna(".")

    # Count alternate alleles
    def count_ones(row):
        return sum((genotype.count('1')) + (genotype.count('2')) for genotype in row[7:])

    merged_df['AC'] = merged_df.apply(count_ones, axis=1)

    # Sort and prepare output
    #merged_df['POS'] = merged_df['POS'].astype(int)
    #merged_df.sort_values(by=['POS'], ascending=[True], inplace=True)
    #merged_df['POS'] = merged_df['POS'].astype(str)
    ######
    # Convert POS to integer
    merged_df['POS'] = merged_df['POS'].astype(int)

    # Ensure CHROM is sorted properly (handling 'chr' and numerical values correctly)
    chrom_order = {chrom: i for i, chrom in enumerate(sorted(set(merged_df['CHROM']), key=lambda x: (not x.isdigit(), x)))}
    merged_df['CHROM_ORDER'] = merged_df['CHROM'].map(chrom_order)

    # Sort by CHROM_ORDER and POS
    merged_df.sort_values(by=['CHROM_ORDER', 'POS'], ascending=[True, True], inplace=True)

    # Drop temporary CHROM_ORDER column
    merged_df.drop(columns=['CHROM_ORDER'], inplace=True)

    # Convert POS back to string (if needed)
    merged_df['POS'] = merged_df['POS'].astype(str)
    ######
    # Ensure sample columns are in correct order
    sample_columns = sorted(list(all_sample_names))
    merged_df = merged_df[['CHROM', 'POS', 'REF', 'ALT', 'END', 'REP_UNIT', 'VAR_ID', 'AC'] + sample_columns]

    # Create output directory if not exists
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"Merged DataFrame saved to {output_path}")
    _check_tab_separated_columns(output_path)

def _process_single_vcf(input_file):
    """
    Process a single VCF file and return a DataFrame with split records
    """
    open_func = gzip.open if input_file.endswith(".gz") else open

    records = []
    columns = None
    sample_names = []

    with open_func(input_file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                if not line.startswith("##") and columns is None:
                    columns = line.strip().split('\t')
                    if columns[0].startswith('#'):
                        columns[0] = 'CHROM'
                    sample_names = columns[9:]
                continue

            fields = line.strip().split('\t')
            
            # Parse VCF record
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            sample_data = fields[9:]

            # Extract key information from INFO field
            info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)

            repeat_unit = info_dict.get('RU', '')
            ref_repeats = int(info_dict.get('REF', 1))

            # Extract and parse ALT field 
            alt_alleles = [allele.strip('<>').replace('STR', '') for allele in alt.split(',') if allele.startswith('<STR')]
            
            # Map allele indexes to repeat counts
            allele_repeat_counts = {str(i + 1): int(repeat_count) for i, repeat_count in enumerate(alt_alleles)} if alt_alleles else {}

            # Process samples
            sample_variants = {name: '.' for name in sample_names}

            for sample_idx, sample_gt in enumerate(sample_data):
                genotype = sample_gt.split(':')[0]

                # Skip records with no variation or homozygous reference
                if genotype in {"./.", "0/0"}:
                    continue

                if not repeat_unit:
                    continue

                repeat_unit_length = len(repeat_unit)

                # Process each allele in the genotype
                alleles = genotype.split('/')
                unique_alleles = set(alleles) - {"0"}  # Exclude reference allele

                for allele in unique_alleles:
                    if allele.isdigit() and (not allele_repeat_counts or allele in allele_repeat_counts):
                        alt_repeats = allele_repeat_counts.get(allele, ref_repeats)
                        total_repeat_length = repeat_unit_length * alt_repeats
                        end_pos = pos + total_repeat_length

                        record = {
                            'CHROM': chrom,
                            'POS': pos,
                            'REF': ref,
                            'ALT': f'STR{alt_repeats}',
                            'END': end_pos,
                            'REP_UNIT': repeat_unit,
                            'VAR_ID': info_dict.get('VARID', '.')
                        }
                        # Update sample variant column
                        sample_name = sample_names[sample_idx]
                        sample_variants[sample_name] = '1'
                        
                        # Add sample variants to the record
                        record.update(sample_variants)
                        
                        records.append(record)

    df = pd.DataFrame(records)
    return df, sample_names

def _check_tab_separated_columns(file_path):
    """
    Check if column names are tab-separated
    """
    with open(file_path, 'r') as f:
        header_line = f.readline().strip()

    if '\t' in header_line:
        print("Column names are tab-separated.")
    else:
        print("Warning: Column names are NOT tab-separated.")

def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ['--files', '--folder']:
        print("Usage: python script.py --files <file1.gz file2.gz ...> --output <output_file>")
        print("Or: python script.py --folder <folder_name> --output <output_file>")
        sys.exit(1)

    if sys.argv[1] == '--files':
        files = sys.argv[2:-2]
    elif sys.argv[1] == '--folder':
        folder = sys.argv[2]
        files = [folder]

    if '--output' in sys.argv:
        output_file = sys.argv[sys.argv.index('--output') + 1]
    else:
        output_file = 'merged_output.txt'

    process_vcf_files(files, output_file)

if __name__ == "__main__":
    main()
