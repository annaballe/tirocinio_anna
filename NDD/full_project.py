import os
import sys
import gzip
import shutil
import pandas as pd
import argparse

def read_sample_file(txt_file):
    """ Reads the sample file and creates two groups based on 0 or 1 labels. """
    group_0 = set()  # Control group
    group_1 = set()  # Case group

    with open(txt_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            sample_name, group = parts
            if group == "0":
                group_0.add(sample_name)
            elif group == "1":
                group_1.add(sample_name)

    return group_0, group_1

def move_files(folder, group_0, group_1):
    """ Moves VCF files into two separate folders based on the sample names in affected.txt. """
    # Create directories for controls and cases
    controls_folder = os.path.join(folder, "controls_0")
    cases_folder = os.path.join(folder, "cases_1")
    os.makedirs(controls_folder, exist_ok=True)
    os.makedirs(cases_folder, exist_ok=True)

    # Move the VCF files to the appropriate folder based on sample names
    for filename in os.listdir(folder):
        if filename.endswith(".vcf.gz"):
            sample_name = filename.split('.')[0]
            src_path = os.path.join(folder, filename)

            if sample_name in group_0:
                shutil.move(src_path, os.path.join(controls_folder, filename))
            elif sample_name in group_1:
                shutil.move(src_path, os.path.join(cases_folder, filename))

def process_vcf_files(file_paths, output_file, output_dir):
    """ Processes VCF files and merges them into a single output file. """
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

    merged_df.fillna(".", inplace=True)

    def count_genotypes(row):
        total_count = 0
        for genotype in row[7:]:
            if genotype == "0/1":
                total_count += 1
            elif genotype == "1/1":
                total_count += 2
            elif genotype == "1/2":
                total_count += 1
        return total_count

    merged_df['AC'] = merged_df.apply(count_genotypes, axis=1)

    merged_df['POS'] = merged_df['POS'].astype(int)

    chrom_order = {chrom: i for i, chrom in enumerate(sorted(set(merged_df['CHROM']), key=lambda x: (not x.isdigit(), x)))}
    merged_df['CHROM_ORDER'] = merged_df['CHROM'].map(chrom_order)

    merged_df.sort_values(by=['CHROM_ORDER', 'POS'], ascending=[True, True], inplace=True)
    merged_df.drop(columns=['CHROM_ORDER'], inplace=True)
    merged_df['POS'] = merged_df['POS'].astype(str)

    sample_columns = sorted(list(all_sample_names))
    merged_df = merged_df[['CHROM', 'POS', 'REF', 'ALT', 'END', 'REP_UNIT', 'VAR_ID', 'AC'] + sample_columns]

    output_path = os.path.join(output_dir, output_file)
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"Merged DataFrame saved to {output_path}")

def _process_single_vcf(input_file):
    """ Processes a single VCF file and returns a DataFrame. """
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
            chrom, pos, ref, alt, info = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
            sample_data = fields[9:]

            info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
            repeat_unit = info_dict.get('RU', '')
            ref_repeats = int(info_dict.get('REF', 1))

            alt_alleles = [allele.strip('<>').replace('STR', '') for allele in alt.split(',') if allele.startswith('<STR')]
            allele_repeat_counts = {str(i + 1): int(repeat_count) for i, repeat_count in enumerate(alt_alleles)} if alt_alleles else {}

            if not repeat_unit or not alt_alleles:
                continue

            # Initialize base record
            base_record = {
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'REP_UNIT': repeat_unit,
                'VAR_ID': info_dict.get('VARID', '.')
            }

            # Process each alternative allele
            for alt_idx, alt_repeats in allele_repeat_counts.items():
                alt_record = base_record.copy()
                alt_record['ALT'] = f'STR{alt_repeats}'
                alt_record['END'] = pos + (len(repeat_unit) * alt_repeats)
                
                # Initialize sample genotypes
                sample_variants = {name: '.' for name in sample_names}

                # Process each sample's genotype
                for sample_idx, sample_gt in enumerate(sample_data):
                    genotype = sample_gt.split(':')[0]
                    if genotype in {"./.", "0/0"}:
                        continue

                    alleles = genotype.split('/')
                    
                    # Handle different genotype cases
                    if genotype == f"{alt_idx}/{alt_idx}":
                        sample_variants[sample_names[sample_idx]] = "1/1"
                    elif alt_idx in alleles:
                        other_allele = alleles[1] if alleles[0] == alt_idx else alleles[0]
                        if other_allele == "0":
                            sample_variants[sample_names[sample_idx]] = "0/1"
                        else:
                            sample_variants[sample_names[sample_idx]] = "0/1"

                alt_record.update(sample_variants)
                records.append(alt_record)

    df = pd.DataFrame(records)
    return df, sample_names

def main():
    parser = argparse.ArgumentParser(description="Sort and process VCF files in subfolders.")
    parser.add_argument("txt_file", help="Path to the sample text file")
    parser.add_argument("--folder", required=True, help="Main folder containing subfolders")
    args = parser.parse_args()

    # Step 1: Read the sample file and get group 0 (controls) and group 1 (cases)
    group_0, group_1 = read_sample_file(args.txt_file)

    # Step 2: Move files into control and case folders inside the VCF folder
    move_files(args.folder, group_0, group_1)

    # Step 3: Get the main folder where the script is located
    main_folder = os.path.dirname(os.path.abspath(args.txt_file))

    # Step 4: Process VCF files in the controls_0 and cases_1 directories and save output in the main folder
    process_vcf_files([os.path.join(args.folder, "controls_0")], "controls_0.txt", main_folder)
    process_vcf_files([os.path.join(args.folder, "cases_1")], "cases_1.txt", main_folder)

if __name__ == "__main__":
    main()
