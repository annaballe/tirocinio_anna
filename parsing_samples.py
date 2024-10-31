import os
import sys
import pandas as pd
import gzip
import subprocess

def process_vcf_files(file_paths, output_file, output_dir='./output'):
    dataframes = []  # List to hold all filtered DataFrames

    for file in file_paths:
        # If a folder is provided, gather all .gz files from it
        if os.path.isdir(file):
            gz_files = [os.path.join(file, f) for f in os.listdir(file) if f.endswith('.gz')]
        else:
            gz_files = [file] if file.endswith('.gz') else []

        for gz_file in gz_files:
            print(f"Processing file: {gz_file}")

            vcf_data = []
            columns = []

            with gzip.open(gz_file, 'rt') as f:
                for line in f:
                    if line.startswith("#"):
                        if not line.startswith("##"):  # Header line
                            columns = line.strip().split('\t')
                            if columns[0].startswith('#'):
                                columns[0] = 'CHROM'
                    else:
                        vcf_data.append(line.strip().split('\t'))

            vcf_df = pd.DataFrame(vcf_data, columns=columns)
            filtered_df = vcf_df.drop(columns=['ID', 'QUAL', 'FILTER', 'FORMAT'], errors='ignore')

            # Extract 'REP UNIT' and 'VAR ID' from the 'INFO' column
            filtered_df[['REP_UNIT', 'VAR_ID']] = filtered_df['INFO'].str.extract(r'RU=(.*?);VARID=(.*);')
            filtered_df = filtered_df.drop(columns=['INFO'], errors='ignore')

            last_column = vcf_df.columns[-1]
            filtered_df[last_column] = filtered_df[last_column].str.split(':').str[0]
            cols = list(filtered_df.columns)
            cols.remove(last_column)
            cols.append(last_column)
            filtered_df = filtered_df[cols]

            dataframes.append(filtered_df)

    if dataframes:
        merged_df = dataframes[0]
        for df in dataframes[1:]:
            merged_df = pd.merge(merged_df, df, how='outer', on=list(merged_df.columns.intersection(df.columns)))

        merged_df = merged_df.fillna(".")

        # Define a function to count occurrences of '1' in the genotypes from the 8th column onward
        def count_ones(row):
            # Count occurrences of '1' in each row starting from the 8th column
            return sum((genotype.count('1'))+(genotype.count('2')) for genotype in row[6:])

        # Apply the function to calculate allele counts and store them in the 'AC' column
        merged_df['AC'] = merged_df.apply(count_ones, axis=1)

        cols = list(merged_df.columns)
        if 'AC' in cols:
            cols.remove('AC')
        cols.insert(6, 'AC')
        merged_df = merged_df[cols]

    output_txt_path = os.path.join(output_dir, output_file)
    merged_df.to_csv(output_txt_path, sep='\t', index=False)
    print(f"Merged DataFrame saved to {output_txt_path}")

    check_tab_separated_columns(output_txt_path)

def check_tab_separated_columns(file_path):
    with open(file_path, 'r') as f:
        header_line = f.readline().strip()

    if '\t' in header_line:
        print("Column names are tab-separated.")
    else:
        print("Warning: Column names are NOT tab-separated.")

def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ['--files', '--folder']:
        print("Usage: python my_script.py --files <file1.gz file2.gz ...> --output <output_file>")
        print("Or: python my_script.py --folder <folder_name> --output <output_file>")
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

    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)

    process_vcf_files(files, output_file, output_dir=output_dir)

if __name__ == "__main__":
    main()
