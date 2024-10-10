import os
import sys
import pandas as pd
import gzip

def process_vcf_files(file_paths, output_file, output_dir='./output'):
    """
    Process multiple VCF .gz files and merge them into a single DataFrame.
    
    Args:
        file_paths: List of paths to VCF .gz files or folder.
        output_file: Name of the output file.
        output_dir: Directory where the output file will be saved.
    """
    dataframes = []  # List to hold all filtered DataFrames

    for file in file_paths:
        # If a folder is provided, gather all .gz files from it
        if os.path.isdir(file):
            gz_files = [os.path.join(file, f) for f in os.listdir(file) if f.endswith('.gz')]
        else:
            # If file names are provided, assume they are in the current directory
            gz_files = [file] if file.endswith('.gz') else []

        # Now process each file (whether found in folder or provided directly)
        for gz_file in gz_files:
            print(f"Processing file: {gz_file}")

            vcf_data = []
            columns = []

            # Open and read the .gz file
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

            # Extract the last column and process it
            last_column = vcf_df.columns[-1]
            filtered_df[last_column] = filtered_df[last_column].str.split(':').str[0]
            # Ensure the last column remains in the last position
            cols = list(filtered_df.columns)
            cols.remove(last_column)
            cols.append(last_column)
            filtered_df = filtered_df[cols]

            # Store the filtered DataFrame to list
            dataframes.append(filtered_df)

    # Merge all DataFrames
    if dataframes:
        merged_df = dataframes[0]
        for df in dataframes[1:]:
            merged_df = pd.merge(merged_df, df, how='outer', on=list(merged_df.columns.intersection(df.columns)))

        merged_df = merged_df.fillna(".")

        # Count occurrences of '0/1' in the merged DataFrame
        def count_zeros_ones(row):
            return sum(val == '0/1' for val in row[6:])  # Adjust index if needed

        merged_df['AC'] = merged_df.apply(count_zeros_ones, axis=1)

        # Rearrange the DataFrame to place 'AC' as the 7th column
        cols = list(merged_df.columns)
        if 'AC' in cols:
            cols.remove('AC')
        cols.insert(6, 'AC')  # Insert 'AC' at the 7th position (index 6)
        merged_df = merged_df[cols]

    # Set output file path
    output_txt_path = os.path.join(output_dir, output_file)
    
    # Save merged DataFrame to output file
    merged_df.to_csv(output_txt_path, sep='\t', index=False)
    print(f"Merged DataFrame saved to {output_txt_path}")

    # Check if the merged output file has tab-separated column names
    check_tab_separated_columns(output_txt_path)

def check_tab_separated_columns(file_path):
    """
    Function to check if the first line (header) of the file has tab-separated columns.
    """
    with open(file_path, 'r') as f:
        header_line = f.readline().strip()  # Read the first line (header)

    # Check if the column names are separated by tabs
    if '\t' in header_line:
        print("Column names are tab-separated.")
    else:
        print("Warning: Column names are NOT tab-separated.")

def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ['--files', '--folder']:
        print("Usage: python my_script.py --files <file1.gz file2.gz ...> --output <output_file>")
        print("Or: python my_script.py --folder <folder_name> --output <output_file>")
        sys.exit(1)

    # Get the list of files or folder from command-line arguments
    if sys.argv[1] == '--files':
        files = sys.argv[2:-2]
    elif sys.argv[1] == '--folder':
        folder = sys.argv[2]
        files = [folder]

    # Get output file name from arguments
    if '--output' in sys.argv:
        output_file = sys.argv[sys.argv.index('--output') + 1]
    else:
        output_file = 'merged_output.txt'  # Default name if not provided

    output_dir = './output'  # Default output directory
    os.makedirs(output_dir, exist_ok=True)

    process_vcf_files(files, output_file, output_dir=output_dir)

if __name__ == "__main__":
    main()
