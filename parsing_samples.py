import os
import pandas as pd
import gzip
import sys

def process_vcf_files(file_paths, output_dir='./output'):
    """
    Process multiple VCF .gz files and merge them into a single DataFrame.
    
    Args:
        file_paths: List of paths to VCF .gz files.
        output_dir: Directory where the output file will be saved.
    """
    dataframes = []  # List to hold all filtered DataFrames

    for gz_file in file_paths:
        # Initialize lists to store VCF data
        vcf_data = []
        columns = []

        # Open and read the .gz file line by line
        with gzip.open(gz_file, 'rt') as f:
            for line in f:
                # Extract column names from the line starting with "#"
                if line.startswith("#"):
                    # Only take the header line for columns
                    columns = line.strip().split('\t') if line.startswith("#") and not line.startswith("##") else columns
                else:
                    # Split the data lines by tabs
                    vcf_data.append(line.strip().split('\t'))

        # Convert the VCF data to a pandas DataFrame
        vcf_df = pd.DataFrame(vcf_data, columns=columns)

        # Filtering by dropping specific columns
        filtered_df = vcf_df.drop(columns=['ID', 'QUAL', 'FILTER', 'FORMAT'], errors='ignore')

        # Extract 'REP UNIT' and 'VAR ID' from the 'INFO' column
        filtered_df[['REP UNIT', 'VAR ID']] = filtered_df['INFO'].str.extract(r'RU=(.*?);VARID=(.*);')

        # Drop the original 'INFO' column
        filtered_df = filtered_df.drop(columns=['INFO'], errors='ignore')

        # Extract the last column and process it
        last_column = vcf_df.columns[-1]  # Get the last column in the DataFrame
        filtered_df[last_column] = filtered_df[last_column].str.split(':').str[0]  # Split by colon and take the first part

        # Ensure the last column remains in the last position
        cols = list(filtered_df.columns)
        cols.remove(last_column)
        cols.append(last_column)
        filtered_df = filtered_df[cols]

        # Add the filtered DataFrame to the list
        dataframes.append(filtered_df)

    # Merge all DataFrames using outer join on all common columns
    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = pd.merge(merged_df, df, how='outer', on=list(merged_df.columns.intersection(df.columns)))

    # Replace NaN (missing values) with dots (".")
    merged_df = merged_df.fillna(".")

    # Path to save the merged DataFrame as a .txt file
    output_txt_path = os.path.join(output_dir, 'merged_output.txt')

    # Save the merged DataFrame to a .txt file using tab as a separator
    merged_df.to_csv(output_txt_path, sep='\t', index=False)

    print(f"Merged DataFrame saved to {output_txt_path}")

def main():
    # Check if the script receives the correct number of arguments
    if len(sys.argv) != 4 or sys.argv[2] != '--output_dir':
        print("Usage: python my_script.py <folder_path> --output_dir <output_directory>")
        sys.exit(1)

    # Get folder path and output directory from arguments
    folder_path = sys.argv[1]
    output_dir = sys.argv[3]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get all .gz files from the specified folder
    gz_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.gz')]

    # Ensure there are .gz files to process
    if not gz_files:
        print("No .gz files found in the specified folder.")
        sys.exit(1)

    # Call the processing function
    process_vcf_files(gz_files, output_dir=output_dir)

if __name__ == "__main__":
    main()
