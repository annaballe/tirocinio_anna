import os
import sys
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def find_and_extract_data(root_directory):
    # Define the file name pattern
    file_pattern = re.compile(r'^\d+\.sv_metrics\.csv$')
    # Define the keywords to search for rows with the categories
    keywords = [
        "Number of deletions (PASS)", 
        "Number of insertions (PASS)",
        "Number of duplications (PASS)"
    ]
    
    # Initialize an empty list to hold extracted rows
    collected_data = []
    
    # Traverse the directory and its subdirectories
    for root, dirs, files in os.walk(root_directory):
        for file in files:
            # Check if the file matches the pattern
            if file_pattern.match(file):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                # Open the file and extract rows containing the keywords
                with open(file_path, 'r') as f:
                    for line in f:
                        if any(keyword in line for keyword in keywords):
                            collected_data.append(line.strip().split(','))  # Split by comma
    
    # Convert collected data to a Pandas DataFrame
    if collected_data:
        cd = pd.DataFrame(collected_data)
        return cd
    else:
        print("No matching data found.")
        return pd.DataFrame()

def create_bar_plot(data):
    # Create a DataFrame with necessary columns
    data.columns = ["Type", "ID", "Category", "Value", "Percentage"]
    
    # Filter only the relevant rows for the three categories
    data_filtered = data[data["Category"].isin([
        "Number of deletions (PASS)", 
        "Number of insertions (PASS)", 
        "Number of duplications (PASS)"
    ])]
    
    # Convert 'Value' to numeric (if necessary) and handle missing values
    data_filtered['Value'] = pd.to_numeric(data_filtered['Value'], errors='coerce')
    data_filtered = data_filtered.dropna(subset=['Value'])  # Drop rows with invalid 'Value'

    # Create a bar plot for the three categories
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Category", y="Value", data=data_filtered, palette="Set2", ci=None)
    
    # Custom x-tick labels for the categories
    custom_labels = ["Deletions", "Insertions", "Duplications"]
    plt.xticks(ticks=range(len(custom_labels)), labels=custom_labels, rotation=0)  # Adjust ticks and labels
    
    # Add individual black dots for each data point
    sns.stripplot(x="Category", y="Value", data=data_filtered, color='black', jitter=True, size=6, dodge=True)
    
    # Add titles and labels
    plt.title("Number of Deletions, Insertions, and Duplications")
    plt.ylabel("Count")
    
    # Remove unnecessary xlabel ("Category") because we are using custom x-tick labels
    plt.tight_layout()

    # Save the plot
    output_plot_file = "deletions_insertions_duplications_bar_plot_with_dots.png"
    plt.savefig(output_plot_file)
    plt.close()
    print(f"Bar plot with dots saved as {output_plot_file}")

# Run functions
if __name__ == "__main__":
    # Check if a directory was provided as an argument
    if len(sys.argv) < 2:
        print("Error: No root directory provided.")
        print("Usage: python your_script.py /path/to/your/directory")
        sys.exit(1)
    
    # Get the directory from the command-line argument
    root_directory = sys.argv[1]
    
    # Validate the directory
    if not os.path.isdir(root_directory):
        print(f"Error: {root_directory} is not a valid directory.")
        sys.exit(1)
    
    print(f"Using root directory: {root_directory}")
    
    # Step 1: Find and extract relevant data
    df = find_and_extract_data(root_directory)
    
    if not df.empty:
        # Step 2: Create the bar plot for the extracted data
        create_bar_plot(df)

