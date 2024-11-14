import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Specify the parent directory containing all folders with CSV files
parent_directory = "/path/to/your/parent_folder"

# Define the column to search and the specific values to look for
target_column = "variant"
target_values = ["total", "another_value"]  # Specify the terms you're interested in

# List to store all gathered data for box plotting
plot_data = []

def process_csv_files(parent_directory, target_column, target_values):
    for root, dirs, files in os.walk(parent_directory):
        # Loop through each file in the current folder
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                try:
                    # Read the CSV file
                    df = pd.read_csv(file_path)
                    
                    # Check if the target column exists in the CSV
                    if target_column not in df.columns:
                        print(f"Warning: Column '{target_column}' not found in file {file_name}")
                        continue  # Skip if target column is missing

                    # Filter rows where the target column matches any of the target values
                    filtered_df = df[df[target_column].isin(target_values)]
                    
                    # If matches found, extract the last column's values
                    if not filtered_df.empty:
                        last_column_name = df.columns[-1]
                        for _, row in filtered_df.iterrows():
                            # Collect information for box plot
                            plot_data.append({
                                "file_name": file_name,
                                "variant": row[target_column],
                                "value": row[last_column_name]
                            })
                
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

def create_box_plot(plot_data, output_plot_path):
    # Convert list of dictionaries to DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # Create a box plot for each variant type
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="variant", y="value", data=plot_df)
    
    # Add title and labels
    plt.title("Box Plot of Values by Variant")
    plt.xlabel("Variant Type")
    plt.ylabel("Value")
    
    # Save or show the plot
    plt.savefig(output_plot_path)
    plt.show()
    print(f"Box plot saved to {output_plot_path}")

# Specify output path for the box plot image
output_plot_path = "/path/to/output_folder/box_plot.png"

# Run the functions
process_csv_files(parent_directory, target_column, target_values)
create_box_plot(plot_data, output_plot_path)
