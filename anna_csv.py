import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Specify the parent directory containing all folders with CSV files
parent_directory = "/path/to/your/parent_folder"  # Adjust this path as needed

# Define the column to search and the specific values to look for
target_column = "variant"
target_values = ["total", "another_value", "value1", "value2", "value3", "value4", "value5", "value6", "value7", "value8", "value9", "value10"]  # Modify as needed for your 12 categories

# List to store all gathered data for box plotting
plot_data = []

# Specify output path for the combined box plot image
output_folder = "./output_anna"  # Path to the output directory (same folder where the script is located)

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Full path for the output plot
output_plot_file = os.path.join(output_folder, "combined_box_plots.png")

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

def create_combined_box_plots(plot_data, output_plot_file):
    # Convert list of dictionaries to DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # Create a figure with subplots (3 rows and 4 columns for 12 box plots)
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(16, 12))  # Adjust the figure size as needed
    axes = axes.flatten()  # Flatten the axes array to iterate easily

    # Create a separate box plot for each variant type in target_values
    for idx, variant in enumerate(target_values):
        variant_df = plot_df[plot_df["variant"] == variant]
        
        if variant_df.empty:
            print(f"No data for variant {variant}")
            continue

        # Create box plot for this variant
        sns.boxplot(x="variant", y="value", data=variant_df, ax=axes[idx])
        
        # Add title and labels to each subplot
        axes[idx].set_title(f"Box Plot of {variant}")
        axes[idx].set_xlabel("Variant Type")
        axes[idx].set_ylabel("Value")

    # Adjust layout for better spacing between plots
    plt.tight_layout()

    # Save the combined plot as a single PNG image
    plt.savefig(output_plot_file)
    plt.close()  # Close the plot to free memory
    print(f"Combined box plot saved to {output_plot_file}")

# Run the functions
process_csv_files(parent_directory, target_column, target_values)
create_combined_box_plots(plot_data, output_plot_file)
