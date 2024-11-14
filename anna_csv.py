import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Specify the parent directory containing all folders with CSV files
parent_directory = "./"  # Change this to the actual path if needed
output_directory = "output_anna"  # Directory where output will be saved

# Define target values to look for in the 3rd column
target_values = ["Biallelic", "Multiallelic", "SNPs", "Ti/Tv ratio"]

# List to store all gathered data for box plotting
plot_data = []

def process_csv_files(parent_directory, target_values):
    for root, dirs, files in os.walk(parent_directory):
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                try:
                    # Read CSV file
                    df = pd.read_csv(file_path, header=None)  # No header row in CSV

                    # Check if the DataFrame has at least 5 columns
                    if df.shape[1] < 5:
                        print(f"File {file_name} has less than 5 columns, skipping.")
                        continue
                    
                    # Filter rows where the 3rd column matches any target value
                    filtered_df = df[df[2].isin(target_values)]  # 3rd column is at index 2
                    
                    # Extract data for the box plot from the 5th column (index 4)
                    for _, row in filtered_df.iterrows():
                        plot_data.append({
                            "file_name": file_name,
                            "variant": row[2],  # 3rd column value
                            "value": row[4]     # 5th column value
                        })
                
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

def create_combined_box_plots(plot_data):
    # Convert list of dictionaries to DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # Create a figure with 5 box plots (5 target values)
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20, 6))
    
    # Loop through each target value to create a box plot
    for idx, variant in enumerate(target_values):
        variant_df = plot_df[plot_df["variant"] == variant]
        
        if variant_df.empty:
            print(f"No data for variant {variant}")
            continue

        # Create box plot for this variant
        sns.boxplot(x="variant", y="value", data=variant_df, ax=axes[idx])
        
        # Set title and labels for each subplot
        axes[idx].set_title(f"Box Plot of {variant}")
        axes[idx].set_xlabel("Variant Type")
        axes[idx].set_ylabel("Value")

    # Adjust layout for better spacing between plots
    plt.tight_layout()

    # Save combined plot to the specified output directory
    output_plot_file = os.path.join(output_directory, "combined_box_plots.png")
    os.makedirs(output_directory, exist_ok=True)  # Ensure output directory exists
    plt.savefig(output_plot_file)
    plt.close()
    print(f"Combined box plot saved to {output_plot_file}")

# Run the functions
process_csv_files(parent_directory, target_values)
create_combined_box_plots(plot_data)
