import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Define parent directory containing all folders with CSV files
parent_directory = "./"  # Change this path as needed

# Define output directory for the box plots
output_directory = "output_anna"
os.makedirs(output_directory, exist_ok=True)

# Define target values to look for in the 2nd field
target_values = ["Total", "Biallelic", "Multiallelic", "SNPs", "Ti/Tv ratio"]

# Dictionary to store data for each target value
plot_data = {key: [] for key in target_values}

def process_csv_files(parent_directory, target_values):
    # Traverse through each directory and file
    for root, dirs, files in os.walk(parent_directory):
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                
                try:
                    with open(file_path, 'r') as f:
                        for line in f:
                            # Split line by commas
                            fields = line.strip().split(',')

                            # Check if line has enough fields (at least 5)
                            if len(fields) < 5:
                                continue

                            # Check if the 2nd field matches any target value
                            if fields[1] in target_values:
                                # Add the value in the 5th field to the corresponding target's list
                                try:
                                    value = float(fields[4])  # Convert to float for box plotting
                                    plot_data[fields[1]].append(value)
                                except ValueError:
                                    print(f"Non-numeric data in 5th field in {file_name}: {fields[4]}")
                                    continue
                
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

def create_box_plots(plot_data, output_directory):
    # Convert dictionary to DataFrame for plotting
    plot_df = pd.DataFrame(
        [(key, value) for key, values in plot_data.items() for value in values],
        columns=["variant", "value"]
    )

    # Create a figure with subplots (1 row, 5 columns for 5 box plots)
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20, 6))
    axes = axes.flatten()

    # Plot each variant type in target_values
    for idx, (variant, values) in enumerate(plot_data.items()):
        if values:
            sns.boxplot(y=plot_df[plot_df["variant"] == variant]["value"], ax=axes[idx])
            axes[idx].set_title(f"{variant} Box Plot")
            axes[idx].set_ylabel("Value")
        else:
            print(f"No data for variant {variant}")

    # Save combined plot as a PNG image
    combined_output_path = os.path.join(output_directory, "combined_box_plots.png")
    plt.tight_layout()
    plt.savefig(combined_output_path)
    plt.close()  # Close the plot to free memory
    print(f"Combined box plot saved to {combined_output_path}")

# Run the functions
process_csv_files(parent_directory, target_values)
create_box_plots(plot_data, output_directory)

