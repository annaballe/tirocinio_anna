import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Define paths and keywords
parent_directory = "/home/anna/Claudia/prova_csv"
target_strings = ["Total", "Biallelic", "Multiallelic", "SNPs", "Het/Hom ratio"]

# Dictionary to store extracted data by category
data_by_category = {key: [] for key in target_strings}
data_by_category["Indels"] = []  # Add Indels as a category

def process_csv_files(parent_directory, target_strings):
    for root, dirs, files in os.walk(parent_directory):
        for file_name in files:
            if file_name.endswith(".csv"):
                file_path = os.path.join(root, file_name)
                try:
                    # Read CSV without assuming specific column names
                    with open(file_path, 'r') as file:
                        total_value = 0
                        snps_value = 0
                        found_total = False
                        found_snps = False

                        for line in file:
                            fields = line.strip().split(',')
                            for i, field in enumerate(fields):
                                if field in target_strings:
                                    # Check if there's a field after the target string
                                    if i + 1 < len(fields):
                                        value = float(fields[i + 1])
                                        data_by_category[field].append(value)
                                        
                                        # Track Total and SNPs for calculating Indels
                                        if field == "Total":
                                            total_value = value
                                            found_total = True
                                        elif field == "SNPs":
                                            snps_value = value
                                            found_snps = True
                        
                        # Add Indels to the data_by_category
                        if found_total and found_snps:
                            indels_value = total_value - snps_value
                            data_by_category["Indels"].append(indels_value)
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

def create_bar_plots_with_dots(data_by_category):
    # Convert dictionary to DataFrame for easier plotting
    plot_data = []
    for category, values in data_by_category.items():
        for value in values:
            plot_data.append({"Category": category, "Value": value})
    
    plot_df = pd.DataFrame(plot_data)

    # Exclude the 'Total' category from the DataFrame
    plot_df = plot_df[plot_df["Category"] != "Total"]

    # Rename categories for plotting
    plot_df["Category"] = plot_df["Category"].replace({
        "SNPs": "SNVs",
        "Indels": "Indels"
    })

    # Use a color palette with different shades of blue
    custom_colors = sns.color_palette("Blues", n_colors=5)

    # Create a bar plot with the hue set to 'Category'
    plt.figure(figsize=(12, 8))

    # Create bar plot with custom colors, setting 'hue' to 'Category'
    ax = sns.barplot(x="Category", y="Value", data=plot_df, hue="Category", palette=custom_colors, ci=None, legend=False)

    # Adding dots for each data point (individual values)
    sns.stripplot(x="Category", y="Value", data=plot_df, color='black', jitter=True, dodge=True, ax=ax, alpha=0.6)

    # Adding titles and labels
    plt.title("Bar Plot of Variant Categories with Dots")
    plt.ylabel("Value")
    plt.xlabel("Category")
    plt.xticks(rotation=45)

    # Adding value labels on top of each bar
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.2f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', fontsize=10, color='black', rotation=0,
                    xytext=(0, 9), textcoords='offset points')

    # Save the plot
    output_plot_file = "combined_bar_plots.png"
    plt.tight_layout()
    plt.savefig(output_plot_file)
    plt.close()
    plt.show()
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
    process_csv_files(root_directory, target_strings)
    
    # Step 2: Plot all categories in a single image (bar plots with dots)
    create_bar_plots_with_dots(data_by_category)

