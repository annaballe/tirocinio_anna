import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Dictionary to store extracted data by category
target_strings = ["Total", "Biallelic", "Multiallelic", "SNPs", "Ti/Tv ratio", "Het/Hom ratio"]
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

def create_combined_plots(data_by_category):
    # Convert dictionary to DataFrame for easier plotting
    plot_data = []
    for category, values in data_by_category.items():
        for value in values:
            plot_data.append({"Category": category, "Value": value})

    plot_df = pd.DataFrame(plot_data)
    # Exclude the 'Total' and 'Het/Hom ratio' categories from the DataFrame
    plot_df = plot_df[~(plot_df["Category"] == "Total")]


    # Separate 'Het/Hom ratio' from the main categories
    main_categories = plot_df[~(plot_df["Category"] == "Het/Hom ratio")]
    het_hom_ratio = plot_df[plot_df["Category"] == "Het/Hom ratio"]

    # Rename categories for plotting
    main_categories.loc[:,"Category"] = main_categories["Category"].replace({
        "SNPs": "SNVs"
    })

    # Use your vibrant color palette
    custom_colors = ['crimson', "#EF4026", 'orangered', "#FFA500", 'purple']

    # Create subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 16), gridspec_kw={'height_ratios': [3, 1]})

    # Plot 1: Main categories bar plot
    sns.barplot(x="Category", y="Value", data=main_categories, hue="Category", palette=custom_colors, errorbar=None, ax=axes[0])
    sns.stripplot(x="Category", y="Value", data=main_categories, color='black', jitter=True, dodge=True, ax=axes[0], alpha=0.6)
    axes[0].set_title("Variant Categories")
    axes[0].set_ylabel("Value")
    axes[0].set_xlabel("Category")
    axes[0].tick_params(axis='x', rotation=45)

    # Adding value labels on top of each bar
    for p in axes[0].patches:
        axes[0].annotate(f'{p.get_height():.2f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                         ha='center', va='center', fontsize=10, color='black', rotation=0,
                         xytext=(0, 9), textcoords='offset points')

    # Plot 2: Het/Hom ratio
    sns.barplot(x="Category", y="Value", data=het_hom_ratio, hue="Category", palette=["#2a9df4"], errorbar=None, ax=axes[1], width=0.5)
    sns.stripplot(x="Category", y="Value", data=het_hom_ratio, color='black', jitter=True, dodge=True, ax=axes[1], alpha=0.6)
    axes[1].set_title("Het/Hom Ratio")
    axes[1].set_ylabel("Value")
    axes[1].set_xlabel("Category")

    # Save the combined plot
    output_plot_file = "variants_plot.png"
    plt.tight_layout()
    plt.savefig(output_plot_file)
    plt.close()
    print(f"Combined plots saved as {output_plot_file}")

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

    # Step 2: Plot all categories in a single image (bar plots with dots and Het/Hom ratio)
    create_combined_plots(data_by_category)
