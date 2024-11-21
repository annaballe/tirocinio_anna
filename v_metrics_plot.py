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
                    with open(file_path, 'r') as file:
                        total_value = None
                        snps_value = None

                        for line in file:
                            fields = line.strip().split(',')
                            for i, field in enumerate(fields):
                                if field in target_strings:
                                    # Check if there's a field after the target string
                                    if i + 1 < len(fields):
                                        try:
                                            value = float(fields[i + 1])
                                            data_by_category[field].append(value)

                                            # Track Total and SNPs for calculating Indels
                                            if field == "Total":
                                                total_value = value
                                            elif field == "SNPs":
                                                snps_value = value
                                        except ValueError:
                                            print(f"Invalid value for {field} in file {file_name}. Skipping.")
                                            
                        # Add Indels to the data_by_category if both Total and SNPs are valid
                        if total_value is not None and snps_value is not None:
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

    # Exclude 'Total' from the main plot
    main_categories = plot_df[~plot_df["Category"].isin(["Total", "Het/Hom ratio", "Ti/Tv ratio"])]
    het_hom_ratio = plot_df[plot_df["Category"] == "Het/Hom ratio"].copy()
    ti_tv_ratio = plot_df[plot_df["Category"] == "Ti/Tv ratio"].copy()

    # Rename categories for plotting
    main_categories = main_categories.copy()
    main_categories["Category"] = main_categories["Category"].replace({
        "SNPs": "SNVs"
    })

    # Define custom color palette
    custom_colors = ['crimson', "#EF4026", 'orangered', "#FFA500", 'purple']

    # Create subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 16), gridspec_kw={'height_ratios': [3, 1]})

    # Plot 1: Main categories bar plot
    sns.barplot(
        x="Category", y="Value", data=main_categories, palette=custom_colors,
        errorbar=None, ax=axes[0]
    )
    sns.stripplot(
        x="Category", y="Value", data=main_categories, color='black',
        jitter=True, dodge=True, ax=axes[0], alpha=0.6
    )
    axes[0].set_title("Variant Categories")
    axes[0].set_ylabel("Value")
    axes[0].set_xlabel("Category")
    axes[0].tick_params(axis='x', rotation=45)

    # Plot 2: Het/Hom and Ti/Tv ratios
    ratios_data = pd.concat([het_hom_ratio, ti_tv_ratio])
    sns.barplot(
        x="Category", y="Value", data=ratios_data,
        palette=["#2a9df4", "#FF6347"], errorbar=None, ax=axes[1]
    )
    sns.stripplot(
        x="Category", y="Value", data=ratios_data, color='black',
        jitter=True, dodge=True, ax=axes[1], alpha=0.6
    )
    axes[1].set_title("Het/Hom and Ti/Tv Ratios")
    axes[1].set_ylabel("Value")
    axes[1].set_xlabel("Category")
    axes[1].tick_params(axis='x', rotation=45)

    # Save the combined plot
    output_plot_file = "variants_plot.png"
    plt.tight_layout()
    plt.savefig(output_plot_file)
    plt.close()
    print(f"Combined plots saved as {output_plot_file}")

# Run functions
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: No root directory provided.")
        print("Usage: python your_script.py /path/to/your/directory")
        sys.exit(1)

    root_directory = sys.argv[1]

    if not os.path.isdir(root_directory):
        print(f"Error: {root_directory} is not a valid directory.")
        sys.exit(1)

    print(f"Using root directory: {root_directory}")

    process_csv_files(root_directory, target_strings)
    create_combined_plots(data_by_category)
