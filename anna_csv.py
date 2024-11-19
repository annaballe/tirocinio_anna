import os
import re
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import pylance

def find_and_extract_data(root_directory):
    # Define the file name pattern
    file_pattern = re.compile(r'^\d+\.vc_metrics\.csv$')
    # Define the keywords to search for in rows
    keywords = [",Total,", ",Biallelic,", ",Multiallelic,", ",SNPs,", ",Ti/Tv ratio,"]
    
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
        df = pd.DataFrame(collected_data)
        return df
    else:
        print("No matching data found.")
        print(pd.DataFrame)
        return pd.DataFrame()

def format_millions(x, _):
    """Custom formatter for y-axis to display values in millions."""
    if x >= 1e6:
        return f"{x/1e6:.1f} M"
    return f"{x:.0f}"

def plot_categories_single_image(dataframe):
    # Define the categories and their corresponding colors
    categories = ["Total", "Biallelic", "Multiallelic", "SNPs", "Ti/Tv ratio"]
    colors = ["blue", "green", "orange", "purple", "red"]
    
    # Initialize a list to hold data for plotting
    plot_data = []
    
    # Iterate through each category and prepare the data for plotting
    for category in categories:
        filtered_df = dataframe[dataframe.apply(lambda row: category in row.to_string(), axis=1)]
        
        if not filtered_df.empty:
            numeric_values = pd.to_numeric(filtered_df.iloc[:, 3], errors='coerce').dropna()
            
            for value in numeric_values:
                plot_data.append({'Category': category, 'Value': value})
    
    if not plot_data:
        print("No valid data found for plotting.")
        return
    
    # Create a new DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # Create subplots for each category
    fig, axes = plt.subplots(nrows=1, ncols=len(categories), figsize=(20, 6))
    fig.tight_layout(pad=5.0)

    for i, category in enumerate(categories):
        ax = axes[i]
        cat_data = plot_df[plot_df["Category"] == category]
        
        # Create boxplot for this category
        sns.boxplot(y="Value", data=cat_data, ax=ax, color=colors[i], width=0.5)  # Adjusted box width
        
        # Overlay a strip plot for individual data points
        sns.stripplot(y="Value", data=cat_data, color="black", alpha=0.6, jitter=True, ax=ax, zorder=2)
        
        # Add title and formatting
        ax.set_title(f"{category}", fontsize=14)
        ax.tick_params(axis='y', labelsize=10)
        
        # Customize y-axis to use M notation
        ax.yaxis.set_major_formatter(FuncFormatter(format_millions))
        
        # Add y-axis label to every plot
        ax.set_ylabel("Values", fontsize=12)

    # Overall title
    fig.suptitle("Box Plots by Category", fontsize=16)
    plt.show()


# Main script
if __name__ == "__main__":
    # Replace with the path to your root directory
    root_directory = "/home/anna/Claudia/prova_csv"
    
    # Step 1: Find and extract relevant data
    dataframe = find_and_extract_data(root_directory)
    
    # Step 2: Plot all categories in a single image
    plot_categories_single_image(dataframe)

