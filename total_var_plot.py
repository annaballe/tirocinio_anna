import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
import sys

# Function to build the database from CSV files in a directory
def build_database(root_directory):
    """
    Traverse the root directory to find CSV files matching the specified pattern,
    extract relevant metrics, and return a DataFrame.
    """
    file_pattern = re.compile(r'^\d+\.vc_metrics\.csv$')  # Pattern for file matching
    categories_to_extract = ["Total"]  # Only extracting the "Total" category

    data = []  # List to hold rows for the DataFrame

    # Walk through the directory and its subdirectories
    for root, _, files in os.walk(root_directory):
        for file in files:
            if file_pattern.match(file):  # Check if the file matches the pattern
                file_path = os.path.join(root, file)
                sample_match = re.match(r'^(\d+)\.vc_metrics\.csv$', file)
                if sample_match:
                    sample = sample_match.group(1)  # Extract the sample number

                # Read the file and process its contents
                with open(file_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split(",")  # Split the line into columns
                        if len(parts) < 4:  # Skip invalid lines
                            continue
                        if parts[2].strip() in categories_to_extract:  # Filter by category
                            category = parts[2].strip()
                            record = float(parts[3].strip())
                            data.append({"Sample": sample, "Category": category, "Record": record})

    # Convert the collected data to a DataFrame
    df = pd.DataFrame(data, columns=["Sample", "Category", "Record"])
    return df

# Function to create a flower plot for the Total category
def flower_plot(df):
    """
    Create a polar bar plot ("flower plot") for the Total category.
    """
    total_data = df[df['Category'] == 'Total']  # Filter data for the Total category
    if total_data.empty:
        print("No data found for category 'Total'. Exiting.")
        return

    # Normalize heights of bars
    max_value = total_data['Record'].max()
    heights = (total_data['Record'] / max_value) * 70 + 30  # Normalize heights
    width = 2 * np.pi / len(total_data.index)  # Equal width for all bars
    angles = np.arange(len(total_data.index)) * width

    # Generate 200 random colors
    random_colors = [mcolors.to_hex(np.random.rand(3)) for _ in range(len(total_data))]

    # Create the polar plot
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, polar=True)

    bars = ax.bar(
        x=angles,
        height=heights,
        width=width,
        bottom=30,  # Set the radial base
        color=random_colors,
        edgecolor="white",
        linewidth=1.5
    )

    # Add labels to the bars
    label_padding = 5
    for bar, angle, label in zip(bars, angles, total_data["Sample"]):
        rotation = np.rad2deg(angle)
        alignment = "left" if angle < np.pi / 2 or angle >= 3 * np.pi / 2 else "right"
        rotation += 180 if alignment == "right" else 0
        ax.text(
            x=angle,
            y=30 + bar.get_height() + label_padding,
            s=label,
            ha=alignment,
            va='center',
            rotation=rotation,
            rotation_mode="anchor"
        )

    # Set plot title
    ax.set_title("Total filtered variants", size=18, pad=20)
    ax.axis('off')  # Remove grid and axis

    plt.tight_layout()
    output_plot_file = "total_plot.png"
    plt.savefig(output_plot_file)
    #plt.close()
    print(f"Bar plot saved as {output_plot_file}")
    plt.show()

# Main script
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
    
    # Step 1: Build the database
    df = build_database(root_directory)
    if not df.empty:
        # Step 2: Create the bar plot for the extracted data
        flower_plot(df)
