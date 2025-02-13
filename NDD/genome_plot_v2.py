import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.patches as mpatches

def load_data(filename):
    """Load and validate input data file."""
    try:
        df = pd.read_csv(filename, sep='\t')
        required_columns = ['CHROM', 'POS']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"Missing required columns: {required_columns}")
        return df
    except Exception as e:
        print(f"❌ Error reading file {filename}: {e}")
        sys.exit(1)

def load_cytobands(filename):
    """Load and validate cytoband file."""
    try:
        cytobands = pd.read_csv(filename, sep='\t', 
                               names=['chrom', 'start', 'end', 'band', 'stain'])
        return cytobands
    except Exception as e:
        print(f"❌ Error reading cytoband file {filename}: {e}")
        sys.exit(1)

def get_cytoband_color(stain):
    """Return color for cytoband staining pattern."""
    color_map = {
        'gneg': 'white',
        'gpos100': 'black',
        'acen': 'navy',
        'gvar': 'gray',
    }
    
    if stain in color_map:
        return color_map[stain]
    elif 'gpos' in stain:
        try:
            intensity = int(stain.replace('gpos', '')) / 100
            return (1 - intensity, 1 - intensity, 1 - intensity)
        except:
            return 'gray'
    return 'gray'

def plot_density(ax, positions, scaled_width, current_y, color, label):
    """Plot density with gaps where there's no data"""
    bins = np.linspace(0, scaled_width, 200)
    density, bin_edges = np.histogram(positions, bins=bins)
    
    # Normalize density
    if np.max(density) > 0:
        density = density / np.max(density)
    
    # Find bins with data
    non_zero_bins = density > 0
    
    # Split into segments where there is data
    segments = []
    current_segment = []
    
    for i, (has_data, bin_start) in enumerate(zip(non_zero_bins, bins[:-1])):
        if has_data:
            current_segment.append((bin_start, density[i]))
        elif current_segment:
            segments.append(current_segment)
            current_segment = []
    
    if current_segment:
        segments.append(current_segment)
    
    # Plot each segment separately
    for segment in segments:
        x_values = [x for x, y in segment]
        y_values = [y for x, y in segment]
        ax.plot(x_values, [current_y + 0.5 + y for y in y_values], 
                color=color, linewidth=1.5, label=label)
        label = None  # Only include label for first segment

def create_karyotype_plot(input_file, output_file, control_file=None):
    """Create genome-wide variant distribution plot."""
    # Define cytoband file directly in the script
    cytoband_file = "cytoBands.txt"  # Path to your cytoband file
    
    # Load and validate input files
    print("Loading data files...")
    variants = load_data(input_file)
    cytobands = load_cytobands(cytoband_file)
    controls = load_data(control_file) if control_file else None
    
    # Define chromosome sizes and order
    chromosome_sizes = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, 
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979, 
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, 
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309, 
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, 
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983, 
        "chr22": 50818468, "chrX": 156040895, "chrY": 57227415
    }
    
    # Set up the plot
    plt.rcParams['figure.figsize'] = [20, 30]
    fig, ax = plt.subplots()
    
    # Plot parameters
    autosomes = [f"chr{i}" for i in range(1, 23)]
    sex_chromosomes = ["chrX", "chrY"]
    ordered_chromosomes = autosomes + sex_chromosomes
    max_size = max(chromosome_sizes.values())
    y_spacing = 5
    current_y = (len(ordered_chromosomes) - 1) * y_spacing
    
    print("Creating plot...")
    for chrom in ordered_chromosomes:
        size = chromosome_sizes[chrom]
        relative_width = size / max_size
        scaled_width = max_size * relative_width
        
        # Plot chromosome baseline
        ax.plot([0, scaled_width], [current_y, current_y], color='black', linewidth=1)
        
        # Plot cytobands
        chrom_bands = cytobands[cytobands['chrom'] == chrom]
        for _, band in chrom_bands.iterrows():
            start = band['start'] * relative_width
            width = (band['end'] - band['start']) * relative_width
            color = get_cytoband_color(band['stain'])
            ax.add_patch(plt.Rectangle((start, current_y - 0.4), width, 0.8, 
                                     facecolor=color, edgecolor='black', linewidth=0.5))
        
        # Plot variant density for affected samples
        chrom_variants = variants[variants['CHROM'] == chrom]
        if len(chrom_variants) > 0:
            positions = chrom_variants['POS'].astype(float) * relative_width
            plot_density(ax, positions, scaled_width, current_y, 'dodgerblue', 'Affected')
        
        # Plot variant density for control samples
        if controls is not None:
            chrom_controls = controls[controls['CHROM'] == chrom]
            if len(chrom_controls) > 0:
                positions = chrom_controls['POS'].astype(float) * relative_width
                plot_density(ax, positions, scaled_width, current_y, 'deeppink', 'Control')
        
        # Add chromosome labels
        ax.text(-0.05 * max_size, current_y, chrom, ha='right', va='center')
        current_y -= y_spacing
    
    # Customize plot appearance
    ax.set_xlim(-0.1 * max_size, max_size * 1.1)
    ax.set_ylim(-y_spacing, (len(ordered_chromosomes) + 1) * y_spacing)
    ax.set_xlabel('Genomic Position (bp)')
    ax.set_yticks([])
    
    # Remove unnecessary spines
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    
    # Add title and legends
    plt.title('Genome-wide Variant Distribution', fontsize=50, fontweight='bold', color='navy')
    
    # Create legends
    density_legend = [
        mpatches.Patch(color='dodgerblue', label='Affected Variants'),
        mpatches.Patch(color='deeppink', label='Control Variants')
    ]
    legend1 = plt.legend(handles=density_legend, loc='center right', 
                        title="Variant Density", fontsize=10)
    plt.gca().add_artist(legend1)
    
    cytoband_legend = [
        mpatches.Patch(color=get_cytoband_color(stain), label=stain) 
        for stain in ['gneg', 'gpos100', 'acen', 'gvar', 'gpos75', 'gpos50', 'gpos25']
    ]
    plt.legend(handles=cytoband_legend, loc='center right', title="Cytoband Staining", 
              fontsize=10, bbox_to_anchor=(0.8, 0.5))
    
    # Save plot
    print(f"Saving plot to {output_file}...")
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png', 
                   pad_inches=0.5, facecolor='white')
        print(f"✅ Plot saved successfully")
    except Exception as e:
        print(f"❌ Error saving plot: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) < 3:
        print("Usage: python genome_plot.py <input_file> <output_file> [control_file]")
        sys.exit(1)
    
    # Validate input files
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    control_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    for file in [input_file, control_file]:
        if file and not os.path.exists(file):
            print(f"❌ Input file not found: {file}")
            sys.exit(1)
    
    # Create plot
    create_karyotype_plot(input_file, output_file, control_file)

if __name__ == "__main__":
    main()
