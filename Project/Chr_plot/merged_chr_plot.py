import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
import re

def natural_sort_chromosomes(chromosome_list):
    def convert(text): 
        return int(text) if text.isdigit() else text.lower()
    def alphanum_key(key): 
        return [convert(c) for c in re.split('([0-9]+)', str(key))]
    return sorted(chromosome_list, key=alphanum_key)

def read_variant_data(file_path):
    try:
        df = pd.read_csv(file_path, sep=r'\s+', engine='python')
        # Print range for each chromosome to verify
        for chrom in natural_sort_chromosomes(df['CHROM'].unique()):
            chrom_data = df[df['CHROM'] == chrom]
            print(f"{chrom}: {chrom_data['POS'].min():,} - {chrom_data['END'].max():,}")
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find the file: {file_path}")

def format_ticks(ax, positions, chrom_data):
    # Get the true min and max positions for this chromosome
    min_pos = positions.min()
    max_pos = positions.max()
    
    # Create fixed number of ticks
    n_ticks = 5
    ticks = np.linspace(min_pos, max_pos, n_ticks)
    
    # Set the ticks and labels
    ax.set_xticks(ticks)
    ax.set_xticklabels(['{:,.1f}M'.format(x/1e6) for x in ticks])
    
    # Set limits to true range
    ax.set_xlim(min_pos, max_pos)

def plot_variant_density_multi_chrom(df):
    chromosomes = natural_sort_chromosomes(df['CHROM'].unique())
    
    n_plots = len(chromosomes)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    fig.suptitle('Variant Density Distribution by Chromosome\nShowing True Genomic Positions', fontsize=16, y=1.02)
    
    axes = axes.flatten() if n_rows > 1 else [axes] if n_cols == 1 else axes
    
    for idx, chrom in enumerate(chromosomes):
        chrom_data = df[df['CHROM'] == chrom]
        
        # Use the actual genomic positions
        positions = pd.concat([
            pd.Series(chrom_data['POS'], name='position'),
            pd.Series(chrom_data['END'], name='position')
        ])
        
        # Create the density plot with actual positions
        sns.kdeplot(
            data=positions,
            fill=True,
            color='#3498db',
            alpha=0.5,
            linewidth=2,
            ax=axes[idx]
        )
        
        # Add rug plot with actual positions
        sns.rugplot(
            data=positions,
            color='#2980b9',
            alpha=0.5,
            height=0.1,
            ax=axes[idx]
        )
        
        # Add position range to title
        start_pos = chrom_data['POS'].min()
        end_pos = chrom_data['END'].max()
        axes[idx].set_title(f'Chromosome {chrom}\n{start_pos:,} - {end_pos:,}')
        axes[idx].set_xlabel('Genomic Position (Mb)')
        axes[idx].set_ylabel('Variant Density')
        
        # Format x-axis with true positions
        format_ticks(axes[idx], positions, chrom_data)
        
        # Add grid
        axes[idx].grid(True, alpha=0.3)
        
        # Add variant count
        variant_count = len(chrom_data)
        axes[idx].text(0.02, 0.98, f'n={variant_count}', 
                      transform=axes[idx].transAxes,
                      verticalalignment='top')
    
    # Remove empty subplots
    for idx in range(len(chromosomes), len(axes)):
        fig.delaxes(axes[idx])
    
    # Add explanation
    fig.text(0.5, -0.02, 
             'Each plot shows the distribution of variants using true genomic coordinates.\n' +
             'X-axis represents the actual genomic position in megabases (Mb), Y-axis shows the density of variants.',
             ha='center', va='center', fontsize=10)
    
    plt.tight_layout()
    return plt

def save_plot(plt, output_path):
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create density plots for variant distribution')
    parser.add_argument('input_file', help='Path to the input file')
    parser.add_argument('-o', '--output', default='variant_density_plot_all_chromosomes.png',
                      help='Output file name (default: variant_density_plot_all_chromosomes.png)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Cannot find {args.input_file}")
        print(f"Current working directory is: {os.getcwd()}")
    else:
        try:
            df = read_variant_data(args.input_file)
            plot = plot_variant_density_multi_chrom(df)
            save_plot(plot, args.output)
            print(f"Plot saved as {args.output}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
