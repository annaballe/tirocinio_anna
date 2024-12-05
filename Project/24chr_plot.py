import matplotlib.pyplot as plt
import pandas as pd
import sys
import gzip

# Get file paths from command line arguments
vcf_file = sys.argv[1]
bed_file = sys.argv[2]

def parse_vcf(file_path):
    """Parse a VCF file and extract chromosome, position, and ADSP fields."""
    data = []
    try:
        with gzip.open(file_path, 'rt') as file:
            for line in file:
                if line.startswith("#"):
                    continue
                columns = line.strip().split("\t")
                if len(columns) <= 9:
                    print(f"Malformed line skipped: {line.strip()}")
                    continue
                chrom, pos, info = columns[0], int(columns[1]), columns[7]
                format_fields = columns[8].split(":")
                sample_fields = columns[9].split(":")
                
                adsp = 0  # Default value if ADSP is not present
                if "ADSP" in format_fields:
                    adsp_index = format_fields.index("ADSP")
                    try:
                        adsp = int(sample_fields[adsp_index].split("/")[0])
                    except (ValueError, IndexError):
                        print(f"Invalid ADSP value at position {pos} in {chrom}")

                data.append((chrom, pos, adsp))
    except Exception as e:
        print(f"Error reading VCF file: {e}")
    return pd.DataFrame(data, columns=["chrom", "pos", "adsp"])

def parse_bed(file_path):
    """Parse a BED file and extract chromosome, start, and end positions."""
    try:
        return pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end"])
    except Exception as e:
        print(f"Error reading BED file: {e}")
        return pd.DataFrame(columns=["chrom", "start", "end"])

# Parse the input files
vcf_data = parse_vcf(vcf_file)
bed_data = parse_bed(bed_file)

# Group data by chromosome for faster access
vcf_grouped = vcf_data.groupby("chrom")
bed_grouped = bed_data.groupby("chrom")

# Prepare the figure
fig, axes = plt.subplots(24, 1, figsize=(10, 40), sharex=True)

for idx, ax in enumerate(axes):
    # Define chromosome label
    if idx < 22:
        label = str(idx + 1)
    elif idx == 22:
        label = "X"
    else:
        label = "Y"

    chrom = f"chr{label}"

    # Get VCF and BED subsets for the current chromosome
    vcf_subset = vcf_grouped.get_group(chrom) if chrom in vcf_grouped.groups else pd.DataFrame(columns=["pos", "adsp"])
    bed_subset = bed_grouped.get_group(chrom) if chrom in bed_grouped.groups else pd.DataFrame(columns=["start", "end"])

    # Plot VCF data (red bars)
    if not vcf_subset.empty:
        ax.bar(vcf_subset["pos"], vcf_subset["adsp"], color="red", label="VCF Source")

    # Plot BED data (blue bars)
    if not bed_subset.empty:
        midpoints = (bed_subset["start"] + bed_subset["end"]) / 2
        ax.bar(midpoints, -1, color="blue", label="BED Source")

    # Add a horizontal reference line at y=0
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)

    # Set subplot title
    ax.set_title(f"Chromosome {label}", fontsize=10)

    # Avoid repeating labels in the legend
    if idx == 0:
        ax.legend()

# Add common labels and adjust layout
plt.xlabel("Position on Chromosome", fontsize=14)
plt.tight_layout()
plt.savefig("chromosome_bar_plot.png", dpi=300)
plt.show()
