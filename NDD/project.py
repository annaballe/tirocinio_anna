import subprocess
import sys
import os

def run_full_project(vcf_folder, affected_txt):
    # Run the full_project.py script with the given arguments
    full_project_command = [
        'python', 'full_project.py', affected_txt, '--folder', vcf_folder
    ]
    
    # Run the command and wait for completion
    subprocess.run(full_project_command, check=True)

def run_genome_plot(cases_file, controls_file, plot_output):
    # Run the genome_plot_v2.py script with the given arguments
    genome_plot_command = [
        'python', 'genome_plot_v2.py', cases_file, plot_output, controls_file
    ]
    
    # Run the command and wait for completion
    subprocess.run(genome_plot_command, check=True)

def main():
    if len(sys.argv) != 4:
        print("Usage: python run_pipeline.py <vcf_folder> <affected_txt> <output_plot.png>")
        sys.exit(1)

    vcf_folder = sys.argv[1]
    affected_txt = sys.argv[2]
    plot_output = sys.argv[3]

    # Check if the affected.txt file exists
    if not os.path.isfile(affected_txt):
        print(f"Error: The file {affected_txt} does not exist.")
        sys.exit(1)

    # Run full_project.py to generate the necessary files
    print("Running full_project.py...")
    run_full_project(vcf_folder, affected_txt)

    # Check if the output files are created
    if not os.path.isfile('cases_1.txt') or not os.path.isfile('controls_0.txt'):
        print("Error: Output files from full_project.py were not generated.")
        sys.exit(1)

    # Run genome_plot_v2.py to generate the plot
    print(f"Generating plot: {plot_output}...")
    run_genome_plot('cases_1.txt', 'controls_0.txt', plot_output)

    print("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
