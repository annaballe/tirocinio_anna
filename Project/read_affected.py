import os
import sys

def outcome_patient(affected_txt, output_txt):
    """
    Processes each line from the text file and determines if the patient is affected or not.
    Outputs results to a specified file.
    """
    with open(affected_txt, 'r') as file, open(output_txt, 'w') as out_file:
        for line in file:
            # Strip any extra whitespace
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Split each line by tab to separate the patient ID and the status
            patient_info = line.split('\t')
            
            # Check if the line has exactly 2 fields (ID and status)
            if len(patient_info) == 2:
                patient_id, status = patient_info[0].strip(), patient_info[1].strip()  # Correctly split and assign fields
                
                # Check the status and write the appropriate message to the output file
                if status == '0':
                    out_file.write(f'The patient: {patient_id} is not affected\n')
                elif status == '1':
                    out_file.write(f'The patient: {patient_id} is affected\n')
                else:
                    out_file.write(f'Unknown status for patient: {patient_id}\n')
            else:
                out_file.write(f'Invalid line format: {line}\n')  # Inform of incorrect format

if __name__ == "__main__":
    # Ensure that the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python outcome_patient.py <input_file_name> <output_file_name>")
        sys.exit(1)

    # Extracting input and output filenames from command-line arguments
    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]

    # Get the current working directory
    current_directory = os.getcwd()

    # Construct full paths for the input and output files
    affected_txt = os.path.join(current_directory, input_file_name)
    output_txt = os.path.join(current_directory, output_file_name)

    # Call the function to process the file
    outcome_patient(affected_txt, output_txt)
