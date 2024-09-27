def process_files(atom_file, replacement_file, output_file):
    with open(atom_file, 'r') as af, open(replacement_file, 'r') as rf, open(output_file, 'w') as of:
        # Read all lines from both files
        atom_lines = af.readlines()
        replacement_lines = rf.readlines()

        # Create a dictionary to store replacement data with key as (column 5, column 2)
        replacement_dict = {}
        for line in replacement_lines:
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"Invalid format in replacement file: {line}")
            
            # Construct key from replacement file (column 1, column 2)
            key = (parts[0].strip(), parts[1].strip())
            replacement_dict[key] = parts[2].strip()

        # Print the replacement dictionary for debugging
        print("Replacement Dictionary:")
        for k, v in replacement_dict.items():
            print(f"{k}: {v}")

        # Process each line in the atom file
        for line in atom_lines:
            # Skip lines that do not start with 'ATOM'
            if not line.startswith("ATOM"):
                of.write(line)
                continue

            # Extract the relevant columns by their fixed positions
            atom_serial = line[6:11].strip()      # Column 2: atom serial number
            chain_name = line[21:22].strip()      # Column 5: chain name

            key = (chain_name, atom_serial)

            # Print the key being checked for debugging
            print(f"Checking key: {key}")

            if key in replacement_dict:
                # Replace the B-factor value in the correct column
                new_b_factor = replacement_dict[key].rjust(6)
                updated_line = line[:60] + new_b_factor + line[66:]
                of.write(updated_line)
            else:
                # Log a warning for no match and write the original line
                print(f"Warning: No replacement found for atom file line: {line.strip()}")
                of.write(line)

if __name__ == "__main__":
    # Replace 'atom_file.txt', 'replacement_file.txt', and 'output_file.txt' with your actual file paths.
    process_files('modified_B_factor_full2.pdb', 'input_insertion_py_full.txt', 'B_factor_averaged_full.pdb')
