def modify_column(filename):
    # Read the content of the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Modify the content
    for i in range(4, len(lines)):  # Starting from row 5
        if lines[i].startswith("TER"):
            break  # Stop if the line starts with "TER"
        else:
            # Modify column 22 (index 21) to the letter 'A'
            lines[i] = lines[i][:21] + 'A' + lines[i][22:]

    # Write the modified content to a new file
    with open('modified_' + filename, 'w') as file:
        file.writelines(lines)

# Replace 'input.txt' with the name of your input file
modify_column('B_factor.pdb')
