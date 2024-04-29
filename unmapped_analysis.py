# Ver de donde se pierden las variantes:
def count_first_characters(file_path):
    # Initialize a dictionary to store the counts of each first character
    char_counts = {}

    # Open the file and iterate through its lines
    with open(file_path, 'r') as file:
        for line in file:
            # Extract the first character of the line
            first_char = line.strip()[0]

            # Update the count for this first character in the dictionary
            char_counts[first_char] = char_counts.get(first_char, 0) + 1

    return char_counts

# Example usage:
file_path = 'C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic_map/GSE33528_2.map'  # Replace 'your_file.txt' with the path to your file
char_counts_og = count_first_characters(file_path)

# Print the counts for each first character
for char, count in char_counts_og.items():
    print(f"Character '{char}': {count} occurrences")