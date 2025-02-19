import os
import sys
from collections import defaultdict

def read_fasta(file_path):
    sequences = defaultdict(str)
    current_isolate = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_isolate = line[1:]
            else:
                sequences[current_isolate] += line
    return sequences

def write_fasta(isolate_sequences, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for isolate, sequences in isolate_sequences.items():
        with open(os.path.join(output_dir, f"{isolate}.fas"), 'w') as file:
            for gene, sequence in sequences.items():
                file.write(f">{gene}\n")
                file.write(f"{sequence}\n")
        print(f"Written file for isolate: {isolate}")

def main(input_dir, output_dir):
    isolate_sequences = defaultdict(lambda: defaultdict(str))
    
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fas'):
            file_path = os.path.join(input_dir, file_name)
            gene_name = file_name.split('_')[0]
            print(f"Processing file: {file_path}")
            sequences = read_fasta(file_path)
            for isolate, sequence in sequences.items():
                isolate_sequences[isolate][gene_name] = sequence
    
    write_fasta(isolate_sequences, output_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python isolate.py <input_directory> <output_directory>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_dir, output_dir)
