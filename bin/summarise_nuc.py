import csv
import argparse
import os
import joblib
import pandas as pd
import xgboost as xgb

def clean_fasta_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    cleaned_content = content.replace('_', '.')
    cleaned_content = remove_gene(cleaned_content, 'Tb927.10.14170')
    cleaned_content = sort_fasta_sequences(cleaned_content)
    with open(file_path, 'w') as file:
        file.write(cleaned_content)

def remove_gene(content, gene_id):
    lines = content.split('\n')
    cleaned_lines = []
    skip = False
    for line in lines:
        if line.startswith('>'):
            if gene_id in line:
                skip = True
            else:
                skip = False
        if not skip:
            cleaned_lines.append(line)
    return '\n'.join(cleaned_lines)

def sort_fasta_sequences(content):
    gene_order = [
        'Tb927.10.13980', 'Tb927.10.15300', 'Tb927.9.10660', 'Tb927.8.3810', 'Tb927.10.4720',
        'Tb927.11.1940', 'Tb927.4.980', 'Tb927.10.1820', 'Tb927.6.3660', 'Tb927.1.3230',
        'Tb927.5.1220', 'Tb927.6.1100', 'Tb927.11.3980', 'Tb927.6.2630', 'Tb927.7.6560',
        'Tb927.10.2810', 'Tb927.3.4350', 'Tb927.6.5000', 'Tb927.8.3480'
    ]
    
    sequences = {}
    current_header = None
    current_sequence = []
    
    for line in content.split('\n'):
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_sequence)
            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)
    
    if current_header:
        sequences[current_header] = ''.join(current_sequence)
    
    sorted_content = []
    for gene_id in gene_order:
        for header, sequence in sequences.items():
            if gene_id in header:
                sorted_content.append(header)
                sorted_content.append(sequence)
    
    return '\n'.join(sorted_content)

def process_dna(header, sequence):
    identifier = header[1:]
    print("\n" + "Analysing " + identifier + "\n")

    sequence = sequence.replace('\n', '').upper()

    if len(sequence) % 3 != 0:
        print("Warning: The length of the DNA sequence is not divisible by 3. Some codons may be incomplete and will be ignored.")

    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    valid_codons = []
    for codon in codons:
        if all(nucleotide in 'ATCG' for nucleotide in codon):
            valid_codons.append(codon)
        else:
            print(f"Warning: Invalid codon '{codon}' found and removed from analysis.")

    print("Codons" + "\n" + str(valid_codons) + "\n")

    counts = {'first': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
              'second': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
              'third': {'A': 0, 'T': 0, 'C': 0, 'G': 0}}

    for codon in valid_codons:
        if len(codon) >= 1 and codon[0] in counts['first']:
            counts['first'][codon[0]] += 1
        if len(codon) >= 2 and codon[1] in counts['second']:
            counts['second'][codon[1]] += 1
        if len(codon) >= 3 and codon[2] in counts['third']:
            counts['third'][codon[2]] += 1

    num_codons = len(valid_codons)
    for position in counts:
        for nucleotide in counts[position]:
            counts[position][nucleotide] /= num_codons

    for position in counts:
        for nucleotide in counts[position]:
            counts[position][nucleotide] = round(counts[position][nucleotide], 5)
    print("Nucleotide counts \n" + str(counts) + "\n")

    average_A = round((counts['first']['A'] + counts['second']['A'] + counts['third']['A']) / 3, 5)
    average_T = round((counts['first']['T'] + counts['second']['T'] + counts['third']['T']) / 3, 5)
    average_C = round((counts['first']['C'] + counts['second']['C'] + counts['third']['C']) / 3, 5)
    average_G = round((counts['first']['G'] + counts['second']['G'] + counts['third']['G']) / 3, 5)

    return {'A': average_A, 'T': average_T, 'C': average_C, 'G': average_G}

def main(dna_file, output_dir):
    clean_fasta_file(dna_file)
    
    with open(dna_file, 'r') as file:
        lines = file.read().strip().split('\n')

    results = {}
    header = None
    sequence = []

    for line in lines:
        if line.startswith('>'):
            if header:
                results[header[1:]] = process_dna(header, ''.join(sequence))
            header = line
            sequence = []
        else:
            sequence.append(line)

    if header:
        results[header[1:]] = process_dna(header, ''.join(sequence))

    rows = ['A', 'T', 'C', 'G']
    output_file = os.path.join(output_dir, f"{os.path.basename(dna_file).split('.')[0]}_counts.csv")

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['isolate'] + [f"{identifier}_{isolate}" for identifier in results for isolate in rows])
        writer.writerow([os.path.basename(dna_file).split('.')[0]] + [round(results[identifier][isolate], 5) for identifier in results for isolate in rows])

    print(f"Results saved to {output_file} \n")

    xgb_mod = joblib.load('ml/model/xgboost_model.pkl')
    new_data = pd.read_csv(output_file)
    new_data_features = new_data.drop(columns=['isolate'])
    new = xgb.DMatrix(new_data_features)
   
    predictions = xgb_mod.predict(new)
    predictions_binary = ['Pleomorphic' if pred > 0.5 else 'Monomorphic' for pred in predictions]
    proba = xgb_mod.predict(new)
    print("This isolate is", predictions_binary[0].upper(), "\nProbability =", proba[0], "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a FASTA file from the MONO-TRAC panel of genes to predict developmental competence.')
    parser.add_argument('dna_file', type=str, help='The path to the DNA sequence file in FASTA format')
    parser.add_argument('output_dir', type=str, help='The directory to save the output CSV file')
    args = parser.parse_args()
    main(args.dna_file, args.output_dir)
