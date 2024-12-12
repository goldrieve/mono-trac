import csv
import argparse

def process_dna(header, sequence):
    identifier = header[1:]
    print("\n" + "Analysing " + identifier + "\n")

    sequence = sequence.replace('\n', '').upper()

    if len(sequence) % 3 != 0:
        raise ValueError("The length of the DNA sequence is not divisible by 3.")

    chunks = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    valid_chunks = []

    for chunk in chunks:
        if set(chunk).issubset({'A', 'C', 'G', 'T'}):
            valid_chunks.append(chunk)
        else:
            print(f"Warning: Chunk {chunk} contains invalid characters and will be ignored.")

    print("Codons" + "\n" + str(valid_chunks) + "\n")

    counts = {'first': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
              'second': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
              'third': {'A': 0, 'T': 0, 'C': 0, 'G': 0}}

    for chunk in valid_chunks:
        if len(chunk) >= 1 and chunk[0] in counts['first']:
            counts['first'][chunk[0]] += 1
        if len(chunk) >= 2 and chunk[1] in counts['second']:
            counts['second'][chunk[1]] += 1
        if len(chunk) >= 3 and chunk[2] in counts['third']:
            counts['third'][chunk[2]] += 1

    num_chunks = len(valid_chunks)
    for position in counts:
        for nucleotide in counts[position]:
            counts[position][nucleotide] /= num_chunks

    print("Nucleotide counts \n" + str(counts) + "\n")

    average_A = (counts['first']['A'] + counts['second']['A'] + counts['third']['A']) / 3
    average_T = (counts['first']['T'] + counts['second']['T'] + counts['third']['T']) / 3
    average_C = (counts['first']['C'] + counts['second']['C'] + counts['third']['C']) / 3
    average_G = (counts['first']['G'] + counts['second']['G'] + counts['third']['G']) / 3

    return {'A': average_A, 'T': average_T, 'C': average_C, 'G': average_G}

def main(dna_file):
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

    fieldnames = ['isolate'] + list(results.keys())
    rows = ['A', 'T', 'C', 'G']

    with open('dna_counts.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['isolate'] + [f"{identifier}_{isolate}" for identifier in results for isolate in rows])
        writer.writerow([dna_file.split('.fasta')[0]] + [results[identifier][isolate] for identifier in results for isolate in rows])

    print(f"Results saved to dna_counts.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a DNA sequence from a FASTA file.')
    parser.add_argument('dna_file', type=str, help='The path to the DNA sequence file in FASTA format')
    args = parser.parse_args()
    main(args.dna_file)
