import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define the directory containing the FASTA files
fasta_dir = 'isolate_fasta/'
output_dir = 'workdir'
os.makedirs(output_dir, exist_ok=True)

# Step 1: Concatenate all sequences from each file into a single sequence
concatenated_fasta = os.path.join(output_dir, 'concatenated_sequences.fasta')
with open(concatenated_fasta, 'w') as outfile:
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith('.fas'):
            file_path = os.path.join(fasta_dir, fasta_file)
            concatenated_sequence = ''
            for record in SeqIO.parse(file_path, 'fasta'):
                concatenated_sequence += str(record.seq)
            # Create a new record with the concatenated sequence
            new_record = SeqRecord(
                Seq(concatenated_sequence),
                id=fasta_file,
                description=''
            )
            SeqIO.write(new_record, outfile, 'fasta')

# Step 2: Align the sequences using MAFFT
aligned_fasta = os.path.join(output_dir, 'aligned_sequences.fasta')
subprocess.run(['mafft', '--anysymbol', '--auto', concatenated_fasta], stdout=open(aligned_fasta, 'w'))

# Step 3: Generate the phylogenetic tree using FastTree
tree_file = os.path.join(output_dir, 'phylogenetic_tree.tree')
subprocess.run(['fasttree', '-nt', aligned_fasta], stdout=open(tree_file, 'w'))

from ete3 import Tree, TreeStyle

# Read the Newick string from the file
with open(tree_file, 'r') as file:
    newick_str = file.read().strip()

# Create a tree from the Newick string
tree = Tree(newick_str)

# Customize the tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"  # Circular tree layout
ts.root_opening_factor = 1.0  # Make the tree unrooted

# Render the tree to a file
output_file = '/Users/goldriev/mono-trac/bin/tree.png'
tree.render(output_file, tree_style=ts)