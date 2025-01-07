from ete3 import Tree, TreeStyle

# Read the Newick string from the file
with open('workdir/phylogenetic_tree.tree', 'r') as file:
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