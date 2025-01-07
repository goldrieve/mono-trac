from ete3 import Tree, TreeStyle, NodeStyle

# Read the Newick string from the file
with open('workdir/phylogenetic_tree.tree', 'r') as file:
    newick_str = file.read().strip()

# Create a tree from the Newick string
tree = Tree(newick_str)

# Strip .fas from all leaf names
for leaf in tree.iter_leaves():
    if leaf.name.endswith('.fas'):
        leaf.name = leaf.name[:-4]
        # Define a dictionary with leaf names and their corresponding colors
        leaf_colors = {
            "submitted_isolate.fasta": "black",
            "OVI": "#F0E68C",
            "927": "#FFD1DC",
            "940": "#F0E68C",
            "280104": "#8FBC8F",
            "108AT": "#D3D3D3",
            "108BT": "#D3D3D3",
            "15BT-relapse": "#D3D3D3",
            "340AT": "#D3D3D3",
            "348BT": "#D3D3D3",
            "57AT": "#D3D3D3",
            "ABBA": "#E6E6FA",
            "AGAL-CI-78-TCH312": "#FFD1DC",
            "Alfort": "#8FBC8F",
            "American-Strain": "#8FBC8F",
            "AnTat-12-1S": "#ADD8E6",
            "AnTat-3-1": "#8FBC8F",
            "AnTat-3-3": "#8FBC8F",
            "ATCC-30019": "#8FBC8F",
            "ATCC-30023": "#8FBC8F",
            "ATCC30019": "#8FBC8F",
            "BIM-AnTat-8-1-P8": "#D3D3D3",
            "Bosendja": "#D3D3D3",
            "BoTat-1-1": "#FFA07A",
            "Canadian-Strain": "#8FBC8F",
            "Colombia": "#8FBC8F",
            "Dodola_943": "#F0E68C",
            "E28": "#FFA07A",
            "EATRO2340": "#ADD8E6",
            "EATRO3": "#ADD8E6",
            "FEO": "#E6E6FA",
            "FEO-AnTat-16-1": "#E6E6FA",
            "GMOM-ZM-83-TRPZ-317": "#FFD1DC",
            "GPAP-CI-82-KP10-29": "#FFD1DC",
            "HAMBURG": "#8FBC8F",
            "IVMt1": "#FF6347",
            "Jua": "#D3D3D3",
            "Kazakstan": "#8FBC8F",
            "Kenya": "#8FBC8F",
            "LIGO": "#E6E6FA",
            "LiTat-1-5-P9": "#D3D3D3",
            "LOGRA": "#D3D3D3",
            "MBA": "#D3D3D3",
            "MBO-NG-74-R10": "#FFD1DC",
            "MBOT-GM-77-GB2": "#FFD1DC",
            "MCAM-ET-2013-MU-01": "#8FBC8F",
            "MCAM-ET-2013-MU-02": "#8FBC8F",
            "MCAM-ET-2013-MU-04": "#8FBC8F",
            "MCAM-ET-2013-MU-05": "#8FBC8F",
            "MCAM-ET-2013-MU-09": "#8FBC8F",
            "MCAM-ET-2013-MU-14": "#FF6347",
            "MCAM-ET-2013-MU-17": "#8FBC8F",
            "MCAP-CI-91-BALEA-2": "#FFD1DC",
            "Merzouga-56": "#8FBC8F",
            "MHOM-SD-82-BIYAMINA": "#ADD8E6",
            "MHOM-SD-82-MUSIKIA-cloneA": "#D3D3D3",
            "MHOM-ZM-80-TRPZ-23": "#ADD8E6",
            "MHOM-ZM-83-TRPZ-349": "#ADD8E6",
            "MSUS-CI-78-TSW-157": "#E6E6FA",
            "MSUS-CI-78-TSW382": "#FFD1DC",
            "MSUS-CI-82-TSW62": "#FFD1DC",
            "MSUS-CI-83-TSW-11": "#FFD1DC",
            "MU09": "#8FBC8F",
            "MU10": "#FF6347",
            "Nabe": "#D3D3D3",
            "NDMI": "#D3D3D3",
            "NKOUA": "#D3D3D3",
            "OUSOU": "#D3D3D3",
            "Pakwa": "#D3D3D3",
            "Philippines": "#8FBC8F",
            "RoTat_1.2": "#8FBC8F",
            "ROUPO-VAVOUA--80-MURAZ-14": "#D3D3D3",
            "Rumphi": "#ADD8E6",
            "STIB-816": "#8FBC8F",
            "STIB-851": "#ADD8E6",
            "STIB247": "#FFD1DC",
            "STIB386": "#D3D3D3",
            "STIB805": "#8FBC8F",
            "STIB818": "#8FBC8F",
            "STIB900": "#ADD8E6",
            "SVP": "#8FBC8F",
            "Te-Ap-N-D1": "#F0E68C",
            "EATRO_1125_AnTat_1.1_90:13": "#FFD1DC",
            "Zagora-I-17": "#8FBC8F",
        }

# Customize the tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"  # Circular tree layout
ts.root_opening_factor = 1.0  # Make the tree unrooted

# Apply colors to the leaves
for leaf in tree.iter_leaves():
    nstyle = NodeStyle()
    if leaf.name in leaf_colors:
        nstyle["bgcolor"] = leaf_colors[leaf.name]
    leaf.set_style(nstyle)

# Render the tree to a file
import os

output_file = os.path.expanduser('~/Desktop/tree.png')
tree.render(output_file, tree_style=ts)
