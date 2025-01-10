# Mono-Trac

Mono-Trac is a Streamlit-based web application for analyzing DNA sequences, predicting isolate types using an XGBoost model, and visualizing phylogenetic trees. The app provides various functionalities including nucleotide counts, summary statistics, and phylogenetic tree generation.

## Features

- **FASTA File Upload**: Upload a FASTA file containing DNA sequences.
- **Nucleotide Counts**: View histograms of nucleotide counts.
- **Prediction**: Predict isolate types using an XGBoost model.
- **Verbose Working**: Detailed analysis of DNA sequences.
- **Phylogenetic Tree**: Generate and visualize phylogenetic trees.
- **Tool Versions**: Display versions of the tools used in the app.

## Installation

1. **Clone the repository**:
   ```sh
   git clone https://github.com/yourusername/mono-trac.git
   cd mono-trac
   
2. **Create a virtual environment:

3. ** Install the dependencies:

4. ** Install additional tools:

MAFFT: Multiple sequence alignment tool
FastTree: Efficiently estimate maximum-likelihood phylogenetic trees
Usage
Run the Streamlit app:

Navigate to the app in your web browser:

Upload a FASTA file:

Choose a FASTA file from your local machine.
Click "Upload" to submit the file.
Navigate through the app:

Submit Country: Highlight a country on the map.
Prediction: View the prediction results for the uploaded isolate.
Verbose Working: Detailed analysis of the DNA sequences.
Nucleotide Counts: View histograms of nucleotide counts.
Phylogenetic Tree: Generate and visualize a phylogenetic tree.
Tool Versions: Display versions of the tools used in the app.
File Structure
app_full.py: Main Streamlit app script.
xgboost_model.pkl: Pre-trained XGBoost model.
meta_data.csv: Metadata for the DNA sequences.
combined.csv: Combined nucleotide counts data.
requirements.txt: List of Python dependencies.
Dependencies
Python 3.7+
Streamlit
Pandas
Matplotlib
Numpy
PyDeck
Joblib
XGBoost
Biopython
ETE3
Additional Tools
MAFFT
FastTree
License
This project is licensed under the MIT License. See the LICENSE file for details.

Acknowledgements
Streamlit
XGBoost
Biopython
ETE Toolkit
Contact
For any questions or feedback, please contact yourname@example.com.