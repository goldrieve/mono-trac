# Mono-Trac

Mono-Trac is (currently) a Streamlit-based web application for analysing DNA sequences from the mono-trac pannel, predicting isolate developmental competence and outbreak potential using an XGBoost model, and visualising the isolates phylogeny. The app is currently in development, we are working on writing the NextFlow pipeline to generate consenus sequences of the mono-trac genes, required as input for this app, from raw Nanopore reads

## Installation

1. **Clone the repository**:
   ```sh
   git clone https://github.com/goldriev/mono-trac.git
   cd mono-trac
   
2. **Create a virtual environment**:
   ```sh
   conda env create -f environment.yml
   conda activate mono-trac

3. **Run the app locally**:
   ```sh
   streamlit run bin/app_fully.py

## Usage

Executing the streamlit command above will load a new tab in your browser. Navigate to the app and work your way through the following pages.

1. **Upload a FASTA file:**

Choose a FASTA file from your local machine. The fasta file shoudl contain the 19 gene sequences from the mono-trac pannel, generated using the mono-trac amplicon pipeline.

Click "Upload" to submit the file.

2. **Submit Country:** 

Input a country name and visualise the geographic location in the context of previous isolates.

3. **Prediction:** 

View the outbreak potential prediction for the uploaded isolate.

4. **Verbose Working:** 

Detailed analysis of the mono-trac sequences.

5. **Nucleotide Counts:** 

View histograms of nucleotide usage counts.

6. **Phylogenetic Tree:** 

Generate and visualise a phylogenetic tree which places your isolate in the context of previous T. brucei genome data.

7. **Tool Versions:** 

Display versions of the tools used in the app.

## Contact

For any questions or feedback, please contact guy.oldrieve@ed.ac.uk

## License

This project is licensed under the GPL-3.0 license. See the LICENSE file for details.