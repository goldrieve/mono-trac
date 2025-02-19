import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import pydeck as pdk
import csv
import joblib
import re
import xgboost as xgb
import subprocess
import matplotlib
import Bio
import ete3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree, TreeStyle, NodeStyle

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
        'Tb927.6.1100', 'Tb927.11.3980', 'Tb927.6.2630', 'Tb927.7.6560',
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
    sequence = sequence.replace('\n', '').upper()

    if len(sequence) % 3 != 0:
        st.warning("Warning: The length of the DNA sequence is not divisible by 3. Some codons may be incomplete and will be ignored.")

    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    
    valid_codons = []
    for codon in codons:
        if all(nucleotide in 'ATCG' for nucleotide in codon):
            valid_codons.append(codon)

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

    average_A = round((counts['first']['A'] + counts['second']['A'] + counts['third']['A']) / 3, 5)
    average_T = round((counts['first']['T'] + counts['second']['T'] + counts['third']['T']) / 3, 5)
    average_C = round((counts['first']['C'] + counts['second']['C'] + counts['third']['C']) / 3, 5)
    average_G = round((counts['first']['G'] + counts['second']['G'] + counts['third']['G']) / 3, 5)

    return {'A': average_A, 'T': average_T, 'C': average_C, 'G': average_G}

@st.cache_data
def process_fasta(uploaded_file):
    file_path = f"/tmp/{uploaded_file.name}"
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    clean_fasta_file(file_path)
    
    with open(file_path, 'r') as file:
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
    output_file = f"/tmp/{os.path.basename(file_path).split('.')[0]}_counts.csv"

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['isolate'] + [f"{identifier}_{isolate}" for identifier in results for isolate in rows])
        writer.writerow([os.path.basename(file_path).split('.')[0]] + [round(results[identifier][isolate], 5) for identifier in results for isolate in rows])

    mod = joblib.load('ml/model/xgboost_model.pkl')
    new_data = pd.read_csv(output_file)
    new_data['lat'] = 0
    new_data['lon'] = 0
    columns = ['lat', 'lon'] + [col for col in new_data.columns if col not in ['lat', 'lon']]
    new_data = new_data[columns]
    new_data_features = new_data.drop(columns=['isolate'])
    dmatrix_data = xgb.DMatrix(new_data_features)
    predictions = mod.predict(dmatrix_data)
    predictions_binary = ['Pleomorphic' if pred > 0.5 else 'Monomorphic' for pred in predictions]
    raw_scores = mod.predict(dmatrix_data, output_margin=True)
    proba = 1 / (1 + np.exp(-raw_scores))
    proba_dict = {'Pleomorphic': round(proba[0], 3), 'Monomorphic': round(1 - proba[0], 3)}
    proba_str = f"Pleomorphic: {proba_dict['Pleomorphic']:.3f}, Monomorphic: {proba_dict['Monomorphic']:.3f}"
    return results, predictions_binary[0], proba_str, output_file

def process_data(uploaded_file):
    # Load data
    data = pd.read_csv(uploaded_file)
    
    # Path to the files to concatenate with
    combined_counts_path = os.path.join('ml/ml_data/combined/combined.csv')
    meta_data_path = os.path.join('ml/meta_data.csv')
    
    # Load the combined counts data
    combined_counts_data = pd.read_csv(combined_counts_path)
    
    # Load the meta data
    meta_data = pd.read_csv(meta_data_path)
    
    # Merge combined counts data with meta data on the 'isolate' column
    combined_counts_data = pd.merge(combined_counts_data, meta_data, on='isolate', how='inner')
    
    # Concatenate the uploaded file data with the combined counts data
    concat_data = pd.concat([data, combined_counts_data], ignore_index=True).fillna('?')
    
    # Save the combined data to a CSV file
    combined_data_path = 'combined_data.csv'
    concat_data.to_csv(combined_data_path, index=False)
    
    # Generate summary statistics for nucleotide counts
    summary_stats = concat_data.describe()
    
    # Plot histograms for each numeric column
    numeric_columns = concat_data.select_dtypes(include=[np.float64]).columns
    histogram_paths = {}
    for column in numeric_columns:
        plt.figure()
        competence_values = concat_data['Competence'].unique()
        for competence in competence_values:
            subset = concat_data[concat_data['Competence'] == competence]
            if not subset.empty and len(subset[column]) > 1:
                try:
                    data_values = subset[column].dropna().values
                    if len(data_values) > 1:
                        # Plot histogram
                        plt.hist(data_values, bins=30, alpha=0.5, label=f'{competence}')
                except Exception as e:
                    st.warning(f"Could not create histogram for {column} - {competence}: {str(e)}")
                    continue
        
        # Add line for the uploaded file values
        if column in data.columns:
            uploaded_values = data[column].dropna().values
            if len(uploaded_values) > 0:
                plt.axvline(uploaded_values, color='r', linestyle='dashed', linewidth=1)
        
        plt.title(f'{column}')
        plt.xlabel(column)
        plt.ylabel('Frequency')
        plt.legend()
        
        # Save the histogram
        uploaded_file_name = os.path.splitext(os.path.basename(uploaded_file))[0]
        output_dir = f'{uploaded_file_name}_output'
        os.makedirs(output_dir, exist_ok=True)
        histogram_path = os.path.join(output_dir, f'{column}_histogram.png')
        plt.savefig(histogram_path)
        histogram_paths[column] = histogram_path
        plt.close()
    
    return summary_stats, histogram_paths, combined_data_path

def set_page(page_name):
    st.session_state.page = page_name

# Set the title of the app
st.title('mono-trac')

# Sidebar for navigation
if 'results' not in st.session_state:
    page = st.sidebar.selectbox("", [""])
else:
    page = st.sidebar.selectbox("Results", ["Submit Country", "Prediction", "Summary Statistics", "Nucleotide Counts", "Verbose ML Working", "Phylogenetic Tree", "Tool Versions"])

# File uploader
if 'results' not in st.session_state:
    uploaded_file = st.file_uploader("Choose a FASTA file", type=["fasta", ".fa", ".fas"])
    if uploaded_file is not None:
        st.session_state.uploaded_file = uploaded_file

    if uploaded_file is not None:
        results, prediction, proba, output_file = process_fasta(uploaded_file)
        st.session_state.results = results
        st.session_state.prediction = prediction
        st.session_state.proba = proba
        st.session_state.output_file = output_file
        st.success("FASTA file submission successful. Click here to proceed.", icon="✅")
        st.button("Proceed to Submit Country", on_click=set_page, args=("Submit Country",))

if 'results' in st.session_state:
    results = st.session_state.results
    prediction = st.session_state.prediction
    proba = st.session_state.proba
    output_file = st.session_state.output_file
    
    if page == "Prediction":
        st.subheader(f"{prediction.upper()}")
        st.write(f"Based on an XGBoost model, built using nucleotide usage for genes known to mutate in previous monomorphic outreaks, we predict your isolate is {prediction}.")
        st.write(f"Probability: {proba}")
        st.write(f"Please continue to the following tabs for a more detailed breakdown of the data.")

    elif page == "Summary Statistics":
        # Display summary statistics
        st.subheader('Summary Statistics')
        if 'summary_stats' not in st.session_state:
            summary_stats, histogram_paths, combined_data_path = process_data(output_file)
            st.session_state.summary_stats = summary_stats
            st.session_state.histogram_paths = histogram_paths
            st.session_state.combined_data_path = combined_data_path
        
        summary_stats = st.session_state.summary_stats
        histogram_paths = st.session_state.histogram_paths
        combined_data_path = st.session_state.combined_data_path
        
        st.write(summary_stats)
        numeric_columns = pd.read_csv(combined_data_path).select_dtypes(include=[np.number]).columns
        st.line_chart(pd.read_csv(combined_data_path).set_index('isolate')[numeric_columns])
        
    elif page == "Nucleotide Counts":
        # Display histograms in a dropdown tab
        st.subheader('Nucleotide Counts')
        if 'summary_stats' not in st.session_state:
            summary_stats, histogram_paths, combined_data_path = process_data(output_file)
            st.session_state.summary_stats = summary_stats
            st.session_state.histogram_paths = histogram_paths
            st.session_state.combined_data_path = combined_data_path
        
        summary_stats = st.session_state.summary_stats
        histogram_paths = st.session_state.histogram_paths
        combined_data_path = st.session_state.combined_data_path
        
        selected_column = st.selectbox('Select a column to view its histogram.\nRed line indicates the new sample.', list(histogram_paths.keys()))
        st.image(histogram_paths[selected_column])
       
    elif page == "Verbose ML Working":
        st.subheader('Verbose XGBoost Working')
        if 'results' in st.session_state:
            results = st.session_state.results
            output_file = st.session_state.output_file
            
            for identifier, summary in results.items():
                st.text(f"\nAnalysing {identifier}\n")
                st.text(f"Nucleotide counts\n{summary}\n")
            st.text(f"Prediction: {prediction}")
            st.text(f"Results saved to {output_file}\n")
            
            st.subheader('XGBoost Model')
            try:
                with open('ml/model/xgb.log', 'r') as file:
                    xgb_log_contents = file.read()
                    st.code(xgb_log_contents, language='text')
            except FileNotFoundError:
                st.warning("XGBoost log file not found.")
        else:
            st.warning("No results available. Please upload a FASTA file first.")

    elif page == "Phylogenetic Tree":
        st.subheader('Phylogenetic Tree')
        
        if 'results' in st.session_state:
            uploaded_file = st.session_state.uploaded_file

            # Define the directory containing the FASTA files
            fasta_dir = 'isolate_fasta/'
            os.makedirs(fasta_dir, exist_ok=True)
            uploaded_file_path = os.path.join(fasta_dir, 'submitted_isolate.fasta')
            with open(uploaded_file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            output_dir = 'workdir'
            os.makedirs(output_dir, exist_ok=True)

            concatenated_fasta = os.path.join(output_dir, 'concatenated_sequences.fasta')
            with open(concatenated_fasta, 'w') as outfile:
                for fasta_file in os.listdir(fasta_dir):
                    if fasta_file.endswith('.fas') or fasta_file.endswith('.fasta'):
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

            aligned_fasta = os.path.join(output_dir, 'aligned_sequences.fasta')
            subprocess.run(['mafft', '--anysymbol', '--auto', concatenated_fasta], stdout=open(aligned_fasta, 'w'))

            tree_file = os.path.join(output_dir, 'phylogenetic_tree.tree')
            subprocess.run(['fasttree', '-nt', aligned_fasta], stdout=open(tree_file, 'w'))

            with open('workdir/phylogenetic_tree.tree', 'r') as file:
                newick_str = file.read().strip()
           
            tree = Tree(newick_str)

            for leaf in tree.iter_leaves():
                if leaf.name.endswith('.fas'):
                    leaf.name = leaf.name[:-4]
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

            output_file = os.path.expanduser('~/Desktop/tree.png')
            tree.render(output_file, tree_style=ts)
            st.write(f"Submitted sample = black")
            st.write(f"T.b.brucei = pink, T.b.gambiense = grey, T.b. gambiense II = purple, T.b.rhodesiense = blue, T.b.evansi type A = green, T.b.evansi type B = red, T.b.equiperdum type OVI = gold, T.b.equiperdum type BoTat = orange")
            st.image(output_file, caption='Phylogenetic Tree')
        
if page == "Submit Country":
    st.subheader('Legend \n Blue = pleomorphic, Red = monomorphic, Green = Submitted Country')
    meta_data = pd.read_csv('ml/meta_data.csv')
    
    if 'Competence' in meta_data.columns:
        # Create a fixed color map for competence levels
        competence_levels = meta_data['Competence'].unique()
        color_map = {level: [0, 0, 255] if i % 2 == 0 else [255, 0, 0] for i, level in enumerate(competence_levels)}
        
        # Add color column to the country data
        country_data = meta_data[['Country', 'lat', 'lon', 'Competence', 'isolate', 'Type']].dropna()
        country_data['color'] = country_data['Competence'].map(color_map)
        
        # Jitter the points if there are multiple points at the same lat and lon
        jitter_amount = 1
        country_data['lat'] += np.random.uniform(-jitter_amount, jitter_amount, size=len(country_data))
        country_data['lon'] += np.random.uniform(-jitter_amount, jitter_amount, size=len(country_data))
        
        # Add the submitted country with a unique color
        selected_country = st.sidebar.text_input("Enter the country name to highlight on the map", key="highlight_country")
        if selected_country:
            countries_df = pd.read_csv('bin/countries.csv')
            submitted_country_data = countries_df[countries_df['Country'] == selected_country]
            if not submitted_country_data.empty:
                submitted_country_data = submitted_country_data[['Country', 'lat', 'lon']].dropna()
                submitted_country_data['Competence'] = 'Submitted Country'
                submitted_country_data['color'] = [[0, 255, 0]]  # Green for the submitted country
                submitted_country_data['isolate'] = 'Submitted Country'
                country_data = pd.concat([country_data, submitted_country_data], ignore_index=True)
            if not submitted_country_data.empty:
                st.success("Country submitted. Click on the 'Results' menu on the left.", icon="✅")
        
        # Create a PyDeck layer for the map
        layer = pdk.Layer(
            'ScatterplotLayer',
            data=country_data,
            get_position='[lon, lat]',
            get_color='color',
            get_radius=50000,
            pickable=True,
            radius_scale=10,  
            radius_min_pixels=1,  
            radius_max_pixels=5,
            get_tooltip='isolate'
        )
        
        # Set the viewport location
        view_state = pdk.ViewState(
            latitude=country_data['lat'].mean(),
            longitude=country_data['lon'].mean(),
            zoom=1
        )
        
        # Render the map
        st.pydeck_chart(pdk.Deck(layers=[layer], initial_view_state=view_state, tooltip={"text": "{isolate} \n {Type}"}))
    else:
        st.warning("The meta_data.csv file does not contain the required column: 'Competence'.")
    

def get_tool_versions():
            tools = {
                "Python": subprocess.run(["python", "--version"], capture_output=True, text=True).stdout.strip(),
                "Streamlit": st.__version__,
                "Pandas": pd.__version__,
                "Matplotlib": matplotlib.__version__,
                "Numpy": np.__version__,
                "PyDeck": pdk.__version__,
                "Joblib": joblib.__version__,
                "XGBoost": xgb.__version__,
                "Biopython": Bio.__version__,
                "ETE3": ete3.__version__,
                "Mafft": subprocess.run("mafft --version", capture_output=True, text=True, shell=True).stdout.strip().split('\n')[0].split(' ')[0],
                "FastTree": subprocess.run("fasttree -version", capture_output=True, text=True, shell=True).stdout.strip().split('\n')[0].split(' ')[0]
            }
            return tools

if page == "Tool Versions":
            st.subheader('Tool Versions')
            tool_versions = get_tool_versions()
            tool_versions_df = pd.DataFrame(tool_versions.items(), columns=["Tool", "Version"])
            st.table(tool_versions_df)