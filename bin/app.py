import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from geopy.geocoders import Nominatim

@st.cache_data
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
    
    # Plot histograms for each nucleotide count
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
        uploaded_file_name = os.path.splitext(uploaded_file.name)[0]
        output_dir = f'{uploaded_file_name}_output'
        os.makedirs(output_dir, exist_ok=True)
        histogram_path = os.path.join(output_dir, f'{column}_histogram.png')
        plt.savefig(histogram_path)
        histogram_paths[column] = histogram_path
        plt.close()
    
    return summary_stats, histogram_paths, combined_data_path

# Set the title of the app
st.title('MONO-TRAC Report')

# Sidebar for navigation
page = st.sidebar.selectbox("Choose a page", ["Upload CSV", "Summary Statistics", "Nucleotide Counts", "Interactive Map"])

# File uploader
uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # Process the uploaded file
    summary_stats, histogram_paths, combined_data_path = process_data(uploaded_file)
    
    if page == "Upload CSV":
        st.subheader('Upload CSV')
        st.write("File uploaded successfully.")
    
    elif page == "Summary Statistics":
        # Display summary statistics
        st.subheader('Summary Statistics')
        st.write(summary_stats)
    
    elif page == "Nucleotide Counts":
        # Display histograms in a dropdown tab
        st.subheader('Nucleotide Counts')
        selected_column = st.selectbox('Select a column to view its histogram.\nRed line indicates the new sample.', list(histogram_paths.keys()))
        st.image(histogram_paths[selected_column])
    
    elif page == "Interactive Map":
        st.subheader('Interactive World Map')
        
else:
    st.write("Please upload a CSV file to proceed.")

if page == "Interactive Map":
    st.subheader('Isolate Locations')
    locations = read.csv('ml/meta_data.csv')
    
    geolocator = Nominatim(user_agent="myGeopyApp")
    location = geolocator.geocode("England")
    country_data = pd.DataFrame(demo_data)
    
    # Display the map with dots for each country
    st.map(country_data)
