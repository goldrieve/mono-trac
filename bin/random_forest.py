import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.ensemble import RandomForestClassifier
import joblib
import os

# Specify the input directory containing the CSV files
input_dir = 'ml/ml_data/'

# Specify the output file path
output_file = 'ml/ml_data/combined/combined.csv'

# List to hold individual DataFrames
dataframes = []

# Iterate over all files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.csv'):
        file_path = os.path.join(input_dir, filename)
        df = pd.read_csv(file_path)
        dataframes.append(df)

# Concatenate all DataFrames
combined_df = pd.concat(dataframes, ignore_index=True)

# Save the combined DataFrame to a new CSV file
combined_df.to_csv(output_file, index=False)

meta_data = pd.read_csv('ml/meta_data.csv')

# Merge combined_df with meta_data on the column 'isolate'
df = pd.merge(meta_data, combined_df, on='isolate')

df = df.drop(columns=['Country'])

# Display the merged DataFrame
df.head()

data = df

# Perform basic QC
print("Data Head:")
print(data.head())
print("\nData Info:")
print(data.info())
print("\nData Description:")
print(data.describe())

# Assume the target column is named 'target' and the rest are features
X = data.drop(columns=['isolate', 'Subspecies', 'Type', 'Subgenus', 'Competence'])
y = data['Competence']

# Convert target column to numerical values
y = y.map({'Pleomorphic': 1, 'Monomorphic': 0})

# Handle missing values by filling them with the mean of each column
X = X.fillna(0)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize the RandomForestClassifier
rf = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
rf.fit(X_train, y_train)

# Save the model
joblib.dump(rf, 'ml/model/random_forest_model.pkl')

# Evaluate the model
y_pred = rf.predict(X_test)

print("\nAccuracy Score:")
print(accuracy_score(y_test, y_pred))
print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred))
