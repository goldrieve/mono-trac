import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import joblib
import os
import xgboost as xgb

# Specify the input directory containing the CSV files
input_dir = 'ml/ml_data'

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

# Assume the target column is named 'target' and the rest are features
X = data.drop(columns=['isolate', 'Subspecies', 'Type', 'Subgenus', 'Competence', 'colour', 'Tb927.5.1220_A', 'Tb927.5.1220_T', 'Tb927.5.1220_C', 'Tb927.5.1220_G'])

# Perform basic QC
print("Data Head:")
print(X.head())
print("\nData Info:")
print(X.info())
print("\nData Description:")
print(X.describe())

y = data['Competence']

# Convert target column to numerical values
y = y.map({'Pleomorphic': 1, 'Monomorphic': 0})

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)
dtest = xgb.DMatrix(X_test, label=y_test)

# Set up parameters for xgboost
params = {
    'objective': 'binary:logistic',  # For binary classification
}

bst = xgb.train(params, dtrain, num_boost_round=100, evals=[(dtest, 'test')], early_stopping_rounds=10)

# Save the model
bst.save_model('ml/model/xgboost_model.json')

# Evaluate the model
y_pred = bst.predict(dtest)
y_pred_binary = [1 if pred > 0.5 else 0 for pred in y_pred]

print("\nAccuracy Score:")
print(accuracy_score(y_test, y_pred_binary))
print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred_binary))
print("\nClassification Report:")
print(classification_report(y_test, y_pred_binary))

# Save the model using joblib
joblib.dump(bst, 'ml/model/xgboost_model.pkl')