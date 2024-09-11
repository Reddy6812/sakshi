import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, confusion_matrix, accuracy_score
import joblib
import logging
from pathlib import Path
import numpy as np
import yaml
module_dir = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def load_data(csv_file):
    df = pd.read_csv(csv_file)
    
    # Ensure Fingerprint is split into individual bits
    df['Fingerprint'] = df['Fingerprint'].apply(lambda x: [int(bit) for bit in x])
    
    # Expand the Fingerprint column into individual columns
    fingerprint_df = pd.DataFrame(df['Fingerprint'].tolist())
    df = pd.concat([df.drop(columns=['Fingerprint']), fingerprint_df], axis=1)
    
    return df

def train_and_evaluate_model(df):
    X = df.drop(columns=['Name', 'Docking Score'])
    y = df['Docking Score']
    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # Train the Random Forest model
    logging.info("Training Random Forest model...")
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred = model.predict(X_test)
    
    # Calculate the accuracy and R^2 score
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    # For regression, accuracy might not be applicable, so we'll print R^2 and MSE
    logging.info(f"Model evaluation complete. MSE: {mse}, R2: {r2}")
    
    # Save the model
    joblib.dump(model, "random_forest_model.joblib")
    logging.info("Random Forest model saved as 'random_forest_model.joblib'")
    
    '''
    # To calculate a confusion matrix, let's first round the predictions
    y_pred_rounded = np.round(y_pred)
    y_test_rounded = np.round(y_test)

    # Confusion matrix
    cm = confusion_matrix(y_test_rounded, y_pred_rounded)
    accuracy = accuracy_score(y_test_rounded, y_pred_rounded)
    
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    accuracy = accuracy_score(y_test, y_pred)
    
    
    logging.info(f"Confusion Matrix:\n{cm}")
    logging.info(f"Accuracy: {accuracy}")
    '''
    
if __name__ == "__main__":
    config_path = module_dir/'../../../configure.yaml'  # Path to your configuration file
    config = read_config(config_path)

    # Load the data from the CSV file
    df = load_data(config['task']['ml']['random_forest']['trainingcsv'])
    
    # Train and evaluate the model
    train_and_evaluate_model(df)
