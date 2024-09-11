'''
import pandas as pd
import joblib
import logging
from pathlib import Path
import yaml

module_dir = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def load_data_for_prediction(csv_file):
    df = pd.read_csv(csv_file)
    
    # Ensure Fingerprint is split into individual bits
    df['Fingerprint'] = df['Fingerprint'].apply(lambda x: [int(bit) for bit in x])
    
    # Expand the Fingerprint column into individual columns
    fingerprint_df = pd.DataFrame(df['Fingerprint'].tolist())
    df = pd.concat([df.drop(columns=['Fingerprint']), fingerprint_df], axis=1)
    
    return df

def make_predictions(model, df):
    # Prepare the input data by dropping non-feature columns
    X = df.drop(columns=['Name', 'Docking Score'])
    
    # Make predictions
    predictions = model.predict(X)
    
    # Combine the results with the original data
    df['Predicted Docking Score'] = predictions
    
    return df

if __name__ == "__main__":
    config_path = module_dir/'../../configure.yaml'  # Path to your configuration file
    config = read_config(config_path)

    # Load the model
    model_path = module_dir/'../training/random_forest/random_forest_model.joblib'
    logging.info(f"Loading model from {model_path}...")
    model = joblib.load(model_path)
    logging.info("Model loaded successfully.")

    # Load the data from the CSV file
    df = load_data_for_prediction(config['operations']['add_fingerprints']['csv_file'])
    
    # Make predictions using the loaded model
    df_with_predictions = make_predictions(model, df)
    
    # Save the predictions to a new CSV file
    output_file = module_dir / 'predictions.csv'
    df_with_predictions.to_csv(output_file, index=False)
    logging.info(f"Predictions saved to {output_file}")
'''
import pandas as pd
import joblib
import logging
from pathlib import Path
import yaml

module_dir = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def load_data_for_prediction(csv_file):
    df = pd.read_csv(csv_file)
    
    # Ensure Fingerprint is split into individual bits
    df['Fingerprint'] = df['Fingerprint'].apply(lambda x: [int(bit) for bit in x])
    
    # Expand the Fingerprint column into individual columns
    fingerprint_df = pd.DataFrame(df['Fingerprint'].tolist())
    df = pd.concat([df.drop(columns=['Fingerprint']), fingerprint_df], axis=1)
    
    return df

def make_predictions(model, df):
    '''
    # Prepare the input data by dropping non-feature columns
    X = df.drop(columns=['Name', 'Docking Score'])
    '''
    # Use only the 'Fingerprint' column for predictions
    X = df[['Fingerprint']]
    # Make predictions
    predictions = model.predict(X)
    
    # Combine the results with the original data
    df['Predicted Docking Score'] = predictions
    
    return df

if __name__ == "__main__":
    config_path = module_dir/'../../../configure.yaml'  # Path to your configuration file
    config = read_config(config_path)

    # Load the model
    model_path = module_dir/'../../training/random_forest/random_forest_model.joblib'
    logging.info(f"Loading model from {model_path}...")
    model = joblib.load(model_path)
    logging.info("Model loaded successfully.")

    # Load the data from the CSV file
    df = load_data_for_prediction(config['operations']['extract_scores']['output_csv'])
    
    # Make predictions using the loaded model
    df_with_predictions = make_predictions(model, df)
    
    # Select only the required columns and save the predictions to a new CSV file
    output_file = module_dir / 'predictions1.csv'
    df_with_predictions[['Name', 'Docking Score', 'Predicted Docking Score']].to_csv(output_file, index=False)
    logging.info(f"Predictions saved to {output_file}")
