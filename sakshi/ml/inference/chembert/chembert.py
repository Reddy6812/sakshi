import torch
from transformers import RobertaTokenizer, RobertaForSequenceClassification
import pandas as pd
import logging
from pathlib import Path
module_dir = Path(__file__).parent
import yaml

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def load_smiles(smiles_file_path):
    """
    Load SMILES strings from a specified file.
    """
    smiles_data = pd.read_csv(smiles_file_path, sep='\s+', header=None)  # Using sep='\s+' as recommended
    smiles_list = smiles_data.iloc[:, 0].tolist()  # Assuming SMILES strings are in the first column
    return smiles_list

def load_saved_model_and_tokenizer(model_path, tokenizer_name="seyonec/ChemBERTa-zinc-base-v1"):
    """
    Load the saved ChemBERTa model from the specified file and initialize the tokenizer.
    """
    tokenizer = RobertaTokenizer.from_pretrained(tokenizer_name)
    model = RobertaForSequenceClassification.from_pretrained(tokenizer_name, num_labels=1)

    # Load the model's state_dict from the .pt file
    model.load_state_dict(torch.load(model_path))
    return tokenizer, model

def preprocess_smiles(smiles_list, tokenizer):
    """
    Preprocess SMILES strings using the tokenizer to prepare them for model prediction.
    """
    encodings = tokenizer(smiles_list, truncation=True, padding=True, max_length=512, return_tensors='pt')
    return encodings

def predict_docking_scores(model, encodings, device):
    """
    Make predictions using the loaded model and preprocessed encodings.
    """
    model.eval()
    model.to(device)
    encodings = {key: val.to(device) for key, val in encodings.items()}

    with torch.no_grad():
        outputs = model(**encodings)
        predictions = outputs.logits.squeeze().cpu().numpy()

    return predictions

def main():
    # Load configuration
    config_path = "/blue/yanjun.li/vi.gade1/seabra-li/configure.yaml"
    config = read_config(config_path)

    # Load SMILES from the specified file
    smiles_file_path = config['operations']['parse_pdbqt']['output_file']
    smiles_list = load_smiles(smiles_file_path)
    print(module_dir)
    # Path to the saved model file
    model_path = "/blue/yanjun.li/vi.gade1/seabra-li/chemberta_model.pt"
    
    # Load the saved model and tokenizer
    tokenizer, model = load_saved_model_and_tokenizer(model_path)
    
    # Check if GPU is available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Preprocess the input SMILES strings
    encodings = preprocess_smiles(smiles_list, tokenizer)
    
    # Make predictions using the loaded model
    predictions = predict_docking_scores(model, encodings, device)
    
    # Output the predictions
    for smi, pred in zip(smiles_list, predictions):
        logging.info(f"SMILES: {smi}, Predicted Docking Score: {pred}")

    # Save the predictions to a CSV file
    results_df = pd.DataFrame({'SMILES': smiles_list, 'Predicted Docking Score': predictions})
    output_file = Path(config['operations']['extract_scores']['output_csv']).with_name('predicted_docking_scores.csv')
    results_df.to_csv(output_file, index=False)
    logging.info(f"Predictions saved to '{output_file}'")
    print(module_dir)

if __name__ == "__main__":
    main()
