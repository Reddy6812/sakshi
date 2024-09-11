import torch
from torch.utils.data import DataLoader, Dataset
from transformers import RobertaTokenizer, RobertaForSequenceClassification, Trainer, TrainingArguments
import pandas as pd
import logging
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import yaml
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MoleculeDataset(Dataset):
    def __init__(self, encodings, labels):
        self.encodings = encodings
        self.labels = labels

    def __getitem__(self, idx):
        item = {key: torch.tensor(val[idx]) for key, val in self.encodings.items()}
        item['labels'] = torch.tensor(self.labels[idx], dtype=torch.float32)
        return item

    def __len__(self):
        return len(self.labels)

def compute_metrics(pred):
    labels = pred.label_ids
    preds = pred.predictions.squeeze()

    mse = mean_squared_error(labels, preds)
    r2 = r2_score(labels, preds)
    
    # Calculate accuracy
    rounded_preds = (preds > 0.5).astype(int)  # Adjust threshold as needed for classification
    rounded_labels = (labels > 0.5).astype(int)
    accuracy = (rounded_preds == rounded_labels).mean()
    
    return {"mse": mse, "r2": r2, "accuracy": accuracy}

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)
        
def load_data_from_smiles(config):
    smiles_dir = Path(config['operations']['parse_pdbqt']['output_dir'])
    smiles_files = list(smiles_dir.glob("*.smi"))
    smiles = []
    labels = []

    # Load the docking score CSV file specified in the configuration
    docking_score_file = Path(config['operations']['extract_scores']['output_csv'])
    df = pd.read_csv(docking_score_file)

    # Load SMILES and corresponding docking scores
    for smiles_file in smiles_files:
        with open(smiles_file, "r") as f:
            content = f.read().strip()
            
            # Split the SMILES string and the name (assuming they are separated by a tab or space)
            smiles_str, name = content.split()
            
            # Extract the docking score for the current molecule
            docking_score_row = df.loc[df['Name'] == name, 'Docking Score']

            if docking_score_row.empty:
                #logging.warning(f"Docking score not found for {name}. Skipping this entry.")
                continue

            docking_score = docking_score_row.values[0]
            #logging.info(f"Read SMILES: {smiles_str}, Name: {name}, Docking Score: {docking_score}")
            
            # Only add the SMILES and docking score if both are available
            smiles.append(smiles_str)
            labels.append(docking_score)

    logging.info(f"Total SMILES processed: {len(smiles)}")
    logging.info(f"Total labels collected: {len(labels)}")
    
    return smiles, labels

def main():
    # Load configuration
    config_path = "configure.yaml"
    config = read_config(config_path)

    # Load data
    smiles, labels = load_data_from_smiles(config)

    # Train-test split
    train_smiles, test_smiles, train_labels, test_labels = train_test_split(smiles, labels, test_size=0.2, random_state=42)

    # Initialize tokenizer and model
    tokenizer = RobertaTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
    model = RobertaForSequenceClassification.from_pretrained("seyonec/ChemBERTa-zinc-base-v1", num_labels=1)

    # Encode the SMILES strings
    train_encodings = tokenizer(list(train_smiles), truncation=True, padding=True, max_length=512)
    test_encodings = tokenizer(list(test_smiles), truncation=True, padding=True, max_length=512)

    # Create datasets
    train_dataset = MoleculeDataset(train_encodings, train_labels)
    test_dataset = MoleculeDataset(test_encodings, test_labels)

    # Check if GPU is available
    device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    model.to(device)

    # Training arguments
    training_args = TrainingArguments(
        output_dir='./results',          # output directory
        num_train_epochs=5,              # total number of training epochs
        per_device_train_batch_size=8,   # batch size per device during training
        per_device_eval_batch_size=8,    # batch size for evaluation
        warmup_steps=500,                # number of warmup steps for learning rate scheduler
        weight_decay=0.01,               # strength of weight decay
        logging_dir='./logs',            # directory for storing logs
        logging_steps=10,
        evaluation_strategy="epoch"
    )

    # Initialize the Trainer
    trainer = Trainer(
        model=model,                         # the instantiated Transformers model to be trained
        args=training_args,                  # training arguments, defined above
        train_dataset=train_dataset,         # training dataset
        eval_dataset=test_dataset,           # evaluation dataset
        compute_metrics=compute_metrics,     # the callback that computes metrics of interest
    )

    # Train the model
    logging.info("Starting model training...")
    trainer.train()

    # Evaluate the model
    logging.info("Evaluating model performance...")
    metrics = trainer.evaluate()
    logging.info(f"Model evaluation metrics: {metrics}")

    # Save the model using torch.save
    model_save_path = "chemberta_model.pt"
    torch.save(model.state_dict(), model_save_path)
    logging.info(f"Model saved to '{model_save_path}'")

if __name__ == "__main__":
    main()
