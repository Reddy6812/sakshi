import pandas as pd

def add_smiles_to_results(results_csv, input_smi, output_csv):
    try:
        # Read the results.csv file
        results_df = pd.read_csv(results_csv)

        # Read the onefile.smi file
        smiles_df = pd.read_csv(input_smi, sep='\t', header=None, names=['SMILES', 'Name'])

        # Merge the DataFrames on the 'Name' column
        merged_df = pd.merge(results_df, smiles_df, on='Name', how='left')

        # Save the updated DataFrame to a new CSV file
        merged_df.to_csv(output_csv, index=False)
        print(f"Updated data written to {output_csv}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

def add_variant_column(output_csv, output2gan):
    try:
        # Read the updated_results.csv file
        df = pd.read_csv(output_csv)

        # Split the 'Name' column into 'Base Name' and 'Variant' columns
        df[['Base Name', 'Variant']] = df['Name'].str.split('-', expand=True)
        
        # Convert the 'Variant' column to an integer type if needed
        df['Variant'] = df['Variant'].astype(int)

        # Save the updated DataFrame to a new CSV file
        df.to_csv(output2gan, index=False)
        print(f"Processed data written to {output_csv}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
        
if __name__ == "__main__":
    results_csv = "/blue/yanjun.li/vi.gade1/seabra-li/data/fingerprints_docking_scores.csv"  # Path to results.csv
    input_smi = "/blue/yanjun.li/vi.gade1/seabra-li/data/smiles1/onefile.smi"  # Path to onefile.smi
    output_csv = "/blue/yanjun.li/vi.gade1/seabra-li/data/updated_results.csv"  # Path to save updated results
    output2gan = "/blue/yanjun.li/vi.gade1/seabra-li/data/input2gan.csv"
    add_smiles_to_results(results_csv, input_smi, output_csv)
    add_variant_column(output_csv, output2gan)
    