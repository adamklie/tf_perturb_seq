import os
import requests
import pandas as pd
import argparse
from typing import List, Tuple

def download_fastq_files(df, access_key: str, secret_key: str):
    """
    Download fastq files for pairs of R1_path and R2_path accessions into a fastq_files directory
    
    Args:
        df: Pandas DataFrame containing R1_path and R2_path columns with accession IDs
        access_key: IGVF API access key
        secret_key: IGVF API secret key
    """
    # Create fastq_files directory if it doesn't exist
    fastq_dir = "fastq_files"
    os.makedirs(fastq_dir, exist_ok=True)
    
    base_url = "https://api.data.igvf.org/sequence-files"
    accession_pairs = list(zip(df['R1_path'], df['R2_path']))
    
    for r1, r2 in accession_pairs:
        # Download R1
        r1_url = f"{base_url}/{r1}/@@download/{r1}.fastq.gz"
        r1_response = requests.get(r1_url, auth=(access_key, secret_key))
        if r1_response.status_code == 200:
            output_path = os.path.join(fastq_dir, f"{r1}.fastq.gz")
            with open(output_path, 'wb') as f:
                f.write(r1_response.content)
            print(f"Downloaded {r1} to {fastq_dir}")
        else:
            print(f"Failed to download {r1}: {r1_response.status_code}")
        
        # Download R2
        r2_url = f"{base_url}/{r2}/@@download/{r2}.fastq.gz"
        r2_response = requests.get(r2_url, auth=(access_key, secret_key))
        if r2_response.status_code == 200:
            output_path = os.path.join(fastq_dir, f"{r2}.fastq.gz")
            with open(output_path, 'wb') as f:
                f.write(r2_response.content)
            print(f"Downloaded {r2} to {fastq_dir}")
        else:
            print(f"Failed to download {r2}: {r2_response.status_code}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Download FASTQ files from IGVF API')
    parser.add_argument('--sample', required=True, 
                    help='Path to the TSV file containing R1_path and R2_path columns')
    parser.add_argument('--access-key', required=True,
                    help='IGVF API access key')
    parser.add_argument('--secret-key', required=True,
                    help='IGVF API secret key')
    
    args = parser.parse_args()
    
    # Read the per-sammple TSV file
    try:
        df = pd.read_csv(args.sample, sep='\t')
        if 'R1_path' not in df.columns or 'R2_path' not in df.columns:
            raise ValueError("TSV file must contain 'R1_path' and 'R2_path' columns")
        
        # Download the files
        download_fastq_files(df, args.access_key, args.secret_key)
        
    except FileNotFoundError:
        print(f"Error: Could not find the file '{args.sample}'")
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{args.sample}' is empty")
    except Exception as e:
        print(f"An error occurred: {str(e)}")