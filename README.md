# transcriptomics-data-query-and-retrieval
A collection of convenience functions to search, retrieve, and prepare transcriptomics data from GEO and GDC.

## Prerequisites
- Python 3.7 or higher
- Conda is optional

## Setup

Run the below commands in any terminal (Linux, Mac, Windows Powershell or CMD). Replace dummy email with your email which will be submitted with your queries to the NCBI API.
```zsh
git clone https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval.git
cd transcriptomics-data-query-and-retrieval
echo YOUR_EMAIL@EXAMPLE.COM > email_for_ncbi_tracking.txt

# If not using conda:
pip install -r requirements.txt

# If using conda, use commented commands below

# Option 1: Create new Conda environment
# conda env create -f environment.yml
# conda activate transcriptomics_data_query

# Option 2: Update existing Conda environment
# conda env update --file environment.yml --name existing_env_name
# conda activate existing_env_name
```
