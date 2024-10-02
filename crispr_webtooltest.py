import requests, sys,time,json,subprocess
import pandas as pd
from Bio.Seq import Seq


# Function to get the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Function to initiate ID mapping
def initiate_id_mapping(ids, from_db, to_db):
    url = "https://rest.uniprot.org/idmapping/run"
    response = requests.post(url, data={'ids': ids, 'from': from_db, 'to': to_db})
    response.raise_for_status()
    return response.json()["jobId"]

# Function to check ID mapping results
def check_id_mapping_results(job_id):
    url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

# Function to convert a single Ensembl ID using a GET request
def convert_single_id(ensembl_id):
    url = f"https://biotools.fr/human/ensembl_symbol_converter/?api=1&id={ensembl_id}"
    response = requests.get(url)
    return response.json()

# Function to convert multiple Ensembl IDs using a POST request
def convert_multiple_ids(ensembl_ids):
    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    ids_json = json.dumps(ensembl_ids)
    body = {'api': 1, 'ids': ids_json}
    response = requests.post(url, data=body)
    return response.json()

# Get the Gene IDs from user input
def main(gene_ids):
    print(f"Received Gene IDs: {gene_ids}")
    print("Processing CRISPR data...")

    # Progress update: 10% (after input is received)
    print("Progress: 10%",flush=True)

    job_id = initiate_id_mapping(gene_ids, "GeneCards", "UniProtKB")
    print(f"Job ID: {job_id}")

# Progress update: 20% (after initiating ID mapping)
    print("Progress: 20%",flush=True)

    while True:
        results = check_id_mapping_results(job_id)
        if "results" in results:
            break
        print("Waiting for results...")
        time.sleep(5)

    uniprot_accession_codes = [result["to"] for result in results["results"]]
    print("UniProt Accession Codes:", uniprot_accession_codes)

    if uniprot_accession_codes:
        requestURL = f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_accession_codes[0]}"
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()

    

        response_data = r.json()
        relevant_info = []
        for gene in response_data.get("gnCoordinate", []):
            ensembl_gene_id = gene.get("ensemblGeneId")
            genomic_location = gene.get("genomicLocation", {}) 
            if "exon" in genomic_location:
                for exon in genomic_location["exon"]:
                    exon_id = exon.get("id")
                    chromosome = genomic_location.get("chromosome")
                    start = exon.get("genomeLocation", {}).get("begin", {}).get("position")
                    end = exon.get("genomeLocation", {}).get("end", {}).get("position")
                    relevant_info.append({
                        "ensembl_gene_id": ensembl_gene_id,
                        "exon_id": exon_id,
                        "chromosome": chromosome,
                        "start": start,
                        "end": end
                    })



        if relevant_info:
            last_exon_info = max(relevant_info, key=lambda x: x['start'])
            print("Last exon information:")
            print(last_exon_info)
            return last_exon_info  # Returning only the last exon information
        else:
            print("No relevant exon information found for valid gene IDs.")
            return None
    else:
        print("No UniProt accession codes found.")
        return None

# Progress update: 30% (after retrieving UniProt accession codes)
    print("Progress: 30%",flush=True)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        gene_ids = sys.argv[1]
    else:
        gene_ids = input("Enter the Gene IDs (comma-separated): ")
    last_exon_info = main(gene_ids)

# Fetch the DNA sequence for the last exon
if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    start = exon_end - 500
    end = exon_end + 350
    
    if exon_end < exon_start:
        start = exon_end - 350
        end = exon_end + 500
        exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    print(f"New coordinates: Chromosome {chromosome}, Start {start}, End {end}")

    ensembl_server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?coord_system_version=GRCh38"
    
    headers = {"Content-Type": "application/json"}
    r = requests.get(ensembl_server + ext, headers=headers)


    if not r.ok:
        r.raise_for_status()
        sys.exit()
 # Progress update: 40% (after fetching protein coordinates)
        print("Progress: 40%",flush=True)

    dna_sequence = r.json().get("seq")
    
    # If the start is greater than the end, print the reverse complement of the extended DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_DNA_extended = reverse_complement(dna_sequence)
        print("Reverse complement of the extended DNA sequence:")
        gdna_sequence = reverse_DNA_extended
    else:
        print(f"DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        gdna_sequence = dna_sequence

    print(gdna_sequence)

    num_base_pairs = len(gdna_sequence)

    # Fetch the DNA sequence for the last exon itself
    exon_ext = f"/sequence/region/human/{chromosome}:{exon_start}..{exon_end}:1?coord_system_version=GRCh38"
    r = requests.get(ensembl_server + exon_ext, headers=headers)
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    exon_dna_sequence = r.json().get("seq")
    
    # If the start is greater than the end, print the reverse complement of the exon DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_exon = reverse_complement(exon_dna_sequence)
        exon_seq = Seq(reverse_complement_exon)
    else:
        exon_seq = Seq(exon_dna_sequence)

    num_exon_base_pairs = len(exon_dna_sequence)

    # Check if the length of the exon sequence is a multiple of three
    if num_exon_base_pairs % 3 != 0:
        # Calculate the number of bases to remove
        bases_to_remove = num_exon_base_pairs % 3
        
        # Remove bases from the start of the sequence
        exon_seq = exon_seq[bases_to_remove:]
    
    # Transcribe and translate the exon DNA sequence to amino acids
    amino_acid_seq = exon_seq.translate(to_stop=True)
    print(f"Amino acid sequence for the last exon ({last_exon_info['exon_id']}):")
    print(amino_acid_seq)

exon_ext = f"/sequence/region/human/{chromosome}:{exon_start}..{exon_end}:1?coord_system_version=GRCh38"
r = requests.get(ensembl_server + exon_ext, headers=headers)

if not r.ok:
    r.raise_for_status()
    sys.exit()

exon_dna_sequence = r.json().get("seq")
# If the start is greater than the end, print the reverse complement of the extended DNA sequence
if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_exon = reverse_complement(exon_dna_sequence)
        print("Reverse complement of the Exon sequence:")
        last_exon_seq = reverse_complement_exon
else:
        print(f"DNA sequence for the last exon({last_exon_info['exon_id']}):")
        last_exon_seq = exon_dna_sequence

print(last_exon_seq)

if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    upstream_start = exon_end - 800
    downstream_end = exon_end + 800
    
    if exon_end < exon_start:
        upstream_start = exon_end - 800
        downstream_end = exon_end + 800
        exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

    ensembl_server = "https://rest.ensembl.org"
    
    # Fetch the whole sequence including upstream and downstream
    whole_ext = f"/sequence/region/human/{chromosome}:{upstream_start}..{downstream_end}:1?coord_system_version=GRCh38"
    headers = {"Content-Type": "application/json"}
    r_whole = requests.get(ensembl_server + whole_ext, headers=headers)
    
    if not r_whole.ok:
        r_whole.raise_for_status()
        sys.exit()

    whole_dna_sequence_800 = r_whole.json().get("seq")
    
# If the start is greater than the end, print the reverse complement of the extended DNA sequence
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_extended = reverse_complement(whole_dna_sequence_800)
        print("Reverse complement of the extended DNA sequence +/- 800 bp:")
        exon_seq_800 = reverse_complement_extended
    else:
        print(f"DNA sequence for the region surrounding the last exon +/- 800 bp({last_exon_info['exon_id']}):")
        exon_seq_800 = whole_dna_sequence_800

print(exon_seq_800)

if last_exon_info:
    chromosome = last_exon_info['chromosome']
    exon_start = last_exon_info['start']
    exon_end = last_exon_info['end']
    upstream_start = exon_end - 50
    downstream_end = exon_end + 32
    
    if exon_end < exon_start:
        upstream_start = exon_end - 32
        downstream_end = exon_end + 50
        #exon_start, exon_end = exon_end, exon_start  # Swap if end is less than start
    
    # Print the new coordinates
    print(f"New coordinates: Chromosome {chromosome}, Upstream Start {upstream_start}, Downstream End {downstream_end}")

    ensembl_server = "https://rest.ensembl.org"
    
    # Fetch the whole sequence including upstream and downstream
    whole_ext = f"/sequence/region/human/{chromosome}:{upstream_start}..{downstream_end}:1?coord_system_version=GRCh38"
    headers = {"Content-Type": "application/json"}
    r_whole = requests.get(ensembl_server + whole_ext, headers=headers)
    
    if not r_whole.ok:
        r_whole.raise_for_status()
        sys.exit()

    whole_dna_sequence = r_whole.json().get("seq")
     
    if last_exon_info['end'] < last_exon_info['start']:
        reverse_complement_CRISPRi = reverse_complement(whole_dna_sequence)
        print(f"Reverse complement of Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = reverse_complement_CRISPRi
    else:
        print(f"Whole DNA sequence for the region surrounding the last exon ({last_exon_info['exon_id']}):")
        CRISPRtgSearch = whole_dna_sequence  

print(CRISPRtgSearch)

with open("sequence.fa", "w") as fa_file:
            fa_file.write(f">{gene_ids}_whole\n")
            fa_file.write(CRISPRtgSearch + "\n")

 # Progress update: 50% (after fetching the last exon info)
print("Progress: 50%",flush=True)

# Define the command and arguments
discover_command = [
    'java', '-Xmx4g', '-jar', 'FlashFry-assembly-1.15.jar', 
    'discover', 
    '--database', 'GRCh38_cas9ngg_database', 
    '--fasta', 'sequence.fa', 
    '--output', 'CRISPRtg.output'
]

# Define the second command and arguments
score_command = [
    'java', '-Xmx4g', '-jar', 'FlashFry-assembly-1.15.jar', 
    'score', 
    '--input', 'CRISPRtg.output', 
    '--output', 'CRISPRtg.output.scored.tsv', 
    '--scoringMetrics', 'doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,moreno2015', 
    '--database', 'GRCh38_cas9ngg_database'
]

# Run the first command
try:
    print("Running discover command...")
    discover_result = subprocess.run(discover_command, capture_output=True, text=True, check=True)
    print("Discover command executed successfully")
    print("Discover Command Output:\n", discover_result.stdout)
    print("Discover Command Errors:\n", discover_result.stderr)
except subprocess.CalledProcessError as e:
    print(f"Discover command failed with return code {e.returncode}")
    print(f"Discover Command Output:\n{e.output}")
    print(f"Discover Command Errors:\n{e.stderr}")
    exit(1)

  # Progress update: 60% (start running FlashFry 'discover' command)
    print("Progress: 60%")

# Run the second command
try:
    print("Running score command...")
    score_result = subprocess.run(score_command, capture_output=True, text=True, check=True)
    print("Score command executed successfully")
    print("Score Command Output:\n", score_result.stdout)
    print("Score Command Errors:\n", score_result.stderr)
except subprocess.CalledProcessError as e:
    print(f"Score command failed with return code {e.returncode}")
    print(f"Score Command Output:\n{e.output}")
    print(f"Score Command Errors:\n{e.stderr}")
    exit(1)

    # Progress update: 70% (after successfully running the discover command)
    print("Progress: 70%",flush=True)

# Define a function to calculate distance from exon for a target sequence
def calculate_distance(row):
    target_seq = row['target']
    orientation = row['orientation']
    print("target seq")
    print(target_seq)

    if orientation == "FWD":
        # Take 5 bp from the end of the target sequence
        segment = target_seq
    else:
        # Reverse complement the target sequence if orientation is RVS
        segment = str(Seq(target_seq).reverse_complement())
    
    if last_exon_info['end'] < last_exon_info['start']:
        rvs_seq = reverse_complement_exon[-15:]
        last_exon_position = reverse_complement_CRISPRi.rfind(rvs_seq)
        last_letter_position = last_exon_position + len(rvs_seq) - 1
        # Find the position of the segment in the whole DNA sequence
        target_position = reverse_complement_CRISPRi.find(segment)
    else:
        fwd_seq = exon_dna_sequence[-15:]
        last_exon_position = whole_dna_sequence.rfind(fwd_seq)
        last_letter_position = last_exon_position + len(fwd_seq) - 1
        # Find the position of the segment in the whole DNA sequence
        target_position = whole_dna_sequence.find(segment)


    #print("SegmentSearch")
    #print(segment)
    #print("This is where CRISPR is targeting")
    #print(target_position)
    #print("last_exon_position")
    #print(last_exon_position)
    #print("last letter position")
    #print(last_letter_position)
    # Calculate the distance from the last letter of the exon sequence
    if orientation == "FWD":
        distance = target_position - last_letter_position + 16
    else:
        distance = target_position - last_letter_position + 5
        #if distance < 0:
            #distance = distance - 1
    print(distance)
    return distance

if last_exon_info['end'] < last_exon_info['start']:
        rvs_seq = reverse_complement_exon[-15:]
        print(rvs_seq)
        last_exon_position = reverse_complement_CRISPRi.rfind(rvs_seq)
        last_letter_position = last_exon_position + len(rvs_seq) - 1
else:
        fwd_seq = exon_dna_sequence[-15:]
        print(fwd_seq)
        last_exon_position = whole_dna_sequence.rfind(fwd_seq)
        last_letter_position = last_exon_position + len(fwd_seq) - 1

   # Progress update: 80% (after calculating the exon positions)
print("Progress: 80%",flush=True)

# Print the positions for debugging
print("Last exon sequence position in whole DNA sequence:", last_exon_position)
print("Last letter position of exon sequence:", last_letter_position)
print("CRISPR Search Input Seq", whole_dna_sequence)

   # Progress update: 90% (after calculating the exon positions)
print("Progress: 90%",flush=True)
try:
    df = pd.read_csv('CRISPRtg.output.scored.tsv', sep='\t')
    print("DataFrame loaded successfully:")
    print(df.head())
    
    # Apply the distance calculation function
    df['Distance from Exon'] = df.apply(calculate_distance, axis=1)

    # Ensure 'Distance from Exon' column is present and doesn't contain None values
    if 'Distance from Exon' in df.columns:
        df['Distance from Exon'] = df['Distance from Exon'].fillna(0)  # Replace None with 0 or appropriate value
        df['Distance from Exon'] = df['Distance from Exon'].abs()  # Apply abs() after handling None values
    else:
        raise KeyError("'Distance from Exon' column not found in the DataFrame.")

    # Sort the DataFrame based on specified criteria
    filtered_df = df.sort_values(
        by=[
            'Hsu2013',
            'DoenchCFD_maxOT',
            'DoenchCFD_specificityscore',
            'Moreno-Mateos2015OnTarget',
            'Doench2014OnTarget',
            'Distance from Exon'
        ],
        ascending=[False, True, False, False, False, True]
    ).head(20)
    
    # Filter the DataFrame based on Distance from Exon
    filtered_df = filtered_df.iloc[filtered_df['Distance from Exon'].abs().argsort()[:20]]

    

    # Select the specified columns
    selected_columns = [
        'target', 
        'orientation', 
        'Doench2014OnTarget', 
        'DoenchCFD_maxOT', 
        'DoenchCFD_specificityscore', 
        'Hsu2013',
        'Moreno-Mateos2015OnTarget',
        'Distance from Exon'
        
    ]
    columns_to_multiply = [
    'Doench2014OnTarget', 
    'DoenchCFD_maxOT', 
    'DoenchCFD_specificityscore', 
    'Moreno-Mateos2015OnTarget'
    ]

    filtered_df[columns_to_multiply] = filtered_df[columns_to_multiply] * 100
    
    # Print the DataFrame with the new column
    print(filtered_df[selected_columns])

    columns_to_round = [
    'Doench2014OnTarget', 
    'DoenchCFD_maxOT', 
    'DoenchCFD_specificityscore', 
    'Hsu2013', 
    'Moreno-Mateos2015OnTarget'
    ]
    filtered_df[columns_to_round] = filtered_df[columns_to_round].round(1)

    filtered_df = filtered_df[selected_columns]
    print(filtered_df[selected_columns])

    # Save the filtered DataFrame to a CSV file
    filtered_csv_file = f'{gene_ids}_CRISPR_tgts.csv'
    filtered_df.to_csv(filtered_csv_file, index=False)
    print(f"Filtered DataFrame saved to {filtered_csv_file}")

    # Progress update: 100% (after processing and saving the filtered DataFrame)
    print("Progress: 100%", flush=True)

except FileNotFoundError:
    print("CRISPRtg.output.scored.tsv file not found.")
except pd.errors.EmptyDataError:
    print("CRISPRtg.output.scored.tsv file is empty.")
except KeyError as e:
    print(f"One or more specified columns are not found in the DataFrame: {e}")
except Exception as e:
    print(f"An error occurred while processing the DataFrame: {e}")


import json

def save_variables(gene_ids, last_exon_seq, gdna_sequence):
    variables = {
        "gene_ids": gene_ids,
        "last_exon_seq": last_exon_seq,
        "gdna_sequence": gdna_sequence
    }
    with open('variables.json', 'w') as file:
        json.dump(variables, file)

# Example usage
save_variables(gene_ids, last_exon_seq, gdna_sequence)

