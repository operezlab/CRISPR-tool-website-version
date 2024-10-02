import pandas as pd
import sys
import json

def process_selected_row(csv_file):
    try:
        # Load the selected row from the CSV file
        df = pd.read_csv(csv_file)
        print(f"Loaded selected row from CSV file: {csv_file}")
        
        # Display the DataFrame (you can replace this with your own processing logic)
        print("Selected Row Data:")
        print(df)
        
        # Ensure the expected columns are present
        expected_columns = ['target', 'orientation', 'Distance from Exon']
        missing_columns = [col for col in expected_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing columns in the CSV file: {', '.join(missing_columns)}")

        # Access and print the selected columns
        selected_target = df['target'].values[0] if 'target' in df.columns else 'N/A'
        selected_orientation = df['orientation'].values[0] if 'orientation' in df.columns else 'N/A'
        selected_distance_from_exon = df['Distance from Exon'].values[0] if 'Distance from Exon' in df.columns else 'N/A'
        
        print(f"Target: {selected_target}")
        print(f"Orientation: {selected_orientation}")
        print(f"Distance from Exon: {selected_distance_from_exon}")

        # Load additional variables from JSON
        def load_variables():
            with open('variables.json', 'r') as file:
                variables = json.load(file)
            return variables["gene_ids"], variables["last_exon_seq"], variables["gdna_sequence"]

        gene_ids, last_exon_seq, gdna_sequence = load_variables()

        def find_and_extract_sequence(gdna, target, upstream, downstream):
            target_index = gdna.find(target)
            if target_index == -1:
                return None, None
            upstream_sequence = gdna[max(0, target_index + len(target) - upstream):target_index + len(target)]
            downstream_sequence = gdna[target_index + len(target) - 1:target_index + len(target) - 1 + downstream]
            return upstream_sequence, downstream_sequence

        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return ''.join(complement[base] for base in reversed(seq))

        gene_name = gene_ids
        sgRNA_sequence = selected_target

        if selected_orientation == "RVS":
            is_sgRNA_inverted = "Y"
        else:
            is_sgRNA_inverted = "N"

        original_sgRNA = sgRNA_sequence

        if is_sgRNA_inverted == "Y":
            sgRNA_sequence = reverse_complement(original_sgRNA)

        if len(last_exon_seq) >= 20:
            protein_coding_exon = last_exon_seq[-20:]
        else:
            protein_coding_exon = last_exon_seq

        left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, 417, 0)

        if is_sgRNA_inverted == "Y":
            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:7], 0, 321)
        else:
            _, right_arm = find_and_extract_sequence(gdna_sequence, sgRNA_sequence[:18], 0, 321)

        left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
        if is_sgRNA_inverted == "Y":
            left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", left_arm)
        else:
            left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", left_arm)
        
        right_ai1_template = "CCTGCGGTGTCTTTGCTTrycatgtGGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCTxGTTTTAGAGCTAGAAATAGCAA"
        if is_sgRNA_inverted == "Y":
            right_ai1_seq = (
                right_ai1_template.replace("r", right_arm)
                .replace("y", original_sgRNA)
                .replace("x", original_sgRNA[:20])
            )
        else:
            right_ai1_seq = (
                right_ai1_template.replace("r", right_arm)
                .replace("y", original_sgRNA)
                .replace("x", original_sgRNA[:20])
            )
        
        # Save output to a .txt file
        processed_output_file = csv_file.replace('.csv', '_processed_output.txt')
        with open(processed_output_file, 'w') as file:
            file.write(f"Target: {selected_target}\n")
            file.write(f"Orientation: {selected_orientation}\n")
            file.write(f"Distance from Exon: {selected_distance_from_exon}\n")
            file.write(f">>{gene_name}-Left-AI1\n{left_ai1_seq}\n\n")
            file.write(f">>{gene_name}-Right-AI1\n{right_ai1_seq}\n\n")
            file.write(f"Left_arm sequence:\n{left_arm}\n")
            file.write(f"Right_arm sequence:\n{right_arm}\n")
        
        print(f"Processed output saved to {processed_output_file}")

    except Exception as e:
        print(f"An error occurred while processing the CSV file: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
        process_selected_row(csv_file)
    else:
        print("No CSV file provided.")
