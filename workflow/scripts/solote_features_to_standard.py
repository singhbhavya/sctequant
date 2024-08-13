import pandas as pd
import os
import sys
import gzip

def parse_gtf(gtf_file):
    gene_map = {}
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            columns = line.split('\t')
            if columns[2] == "gene":
                info = columns[8]
                gene_id = info.split('gene_id "')[1].split('";')[0]
                gene_name = info.split('gene_name "')[1].split('";')[0]
                gene_map[gene_name] = gene_id
    return gene_map

def process_features(features_file, gene_map, output_file):
    with open(features_file, 'r') as file:
        lines = file.readlines()

    output_data = []

    for line in lines:
        columns = line.strip().split('\t')
        gene_id = gene_map.get(columns[0], columns[0])  # Use the gene name if no match in gene_map
        output_data.append([gene_id, columns[1], "Gene Expression"])

    with open(output_file, 'w') as out_file:
        for row in output_data:
            out_file.write('\t'.join(row) + '\n')

def rename_file(original_file):
    directory = os.path.dirname(original_file)
    original_basename = os.path.basename(original_file)
    new_name = os.path.join(directory, 'original_' + original_basename)
    os.rename(original_file, new_name)
    return new_name

def process_files(gtf_file, features_files):
    # Parse GTF file
    print("Parsing GTF file...")
    gene_map = parse_gtf(gtf_file)
    print("GTF file parsed successfully.")

    for features_file in features_files:
        print(f"Processing file: {features_file}")
        # Rename original features file to avoid overwriting
        original_features_file = rename_file(features_file)
        
        # Define the path for the processed file
        processed_features_file = features_file  # The processed file will replace the original file
        
        # Process features file
        temp_output_file = processed_features_file.replace(".tsv", "_temp.tsv")
        process_features(original_features_file, gene_map, temp_output_file)
        
        # Gzip the output file
        with open(temp_output_file, 'rb') as f_in, gzip.open(processed_features_file + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        
        # Remove the temporary file
        os.remove(temp_output_file)

        print(f"File processed and saved as: {processed_features_file}.gz")

    print("Processing complete. All output files have been saved in the same directories as the original files.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py path_to_gtf_file.gtf path_to_features_file_or_list.txt")
        sys.exit(1)

    gtf_file = sys.argv[1]
    input_file_or_list = sys.argv[2]

    # Check if the input is a single features file or a list of files
    if os.path.isfile(input_file_or_list) and input_file_or_list.endswith(".tsv"):
        features_files = [input_file_or_list]
    else:
        with open(input_file_or_list, 'r') as file:
            features_files = [line.strip() for line in file if line.strip()]

    process_files(gtf_file, features_files)

