"""
This script (Assignment 3) analyzes upstream promoter sequences for co-regulated maize genes.
Workflow:
1. Parse genome FASTA files and GFF3 file to extract upstream promoter sequences.
2. Identify transcription factor binding motifs in target gene promoters.
3. Perform the same analysis on five random sets of genes for comparison.
4. Write motif counts for target and random sets to an output file.
"""

import os
import random
import logging
from Bio import SeqIO
import argparse
import re

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Function to load genome data from FASTA files
def load_genome(fasta_dir):
    """
    Loads genome data from a directory containing FASTA files and normalizes chromosome keys.

    Args:
        fasta_dir (str): Path to the directory containing genome FASTA files.

    Returns:
        dict: A dictionary where keys are normalized chromosome identifiers (e.g., "1", "2", etc.)
              and values are Bio.SeqRecord objects containing sequence data.
    """
    genome = {}  # Initialize an empty dictionary to store genome data

    # Iterate through all files in the specified directory
    for file in os.listdir(fasta_dir):
        # Only process files with ".fa" or ".fasta" extensions
        if file.endswith(".fa") or file.endswith(".fasta"):
            filepath = os.path.join(fasta_dir, file)  # Get the full path of the file

            # Parse the FASTA file using Biopython's SeqIO parser
            for record in SeqIO.parse(filepath, "fasta"):
                # Normalize chromosome identifiers by removing prefixes and leading zeros
                # Example: "chromosome_01" -> "1", "chr2" -> "2"
                normalized_key = record.id.replace("chromosome_", "").replace("chr", "").lstrip("0").upper()

                # Add the normalized key and corresponding sequence record to the genome dictionary
                genome[normalized_key] = record

    # Log the number of chromosomes loaded for debugging or confirmation
    logging.info(f"Loaded genome data from {len(genome)} chromosomes.")

    return genome  # Return the populated genome dictionary


# Function to parse GFF3 file for gene data
def parse_gff3(file_path):
    """
    Parses a GFF3 file to extract and normalize gene data.

    Args:
        file_path (str): Path to the GFF3 file.

    Returns:
        dict: A dictionary where keys are normalized gene IDs, and values are dictionaries
              containing gene information (chromosome, start, end, and strand).
    """
    genes = {}  # Initialize an empty dictionary to store gene data

    # Open the GFF3 file for reading
    with open(file_path, 'r') as file:
        for line in file:
            # Skip header or comment lines
            if line.startswith('#'):
                continue
            
            # Split the line into tab-separated columns
            cols = line.strip().split('\t')

            # Ensure the line contains at least 9 columns; otherwise, skip it
            if len(cols) < 9:
                continue

            # Process only lines that represent genes (e.g., `cols[2] == 'gene'`)
            if cols[2] == 'gene':
                # Parse the attributes column (column 9) into a dictionary
                attributes = {
                    item.split('=')[0]: item.split('=')[1]
                    for item in cols[8].split(';') if '=' in item
                }

                # Extract and normalize the gene ID
                if 'ID' in attributes:
                    # Remove prefixes like "gene:" if present
                    gene_id = attributes['ID'].split(':')[1] if ':' in attributes['ID'] else attributes['ID']
                    # Remove version numbers (e.g., ".1") and convert to uppercase for consistency
                    gene_id = gene_id.split('.')[0].upper()

                    # Normalize chromosome IDs by removing prefixes and leading zeros
                    seqid = cols[0].replace("chromosome_", "").replace("chr", "").lstrip("0").upper()

                    # Only include genes located on chromosomes 1–10
                    if seqid in [str(i) for i in range(1, 11)]:
                        # Store gene data in the dictionary
                        genes[gene_id] = {
                            'seqid': seqid,          # Chromosome identifier
                            'start': int(cols[3]),  # Start position
                            'end': int(cols[4]),    # End position
                            'strand': cols[6]       # Strand information ('+' or '-')
                        }

    # Log the total number of genes parsed for debugging or confirmation
    logging.info(f"Parsed {len(genes)} genes from GFF3 file.")

    return genes  # Return the dictionary containing parsed gene data

# Function to extract upstream promoter sequences
def extract_promoter(genome, gene_info, length=500):
    """
    Extracts the upstream promoter sequence for a given gene.

    Args:
        genome (dict): A dictionary containing chromosome sequences.
        gene_info (dict): Information about the gene (chromosome, start, end, strand).
        length (int, optional): Length of the upstream promoter region to extract. Defaults to 500.

    Returns:
        str: The extracted promoter sequence or None if not valid.
    """
    chrom_key = gene_info.get('seqid')
    if chrom_key not in genome:
        logging.warning(f"Chromosome {chrom_key} not found in genome. Skipping...")
        return None

    # Determine start and end positions
    if gene_info['strand'] == '+':
        start = max(0, gene_info['start'] - length)  # Ensure start is not negative
        end = gene_info['start'] - 1
    else:
        start = gene_info['end'] + 1
        end = min(gene_info['end'] + length, len(genome[chrom_key].seq))  # Ensure end is within sequence bounds

    # Extract promoter sequence
    promoter = genome[chrom_key].seq[start:end].upper()

    # Handle ambiguous bases (truncate at first occurrence of 'N')
    if 'N' in promoter:
        first_n_index = promoter.find('N')
        logging.warning(f"Ambiguous bases found. Truncating promoter for gene {gene_info.get('seqid')}.")
        promoter = promoter[:first_n_index]

    # Log completely ambiguous promoters
    if not promoter:  # Promoter is empty after truncation
        logging.warning(f"Completely ambiguous promoter for gene {gene_info.get('seqid')}. Skipping...")
        return None

    # Return reverse complement for negative strand
    return promoter.reverse_complement() if gene_info['strand'] == '-' else promoter


# Function to compile regex patterns for motifs
def compile_motif_patterns(motifs):
    """
    Compiles regex patterns for a list of motifs, converting ambiguous bases (e.g., [AC]) to regex groups.

    Args:
        motifs (list): List of motifs in string format (e.g., "A[CT]G").

    Returns:
        dict: Dictionary with motifs as keys and their compiled regex patterns as values.
    """
    patterns = {}  # Dictionary to store compiled regex patterns

    # Iterate through each motif and compile it into a regex pattern
    for motif in motifs:
        normalized_motif = motif.upper()  # Convert motif to uppercase for consistency
        # Convert ambiguous bases (e.g., [AC]) to regex groups (e.g., (A|C))
        regex = re.sub(r'\[([A-Z]+)\]', lambda x: f"({'|'.join(x.group(1))})", normalized_motif)
        patterns[normalized_motif] = re.compile(regex)  # Compile the regex pattern
        logging.debug(f"Compiled regex for motif: {normalized_motif} -> {regex}")  # Log the compiled regex

    return patterns  # Return the dictionary of compiled patterns

# Function to read gene IDs
def read_genes(file_path):
    """
    Reads and normalizes gene IDs from the input file.

    Args:
        file_path (str): Path to the file containing gene IDs.

    Returns:
        list: A list of normalized gene IDs, converted to uppercase and stripped of version numbers.
    """
    # Check if the file exists; log an error and return an empty list if not
    if not os.path.exists(file_path):
        logging.error(f"Gene file '{file_path}' not found.")
        return []

    # Read and normalize gene IDs
    with open(file_path, 'r') as file:
        genes = [
            line.strip().split('.')[0].upper()  # Remove version numbers (e.g., ".1") and convert to uppercase
            for line in file.readlines()
        ]

    # Log the total number of genes loaded for debugging
    logging.debug(f"Loaded {len(genes)} target genes.")

    return genes  # Return the list of normalized gene IDs


# Function to read promoter motifs
def read_promoters(file_path):
    """
    Reads promoter motifs from the input file.

    Args:
        file_path (str): Path to the file containing promoter motifs.

    Returns:
        list: A list of motifs, each stripped of leading/trailing whitespace.
    """
    # Check if the file exists; log an error and return an empty list if not
    if not os.path.exists(file_path):
        logging.error(f"Promoter file '{file_path}' not found.")
        return []

    # Read and process each line from the promoter motifs file
    with open(file_path, 'r') as file:
        motifs = [line.strip() for line in file.readlines()]  # Strip whitespace from each line

    # Log the total number of motifs loaded for debugging
    logging.debug(f"Loaded {len(motifs)} promoter motifs.")

    return motifs  # Return the list of motifs


# Main function
def main():
    parser = argparse.ArgumentParser(description="Analyze maize promoter sequences.")
    parser.add_argument('fasta_dir', help="Directory containing the maize genome FASTA files.")
    parser.add_argument('gff3', help="Path to the maize genome GFF3 file.")
    parser.add_argument('genes', help="Path to the file with target gene IDs.")
    parser.add_argument('promoters', help="Path to the file with known motifs.")
    parser.add_argument('--promoter_length', type=int, default=500,
                        help="Length of upstream promoter region to extract (default: 500 nt).")
    parser.add_argument('--random_sets', type=int, default=5,
                        help="Number of random gene sets to analyze (default: 5).")
    args = parser.parse_args()

    try:
        genome = load_genome(args.fasta_dir)
        gene_info = parse_gff3(args.gff3)
        target_genes = read_genes(args.genes)
        motifs = read_promoters(args.promoters)
        motif_patterns = compile_motif_patterns(motifs)
    except Exception as e:
        logging.error(e)
        return

    # Validate matching IDs
    unmatched_genes = [gene for gene in target_genes if gene not in gene_info]
    if unmatched_genes:
        logging.warning(f"{len(unmatched_genes)} out of {len(target_genes)} target genes did not match the GFF3 data.")
        logging.warning(f"Examples of unmatched genes: {unmatched_genes[:10]}")
    else:
        logging.info("All target genes matched successfully.")

    target_counts = {}

    for gene in target_genes:
        if gene in gene_info:
            # Extract the promoter sequence
            promoter = extract_promoter(genome, gene_info[gene], args.promoter_length)
            
            if promoter is None:
                logging.warning(f"Promoter not found for gene {gene}.")
                continue
            
            # Write promoter sequence to the debug file
            with open("promoters_debug.txt", "a") as debug_file:
                debug_file.write(f">{gene}\n{promoter}\n")
            
            # Count motifs in the promoter sequence
            counts = {motif: len(pattern.findall(str(promoter))) for motif, pattern in motif_patterns.items()}
            
            # Log if no motifs are found
            if all(count == 0 for count in counts.values()):
                logging.warning(f"No motifs found in the promoter of gene {gene}.")
            
            # Accumulate motif counts across all target genes
            for motif, count in counts.items():
                target_counts[motif] = target_counts.get(motif, 0) + count

            

    random_counts = []
    all_genes = list(gene_info.keys())
    for _ in range(args.random_sets):
        sample_genes = random.sample(all_genes, len(target_genes))
        sample_counts = {}
        for gene in sample_genes:
            if gene in gene_info:
                promoter = extract_promoter(genome, gene_info[gene], args.promoter_length)
                if promoter:
                    counts = {motif: len(pattern.findall(str(promoter))) for motif, pattern in motif_patterns.items()}
                    for motif, count in counts.items():
                        sample_counts[motif] = sample_counts.get(motif, 0) + count
        random_counts.append(sample_counts)

    # Write results to file
    with open("motif_counts.txt", 'w') as output_file:
        # Define the header row
        headers = ["Motif", "Selected_Genes"] + [f"Random_Set_{i+1}" for i in range(args.random_sets)]
        col_widths = {col: max(len(col), 15) for col in headers}  # Define minimum column widths
        for motif in motifs:
            col_widths["Motif"] = max(col_widths["Motif"], len(motif) + 2)

        # Write the header row
        header_row = "".join(f"{col:<{col_widths[col]}}" for col in headers)
        output_file.write(header_row + "\n\n")  # Add a space after the header row

        # Write the rows of data
        for motif in motifs:
            row = [motif, target_counts.get(motif, 0)]
            row += [counts.get(motif, 0) for counts in random_counts]
            row_formatted = "".join(f"{str(value):<{col_widths[header]}}" for value, header in zip(row, headers))
            output_file.write(row_formatted + "\n\n")  # Add an extra newline between rows


    logging.info("Analysis complete :) Results written to 'motif_counts.txt'.")

if __name__ == "__main__":
    main()

# ******************************************************************************************************************************************************
# I used the following command in the terminal to run the script:
# python Assignment3.py "C:\Users\BAROK STORE\OneDrive\Desktop\Assignment3" Zea_mays.B73_RefGen_v4.48.chr.gff3 zea_mays_genes.txt promoters.txt
#*******************************************************************************************************************************************************
# This assignment analyzes upstream promoter sequences for co-regulated maize genes to identify transcription factor binding motifs. It processes genome
# data from FASTA files and extracts gene information from a GFF3 file, normalizing gene IDs and chromosome identifiers for consistency. 
# The script reads a list of target genes and motifs, extracts upstream promoter regions (default 500 bp) for each target gene,
# and counts occurrences of motifs in the promoters. Additionally, it performs the same analysis on random sets of genes for comparison.
# Results, including motif counts for target and random gene sets, are written to a output file (`motif_counts`). 
# The script uses extensive logging for error handling, debugging, and ensuring transparency throughout the workflow.
#*******************************************************************************************************************************************************


