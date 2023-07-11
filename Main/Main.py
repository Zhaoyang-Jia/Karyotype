import json
import argparse

from Events import *
from Output_KT import Output_KT
from KT_to_FASTA import KT_to_FASTA

parser = argparse.ArgumentParser()
parser.add_argument("JSON_file", help="path to the task-description JSON")
args = parser.parse_args()

with open(args.JSON_file) as file:
    instruction = json.load(file)

mode = instruction['mode']

# manual mode
if mode == 'manual':
    KT = Prepare_Raw_KT(instruction['chromosomes'], "../Metadata/Full_Genome_Indices.txt")
    events = instruction['events']
    for event in events:
        if event['type'] == 'deletion':
            deletion(KT, 0, "p", event['start'], event['end'])
        elif event['type'] == 'duplication':
            duplication(KT, 0, "p", event['start'], event['end'])
        elif event['type'] == 'inversion':
            inversion(KT, 0, "p", event['start'], event['end'])

    Output_KT(KT, "test_manual_mode_KT_out.txt")
    KT_to_FASTA(KT, genome_path="../Genomes/GCF_000001405.26_GRCh38_genomic.fasta",
                genome_index_file="../Metadata/Full_Genome_Indices.txt",
                chr_name_file="../Metadata/Chr_Names.txt",
                output_path="test_manual_mode_FASTA_out.fasta")
