import json
import argparse
import random

from Events import *
from Output_KT import Output_KT
from KT_to_FASTA import KT_to_FASTA

parser = argparse.ArgumentParser()
parser.add_argument("JSON_file", help="path to the task-description JSON")
args = parser.parse_args()

with open(args.JSON_file) as file:
    instruction = json.load(file)

KT = Prepare_Raw_KT(instruction['chromosomes'], "../Metadata/Full_Genome_Indices.txt")
mode = instruction['mode']

# manual mode
if mode == 'manual':
    events = instruction['events']
    for event in events:
        event_type = event['type']
        if event['arm'] == 'p':
            event_arm = KT[event['chromosome']].p_arm
        elif event['arm'] == 'q':
            event_arm = KT[event['chromosome']].q_arm
        else:
            raise ValueError('arm selection must be either p or q')

        event_start_index = event['start']  # indices are 0-indexed based on the arm
        event_end_index = event['end']
        if event_start_index < 0:
            raise ValueError('start index must be greater or equal to 0')
        if event_end_index > len(event_arm):
            raise ValueError('end index is greater than the end of the selected arm')

        if event['type'] == 'deletion':
            deletion(event_arm, event['start'], event['end'])
        elif event['type'] == 'duplication':
            duplication(event_arm, event['start'], event['end'])
        elif event['type'] == 'inversion':
            inversion(event_arm, event['start'], event['end'])

    Output_KT(KT, "KT_" + instruction['job_name'] + ".txt")
    KT_to_FASTA(KT, genome_path="../Genomes/GCF_000001405.26_GRCh38_genomic.fasta",
                genome_index_file="../Metadata/Full_Genome_Indices.txt",
                chr_name_file="../Metadata/Chr_Names.txt",
                output_path=instruction['job_name'] + ".fasta")

# automatic mode
# elif mode == 'automatic':
#     event_settings = instruction['event_setting']
#     event_weights = [event_settings[0]['ratio'], event_settings[1]['ratio'], event_settings[2]['ratio']]
#     if sum(event_weights) != 1:
#         print('ratio error, must sum up to 1')
#
#     for event_index in range(instruction['number_of_events']):
#         # choose event
#         current_event = random.choices([0, 1, 2], weights=event_weights)[0]
#         # choose length
#         current_event_length = random.randint(event_settings[current_event]['min_size'],
#                                               event_settings[current_event]['max_size'])
#         # choose arm
#         current_arm = random.choices(["p", "q"])[0]
#
#         # choose start location
#         current_event_start_location = random.randint(10000, 10000000 - current_event_length)
#
#         # perform event
#         if current_event == 0:
#             deletion(KT, 0, current_arm,
#                      current_event_start_location, current_event_start_location + current_event_length)
#         elif current_event == 1:
#             duplication(KT, 0, current_arm,
#                         current_event_start_location, current_event_start_location + current_event_length)
#         elif current_event == 2:
#             inversion(KT, 0, current_arm,
#                       current_event_start_location, current_event_start_location + current_event_length)
#
#         Output_KT(KT, "KT_" + instruction['job_name'] + ".txt")
#         KT_to_FASTA(KT, genome_path="../Genomes/GCF_000001405.26_GRCh38_genomic.fasta",
#                     genome_index_file="../Metadata/Full_Genome_Indices.txt",
#                     chr_name_file="../Metadata/Chr_Names.txt",
#                     output_path=instruction['job_name'] + ".fasta")
