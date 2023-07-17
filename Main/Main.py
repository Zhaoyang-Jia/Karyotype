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

chr_of_events = []
if instruction['chromosomes'][0] == "ALL":
    chr_of_events = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8', 'Chr9', 'Chr10', 'Chr11', 'Chr12',
                     'Chr13', 'Chr14', 'Chr15', 'Chr16', 'Chr17', 'Chr18', 'Chr19', 'Chr20', 'Chr21', 'Chr22', 'ChrX',
                     'ChrY']
else:
    chr_of_events = instruction['chromosomes']
KT = Prepare_Raw_KT(chr_of_events, "../Metadata/Full_Genome_Indices.txt")

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
                chr_name_file="../Metadata/Chr_Names.txt",
                output_path=instruction['job_name'] + ".fasta")
    # KT_to_FASTA(KT, genome_path="../Genomes/test_genome.fasta",
    #             chr_name_file="../Metadata/Chr_Names.txt",
    #             output_path=instruction['job_name'] + ".fasta")

# automatic mode
elif mode == 'automatic':
    event_settings = instruction['event_setting']

    # scale event weights
    sum_weight = 0.0
    for index in range(6):
        sum_weight += event_settings[index]['ratio']
    event_weights = []
    for index in range(6):
        event_weights.append(float(event_settings[index]['ratio']) / sum_weight)

    # perform events
    for event_index in range(instruction['number_of_events']):
        # choose event
        current_event = random.choices(range(6), weights=event_weights)[0]

        # choose chr
        current_chr1 = -1
        current_chr2 = -1
        chr_weights = []
        sum_length = 0.0
        for chromosome in chr_of_events:
            sum_length += len(chromosome)
        for chromosome in chr_of_events:
            chr_weights.append(float(len(chromosome)) / sum_length)
        current_chr1 = random.choices(chr_of_events, chr_weights)[0]
        if current_event == 5:
            # event is inter-chromosomal
            current_chr1_index = chr_of_events.index(current_chr1)
            possible_chr = chr_of_events
            possible_chr.remove(current_chr1)
            chr_weights.pop(current_chr1_index)
            current_chr2 = random.choices(possible_chr, chr_weights)[0]
        elif current_event == 4:
            current_chr2 = current_chr1

        # choose arm
        arm_weights = [float(KT[current_chr1].p_arm_len()) / len(KT[current_chr1]),
                       float(KT[current_chr1].q_arm_len()) / len(KT[current_chr1])]
        current_arm1_value = random.choices(["p", "q"], arm_weights)[0]
        current_arm1 = None
        current_arm2 = None
        if current_arm1_value == 'p':
            current_arm1 = KT[current_chr1].p_arm
        else:
            current_arm1 = KT[current_chr1].q_arm
        if current_event == 5:
            sum_arm_length = len(KT[current_chr2])
            arm_weights = [float(KT[current_chr2].p_arm_len()) / len(KT[current_chr2]),
                           float(KT[current_chr2].q_arm_len()) / len(KT[current_chr2])]
            current_arm2_value = random.choices(["p", "q"], arm_weights)[0]
            if current_arm2_value == 'p':
                current_arm2 = KT[current_chr2].p_arm
            else:
                current_arm2 = KT[current_chr2].q_arm
        elif current_event == 4:
            arm_weights = [float(KT[current_chr1].p_arm_len()) / len(KT[current_chr1]),
                           float(KT[current_chr1].q_arm_len()) / len(KT[current_chr1])]
            current_arm2_value = random.choices(["p", "q"], arm_weights)[0]
            if current_arm2_value == 'p':
                current_arm2 = KT[current_chr1].p_arm
            else:
                current_arm2 = KT[current_chr1].q_arm

        # choose length
        current_event1_length = -1
        current_event2_length = -1
        if current_event in range(4):
            current_event1_length = random.randint(event_settings[current_event]['min_size'],
                                                   event_settings[current_event]['max_size'])
            current_event1_length = min(current_event1_length, len(current_arm1) - 1)
        elif current_event == 4 or current_event == 5:
            current_event1_length = random.randint(event_settings[current_event]['min_size1'],
                                                   event_settings[current_event]['max_size1'])
            current_event1_length = min(current_event1_length, len(current_arm1) - 1)
            current_event2_length = random.randint(event_settings[current_event]['min_size2'],
                                                   event_settings[current_event]['max_size2'])
            current_event2_length = min(current_event2_length, len(current_arm2) - 1)

        # choose start location
        current_event_start_location1 = random.randint(0, len(current_arm1) - current_event1_length - 1)
        current_event_start_location2 = -1
        if current_event2_length != -1:
            current_event_start_location2 = random.randint(0, len(current_arm2) - current_event2_length - 1)

        # perform event
        if current_event == 0:
            deletion(current_arm1, current_event_start_location1,
                     current_event_start_location1 + current_event1_length)
        elif current_event == 1:
            duplication(current_arm1, current_event_start_location1,
                        current_event_start_location1 + current_event1_length)
        elif current_event == 2:
            inversion(current_arm1, current_event_start_location1,
                      current_event_start_location1 + current_event1_length)
        elif current_event == 3:
            duplication_inversion(current_arm1, current_event_start_location1,
                                  current_event_start_location1 + current_event1_length)
        elif current_event == 4 or current_event == 5:
            translocation_reciprocal(current_arm1, current_event_start_location1,
                                     current_event_start_location1 + current_event1_length,
                                     current_arm2, current_event_start_location2,
                                     current_event_start_location2 + current_event2_length)

    # output
    Output_KT(KT, "KT_" + instruction['job_name'] + ".txt")
    KT_to_FASTA(KT, genome_path="../Genomes/GCF_000001405.26_GRCh38_genomic.fasta",
                chr_name_file="../Metadata/Chr_Names.txt",
                output_path=instruction['job_name'] + ".fasta")
