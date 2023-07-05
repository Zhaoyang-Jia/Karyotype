import re


def Get_Telomere(genome_path):
    sequence_dict = {}
    telomere_dict = {}
    with open(genome_path) as fp_read:
        first_segment_met = False
        for line in fp_read:
            if line[0] == '>':
                if first_segment_met:
                    sequence_dict[header] = segment_sequence
                else:
                    first_segment_met = True

                segment_sequence = ""
                header = line[1:].replace('\n', '')
                continue
            else:
                # read in the entire segment
                segment_sequence += line.replace('\n', '')
    # record the last segment
    sequence_dict[header] = segment_sequence

    # locate the telomere regions
    for header in sequence_dict:
        sequence = sequence_dict[header]
        start_match = re.search(r'^N+', sequence)
        end_match = re.search(r'N+$', sequence)
        telomere_dict[header] = (start_match.start(), end_match.start())

    return telomere_dict
