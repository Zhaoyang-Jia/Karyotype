def Get_Telomere(genome_path):
    sequence_dict = {}
    telomere_dict = {}
    with open(genome_path) as fp_read:
        first_segment_met = False
        for line in fp_read:
            if line[0] == '>':
                if first_segment_met:
                    sequence_dict[header] = ''.join(segment_sequence)
                else:
                    first_segment_met = True

                segment_sequence = []
                header = line[1:].replace('\n', '')
                continue
            else:
                # read in the entire segment
                segment_sequence.append(line.replace('\n', ''))
    # record the last segment
    sequence_dict[header] = ''.join(segment_sequence)
    print('step1_done')

    # locate the telomere regions
    for itr_header in sequence_dict:
        sequence = sequence_dict[itr_header]
        start_index = 0
        end_index = len(sequence) - 1

        while start_index < len(sequence) and sequence[start_index] == 'N':
            start_index += 1
        while end_index >= 0 and sequence[end_index] == 'N':
            end_index -= 1

        telomere_dict[itr_header] = (start_index, end_index)

    return telomere_dict
