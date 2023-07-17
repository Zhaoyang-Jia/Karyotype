def read_in_FASTA(genome_path: str, chr_of_interest):
    sequence_dict = {}
    with open(genome_path) as fp_read:
        recording = False
        for line in fp_read:
            if line[0] == '>':
                # header line
                if recording:
                    sequence_dict[header] = ''.join(segment_sequence)
                    recording = False
                if line[1:].replace('\n', '') in chr_of_interest:
                    recording = True
                    segment_sequence = []
                    header = line[1:].replace('\n', '')
            else:
                # sequence line
                if recording:
                    segment_sequence.append(line.replace('\n', ''))
    # record the last segment
    if recording:
        sequence_dict[header] = ''.join(segment_sequence)

    return sequence_dict