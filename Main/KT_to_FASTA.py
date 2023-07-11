from Events import *


# ONLY WORKS ON Chr1 AT THE MOMENT
def KT_to_FASTA(KT, genome_path, genome_index_file, chr_name_file, output_path):
    chr_name_conversion = {}
    full_name_list = []
    with open(chr_name_file) as fp_read:
        for line in fp_read:
            line = line.replace("\n", "").split("\t")
            chr_name_conversion[line[0]] = line[1]
            full_name_list.append(line[1])

    # read in genome
    sequence_dict = {}
    with open(genome_path) as fp_read:
        recording = False
        for line in fp_read:
            if line[0] == '>':
                # header line
                if recording:
                    sequence_dict[header] = ''.join(segment_sequence)
                    recording = False
                if line[1:].replace('\n', '') in full_name_list:
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

    # construct new chromosome sequence
    new_sequence = []
    # telomere 1
    new_sequence.append(sequence_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly'][0:10000])
    # p-arm
    for segment in KT[0]['p']:
        start = segment.start
        end = segment.end
        new_sequence.append(sequence_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly']
                            [start: end + 1])
    # centromere
    new_sequence.append(sequence_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly']
                            [122026459: 125184588])
    # q-arm
    for segment in KT[0]['q']:
        start = segment.start
        end = segment.end
        new_sequence.append(sequence_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly']
                            [start: end + 1])
    # telomere 2
    new_sequence.append(sequence_dict['NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly']
                        [248946421:])

    with open(output_path, "w") as fp_write:
        fp_write.writelines(">1\n")
        fp_write.writelines("".join(new_sequence))
