from Events import *
from read_in_FASTA import read_in_FASTA


def KT_to_FASTA(KT: {str: Chromosome}, genome_path: str, chr_name_file: str, output_path: str):
    chr_name_conversion = {}
    full_name_list = []
    with open(chr_name_file) as fp_read:
        for line in fp_read:
            line = line.replace("\n", "").split("\t")
            chr_name_conversion[line[0]] = line[1]
            full_name_list.append(line[1])

    sequence_dict = read_in_FASTA(genome_path, full_name_list)

    output_dict = {}
    for chr_name, chr_obj in KT.items():
        new_sequence = []
        # telomere 1
        new_sequence.append('N' * chr_obj.t1_len)
        # p-arm
        for segment in chr_obj.p_arm.segments:
            if segment.direction():
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.start: segment.end + 1])
            else:
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.end: segment.start + 1]
                                    [::-1])
        # centromere
        for segment in chr_obj.centromere.segments:
            if segment.direction():
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.start: segment.end + 1])
            else:
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.end: segment.start + 1]
                                    [::-1])
        # q-arm
        for segment in chr_obj.q_arm.segments:
            if segment.direction():
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.start: segment.end + 1])
            else:
                new_sequence.append(sequence_dict[chr_name_conversion[segment.chr]][segment.end: segment.start + 1]
                                    [::-1])
        # telomere 2
        new_sequence.append('N' * chr_obj.t2_len)

        output_dict[chr_name] = ''.join(new_sequence)

    # output
    with open(output_path, 'w') as fp_write:
        for header, sequence in output_dict.items():
            fp_write.writelines(">{}\n".format(header))
            fp_write.writelines(sequence + "\n")
