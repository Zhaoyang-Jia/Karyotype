from Main.read_in_FASTA import read_in_FASTA


def Get_Telomere(genome_path, chr_name_file, telomere_output_path, length_output_path):
    """
    Detect and Return Telomere region, and total length of the chromosome
    e.g. NNNATCGNN -> 3 6 with length 9
    :param genome_path: FASTA genome
    :param chr_name_file: .txt file with (simpler name \t original name)
    :param telomere_output_path: .txt file with (simpler_chr_name \t coding_start_index \t coding_end_index)
    :param length_output_path: .txt file with (simpler_chr_name \t total_length)
    :return: None
    """
    # record the chromosomes of interest
    chr_of_interest = {}
    with open(chr_name_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            chr_of_interest[line[1]] = line[0]

    sequence_dict = read_in_FASTA(genome_path, chr_of_interest.keys())
    telomere_dict = {}

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

    # output
    with open(telomere_output_path, 'w') as fp_write:
        for telomere_header in telomere_dict:
            line = "{}\t{}\t{}\n".format(chr_of_interest[telomere_header],
                                         telomere_dict[telomere_header][0],
                                         telomere_dict[telomere_header][1])
            fp_write.writelines(line)
    with open(length_output_path, 'w') as fp_write:
        for sequence_header in sequence_dict:
            line = "{}\t{}\n".format(chr_of_interest[sequence_header],
                                     len(sequence_dict[sequence_header]))
            fp_write.writelines(line)
