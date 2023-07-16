def Genome_Indexing(telomere_file, centromere_file, length_file, output_path):
    """
    Compose all chromosomes' indices
    :param telomere_file: .txt metadata file containing start and end of telomere indices
    :param centromere_file: .txt metadata file containing start and end of centromere indices
    :param length_file: .txt metadata file containing chromosome's length
    :param output_path: .txt metadata file with (chr# \t length \t sequence_start_after_t1 \t sequence_end_before_C
    \t sequence_start_after_C \t sequence_start_before_t2)
    :return: None
    """
    indices_dict = {}
    with open(length_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            indices_dict[line[0]] = [int(line[1])]
    with open(telomere_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            indices_dict[line[0]].append(int(line[1]))
            indices_dict[line[0]].append(int(line[2]))
    with open(centromere_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            # as NCBI labels the start/end of the centromere index, but we want to label the last/first coding region
            # before/after the centromere, consistent with the telomere labeling
            indices_dict[line[0]].insert(2, int(line[1].replace(',', '')) - 1)
            indices_dict[line[0]].insert(3, int(line[2].replace(',', '')) + 1)

    with open(output_path, 'w') as fp_write:
        for chromosome in indices_dict:
            line = "{chr}\t{len}\t{telo_end}\t{centro_start}\t{centro_end}\t{telo_start}\n"\
                .format(chr=chromosome,
                        len=indices_dict[chromosome][0],
                        telo_end=indices_dict[chromosome][1],
                        centro_start=indices_dict[chromosome][2],
                        centro_end=indices_dict[chromosome][3],
                        telo_start=indices_dict[chromosome][4])
            fp_write.writelines(line)
