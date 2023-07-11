def Prepare_Raw_KT(chr_of_interest, genome_index_file):
    """
    Compose unedited KT
    :param chr_of_interest: list of chromosomes to generate KT
    :param genome_index_file: .txt metadata file containing centromere, telomere, and genome length information
    :return:
    """
    chromosomes = []
    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] in chr_of_interest:
                new_chromosome = {"p": [], "q": []}
                chromosome_origin = line[0]

                # p-arm
                start_index = int(line[2])
                end_index = int(line[3])
                new_chromosome['p'].append([chromosome_origin, start_index, end_index])

                # q-arm
                start_index = int(line[4])
                end_index = int(line[5])
                new_chromosome['q'].append([chromosome_origin, start_index, end_index])
                chromosomes.append(new_chromosome)

    return chromosomes
