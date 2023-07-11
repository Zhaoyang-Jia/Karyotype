def Prepare_Raw_KT(chr_of_interest, genome_index_file, output_path):
    """
    Compose unedited KT
    :param chr_of_interest: list of chromosomes to generate KT
    :param output_path: .txt metadata file containing the unedited KT for each chromosome
    :param genome_index_file: .txt metadata file containing centromere, telomere, and genome length information
    :return:
    """
    segment_index = 1
    segments = []
    KT = {}
    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] in chr_of_interest:
                chromosome_origin = line[0]
                KT[chromosome_origin] = []
                # telomere A
                segment_name = chromosome_origin.replace('Chr', 'T') + 'A'
                start_index = 0
                end_index = int(line[2]) - 1
                segments.append("{}\t{}:{}-{}".format(segment_name, chromosome_origin, start_index, end_index))
                KT[chromosome_origin].append("+" + segment_name)

                # telomere B
                segment_name = chromosome_origin.replace('Chr', 'T') + 'B'
                start_index = int(line[5]) + 1
                end_index = int(line[1]) - 1
                segments.append("{}\t{}:{}-{}".format(segment_name, chromosome_origin, start_index, end_index))
                KT[chromosome_origin].append("+" + segment_name)

                # centromere
                segment_name = chromosome_origin.replace('Chr', 'C')
                start_index = int(line[3]) + 1
                end_index = int(line[4]) - 1
                segments.append("{}\t{}:{}-{}".format(segment_name, chromosome_origin, start_index, end_index))
                KT[chromosome_origin].insert(1, "+" + segment_name)

                # p-arm
                segment_name = str(segment_index)
                start_index = int(line[2])
                end_index = int(line[3])
                segments.append("{}\t{}:{}-{}".format(segment_name, chromosome_origin, start_index, end_index))
                KT[chromosome_origin].insert(1, "+" + segment_name)
                segment_index += 1

                # q-arm
                segment_name = str(segment_index)
                start_index = int(line[4])
                end_index = int(line[5])
                segments.append("{}\t{}:{}-{}".format(segment_name, chromosome_origin, start_index, end_index))
                KT[chromosome_origin].insert(3, "+" + segment_name)
                segment_index += 1

    # output to file
    with open(output_path, 'w') as fp_write:
        for segment in segments:
            fp_write.writelines(segment + "\n")

        fp_write.writelines("\n")

        for chromosome in KT:
            fp_write.writelines(",".join(KT[chromosome]) + "\n")

    # prepare for return

    return segments, KT
