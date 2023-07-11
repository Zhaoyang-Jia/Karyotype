import re


def Get_All_Chr_Names(genome_path, output_path):
    """
    Extract headers for the FASTA genome
    :param genome_path: FASTA genome
    :param output_path: .txt file with headers
    :return: None
    """
    with open(genome_path) as fp_read:
        with open(output_path, 'w') as fp_write:
            for line in fp_read:
                if line[0] == ">":
                    fp_write.writelines(line[1:])


def Extract_Whole_Chr_Names(header_file, output_path):
    """
    Only get the whole chromosome names, and notate with simpler names
    :param header_file: contains headers from the genome FASTA
    :param output_path: .txt file with (simpler name \t original name)
    :return: None
    """
    pattern = r'chromosome ([XY\d]+), GRCh38 Primary Assembly'
    with open(header_file) as fp_read:
        with open(output_path, 'w') as fp_write:
            for line in fp_read:
                match = re.findall(pattern, line)
                if match:
                    new_line = "Chr{}\t{}".format(match[0], line)
                    fp_write.writelines(new_line)
