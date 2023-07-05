def Extract_Chr(genome_path, new_file_path, header_to_extract):
    with open(genome_path) as fp_read:
        with open(new_file_path, 'w') as fp_write:
            segment_found = False
            for line in fp_read:
                if segment_found:
                    if line[0] != ">":
                        fp_write.writelines(line)
                    else:
                        segment_found = False
                if line[0] == ">":
                    hg38_header = line[1:].replace('\n', '')
                    if hg38_header in header_to_extract.keys():
                        segment_found = True
                        fp_write.writelines(">{}\n".format(header_to_extract[hg38_header]))
