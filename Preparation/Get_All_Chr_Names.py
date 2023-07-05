def Get_All_Chr_Names(genome_path):
    with open(genome_path) as fp_read:
        with open('Chr_Names.txt', 'w') as fp_write:
            for line in fp_read:
                if line[0] == ">":
                    fp_write.writelines(line[1:])
