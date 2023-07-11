from Preparation.Get_Chr_Names import Get_All_Chr_Names
from Preparation.Get_Chr_Names import Extract_Whole_Chr_Names
from Preparation.Get_Telomere import Get_Telomere
from Preparation.Genome_Indexing import Genome_Indexing


"""
Steps required to prepare the metadata files
REQUIRE manual input of centromere location
"""

Get_All_Chr_Names("../Genomes/GCF_000001405.26_GRCh38_genomic.fasta", "../Metadata/Header_Names.txt")

Extract_Whole_Chr_Names("../Metadata/Header_Names.txt", "../Metadata/Chr_Names.txt")

Get_Telomere("../Genomes/GCF_000001405.26_GRCh38_genomic.fasta",
             "../Metadata/Chr_Names.txt",
             "../Metadata/Telomere_Indices.txt",
             "../Metadata/Chr_Length.txt")
Genome_Indexing("../Metadata/Telomere_Indices.txt",
                "../Metadata/Centromere_Indices.txt",
                "../Metadata/Chr_Length.txt",
                "../Metadata/Full_Genome_Indices.txt")

# chr_to_extract = {'NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly': 'chr1',
#                   'NC_000002.12 Homo sapiens chromosome 2, GRCh38 Primary Assembly': 'chr2'}
# Extract_Chr("GCF_000001405.26_GRCh38_genomic.fasta",
#             'REF_chr1_chr2.fasta',
#             chr_to_extract)

# chr_to_extract = {'NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly': 'chr1'}
# Extract_Chr("GCF_000001405.26_GRCh38_genomic.fasta",
#             'REF_chr1.fasta',
#             chr_to_extract)

# with open('../REF_chr1.fasta') as fp_read:
#     sequence = []
#     for line in fp_read:
#         if line[0] == ">":
#             continue
#         else:
#             sequence.append(line.replace('\n', ''))
#
#     sequence_str = ''.join(sequence)
