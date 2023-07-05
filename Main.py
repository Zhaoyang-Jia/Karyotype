from Preparation.Get_All_Chr_Names import Get_All_Chr_Names
from Preparation.Extract_Chr import Extract_Chr

# Get_All_Chr_Names("GCF_000001405.26_GRCh38_genomic.fasta")
chr_to_extract = {'NC_000001.11 Homo sapiens chromosome 1, GRCh38 Primary Assembly': 'chr1',
                  'NC_000002.12 Homo sapiens chromosome 2, GRCh38 Primary Assembly': 'chr2'}
Extract_Chr("GCF_000001405.26_GRCh38_genomic.fasta",
            'REF_chr1_chr2.fasta',
            chr_to_extract)

