# Karyotype
Require downloading your own genome, get RefSeq Version FASTA at
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

If a custom genome is used, indexing can be done using Preparation/Prepare_Metadata.py.
However, the parameters need to be changed within the script.

Put the genome in /Genomes, and Main/Main.py should be ready to go

The "test_manual.json" should work directly with Main.py
just call `Main.py test_manual.json`

Current Support:
Intra-chromosomal deletion, inversion, duplication, duplication_inversion, and translocation are implemented and tested.
JSON works for both manual and automatic modes. Outputting mutated KT and FASTA both work.
KT outputs segment indexing, SV history, and final KT.
KT output thoroughly tested.
FASTA output tested on very small scale artificial genome.

Known Problem:
SV history for translocation can be a hard to interpret with the current output.
Translocation_deletion not implemented.
FASTA output need to be tested on a larger scale.
I am aware that the FASTA outputting takes incredibly long (about 20min for all 24 Chr).