# Karyotype
Require downloading your own genome, get RefSeq Version at
https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19

Put the genome in /Genomes, and Main/Main.py should be ready to go

The "test_manual.json" should work directly with Main.py
just call `Main.py test_manual.json`

Current State:
Intra-chromosomal deletion, inversion, and duplication are implemented and tested.
JSON works for manual mode only. Outputting mutated KT and FASTA both works, though
I can't test the FASTA output.

Plan: Adaptation to multiple chromosome should be quite easy as everything are done
in a general manner. A random generator will also need to be implemented to support 
Random Mode.