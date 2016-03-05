Changes
=======

0.1.12
--------
- **bedautosql**: option `--lines` specifying the number of lines to be 
analyzed from an input file.
- **vcfeffect2bed**: report errors to an output file (useful for 
checking whether all variant effects were successfully processed).
- Test suite added for the `flanknfilter` tool.

0.1.11
--------
- Tool `flanknfilter` to filter features from a BED or VCF file by 
having *N*'s in their flanking regions.

0.1.10
------
- Tool `vcfeffect2bed` to get a BED file of genotype effects from an 
snpEff-annotated VCF file.
- Reference alleles added to output of the `snpeff2bed` tool.


0.1.9
-----
- Tool `snpeff2bed` to obtain the BED3+ file of variant effects from an
 snpEff-annotated VCF file.
- Tool `vcfgeno2bed` to obtain the BED3+ file of variant genotypes 
from a VCF file.
 
0.1.8
-----
- Tool `gff2bed` to convert a GFF3 file to the BED format.
- Option for BED and GFF3 readers to check for sorted records. 

0.1.7
-----
- Tool `snpeff2pph` to produce PolyPhen2 input files from 
snpEff-annotated VCF files.
- Tool `gff2to3` to convert a GFF2 file to the GFF3 format.
- Tool `gfftagstat` to investigate attribute tags of a GFF3 file.
- Additional routines to process GFF files.
- Subroutines sorted by name.

0.1.6
-----
- Tool `rmout2bed` to convert RepeatMasker out files to the BED format.
- Routines to process RepeatMasker out files.

0.1.5
-----
- Unified launcher `bioformats` for package tools.
- Tool `bedcolumns` to determine the numbers of BED and extra columns
 in a BED file.
- Tool `bedautosql` to build an autoSql table for a BED file.

0.1.4
-----
- Tool `fastareorder` to reorder sequences in a FASTA file added.

0.1.3
-----
- Sequence renaming tools `renameseq` and `ncbirenameseq` added.

0.1.2
-----
- Parser for VCFtools frequency count format
- Tool to detect gaps in FASTA files

0.1.1
-----
- Parsers for BED, GFF3, LAV and BLAST tabular alignment formats
- Test suite for FASTA routines

0.1.0
-----
- Random nucleotide sequence generator
- Routines to write FASTA files
- Tool to generate random FASTA files

