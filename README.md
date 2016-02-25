#JBrowse Scripts

custom python scripts for the Schmitz lab JBrowse

all scripts can be run with Python3  
build_tracklist_json.py can be run with Python2 or Python3

##add_ortholog_gff.py
adds ortholog information to species GFF files
useful for being able to jump between orthologous genes in the browser

##allc_to_bigwig_pe.py
converts allC files to the appropriate BigWig format to be used in the browser. only includes information for chromosomes listed in the chromosome sizes file  
Prerequistes:
bedGraphToBigWig is a program from UCSC and it needs to be on the computers path; do not include the program name, just the directory it is in  
`export PATH=$PATH:directory/to/bedGraphToBigWig/`

##bed_to_gff.py
converts a BED file to GFF format  
useful for locally (in internet browser) uploading BED files, such as ChIP peaks

##build_tracklist_json.py
converts a CSV or TSV format file to the JSON format needed by JBrowse  

##file_to_bigwig_pe.py
converts BED and BAM files to BigWig format; BigWig format is used for the coverage plots. 
only includes information from chromosomes listed in the chromosome sizes file  
Prerequistes:
bedGraphToBigWig is a program from UCSC and needs to be on the computer's path
bedSort is a program from UCSC and needs to be on the computer's path
bedtools is a software suite for manipulating BED files; also needs to be on the computer's path  
`export PATH=$PATH:directory/to/bedGraphToBigWig/:directory/to/bedSort/:directory/to/bedtools/bin`

##prepare_fasta_for_browser.py
cleans a FASTA file to rename the chromosomes so they are consistent within the browser
optionally, removes all scaffolds (`-no-scaf`) and/or removes chloroplast, mitochondria, and lambda (`-no-clm`). alternatively, can specify to only include some chromosomes

##prepare_gff_for_browser.py
cleans a GFF file to be appropriate for the browser, i.e. removes chromosome, protein, and exon features
renames the chromosomes/scaffolds to be consistent within the browser
optionally removes features from scaffolds or chloroplast/mitochondria chromosomes
optionally separates the non-coding RNAs, repeats, and transposable elements to separate GFF files

#Prerequisites
For file conversions, `bedGraphToBigWig` and `bedSort` need to be on the computer's path. The programs can be downloaded from <http://hgdownload.cse.ucsc.edu/admin/exe/>; chose the appropriate operating system.  
Download the two programs to your computer and change the permissions so they are directly executable
For those unfamiliar with adding programs to the computer's path, follow these steps:
1. In terminal/command-line, navigate to the directory with the two programs
2. Type `chmod 755 bedGraphToBigWig bedSort`, just to make sure the permissions are correct
3. Type `pwd`, which will return the path to the current directory
4. Type `echo 'export PATH=${PATH}:path/from/above' >> ~./bashrc`, where `path/from/above` is the path printed out by `pwd`
5. Type `source ~/.bashrc`
6. To test that the files are now on your path, type `which bedGraphToBigWig`. This should print out the path to your current directory