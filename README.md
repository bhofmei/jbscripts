#JBrowse Scripts

custom python scripts for the Schmitz lab JBrowse

all scripts can be run with Python3  
build_tracklist_json.py can be run with Python2 or Python3

##add_ortholog_gff.py
adds ortholog information to species GFF files
useful for being able to jump between orthologous genes in the browser

##allc_to_bigwig_pe.py
converts allC files to the appropriate BigWig format to be used in the browser  
Prerequistes:
bedGraphToBigWig is a program from UCSC and it needs to be on the computers path; do not include the program name, just the directory it is in
`export PATH=$PATH:directory/to/bedGraphToBigWig/`

##bed_to_gff.py
converts a BED file to GFF format  
useful for locally (in internet browser) uploading BED files, such as ChIP peaks

##build_tracklist_json.py
converts a CSV format file to the JSON format needed by JBrowse  

##file_to_bigwig_pe.py
converts BED and BAM files to BigWig format; BigWig format is used for the coverage plots  
Prerequistes:
bedGraphToBigWig is a program from UCSC and needs to be on the computer's path
bedtools is a software suite for manipulating BED files; also needs to be on the computer's path  
`export PATH=$PATH:directory/to/bedGraphToBigWig/:directory/to/bedtools/bin`

##prepare_gff_for_browser.py
cleans a GFF file to be appropraite for the browser, i.e. removes chromosome, protein, and exon features  
optionally separates the non-coding RNAs and repeats/transposable elements to separate GFF files
