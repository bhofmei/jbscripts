import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 prepare_gff_for_browser.py [-r] [-t] [-o=output_prefix] <gff_file>

def processInputs( gffFileStr, outPre, includeRNA, includeTransposons ):
	
	print( 'Reading {:s}...'.format( gffFileStr ) )
	header, genes, rna, te = processGFF( gffFileStr )
	
	if outPre == None:
		outBase = os.path.basename( gffFileStr )
		rInd = outBase.rfind ('.')
		outPre = outBase[:rInd]
	outFileGenes = outPre + '_browser_genes.gff'
	print( 'Writing {:s}...'.format( outFileGenes ) )
	writeOutput( outFileGenes, header + genes )
	
	if includeRNA:
		outFileRNA = outPre + '_browser_rna.gff'
		print( 'Writing {:s}...'.format( outFileRNA ) )
		writeOutput( outFileRNA, header + rna )
		
	if includeTransposons:
		outFileTE = outPre + '_browser_transposons.gff'
		print( 'Writing {:s}...'.format( outFileTE ) )
		writeOutput( outFileTE, header + te )
	print( 'Done' )

def processGFF( gffFileStr ):
	
	gffFile = open( gffFileStr, 'r' )
	header = ''
	genes = ''
	rna = ''
	te = ''
	
	currGeneLine = None
	writeGene = False
	writeRNA = False
	
	for line in gffFile:
		# blank lines
		if line.rstrip() == '###':
			continue
		# headers
		elif line.startswith( '#' ):
			header += line
		# all other lines
		else:
			lineAr = line.rstrip().split('\t')
			# GFF: (0) chrom (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) frame/phase (8) attributes
			# don't want chromosomes or exons
			if lineAr[2] in [ 'chromosome', 'exon', 'protein', 'intron' ]:
				continue
			# transposons
			if lineAr[2] in [ 'transposable_element', 'transposon_fragment', 'transposable_element_gene' ]:
				writeGene = False
				writeRNA = False
				te += line
			# pseudogenes
			elif lineAr[2] in [ 'pseudogene', 'pseudogenic_transcript', 'pseudogenic_exon' ]:
				writeGene = False
				writeRNA = False
				genes += line
			# genes -> can be mRNA, ncRNA, rRNA, snRNA, snoRNA, tRNA
			elif lineAr[2] == 'gene':
				currGeneLine = line
				writeGene = False
				writeRNA = False
			# collected gene information
			elif currGeneLine != None:
				if lineAr[2] == 'mRNA':
					genes += currGeneLine + line
					currGeneLine = None
					writeGene = True
				else:
					rna += currGeneLine + line
					currGeneLine = None
					writeRNA = True
			# rest of gene stuff
			elif writeGene:
				genes += line
			# possible rest of rna stuff
			elif writeRNA:
				rna += line
	# end for line
	gffFile.close()
	return header, genes, rna, te

def writeOutput( outFileStr, outStr ):
	outFile = open( outFileStr, 'w' )
	outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	includeRNA = False
	includeTransposons = False
	outPre = None
	startInd = 0
	
	for i in range(min(3,len(argv))):
		if argv[i] == '-r':
			includeRNA = True
			startInd += 1
		elif argv[i] == '-t':
			includeTransposons = True
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-h=' ):
			printHelp( )
			exit()
		elif argv[i].startswith('-'):
			print( 'Error: {:s} is not a valid option' )
			exit()
	
	gffFileStr = argv[startInd]
	processInputs( gffFileStr, outPre, includeRNA, includeTransposons )
			
def printHelp():
	print ("Usage: python3 prepare_gff_for_browser.py [-r] [-t] [-o=output_prefix] <gff_file>")
	print( 'Converts an existing GFF file to be optimum for JBrowse' )
	print( 'Required:' )
	print( 'gff_file\tGFF formatted file to be processed' )
	print( 'Optional:' )
	print( '-r\tGFF file has non-coding RNA annotated; writes these to a separate file' )
	print( '-t\tGFF file has transposons annotated; writes these to a separate file' )
	print( '-o=out_prefix\tprefix for output GFF file(s) [default: GFF file name' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
