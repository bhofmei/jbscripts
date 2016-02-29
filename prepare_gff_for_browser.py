import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python prepare_gff_for_browser.py [-no-clean] [-rna] [-tes] [-rpt] [-o=output_prefix] [-no-scaf] [-no-clm] <gff_file>

def processInputs( gffFileStr, outPre, includeRNA, includeTransposons, includeRepeats, useScaffolds, useCLM, needClean ):
	
	print( 'Reading {:s}...'.format( gffFileStr ) )
	genes, rna, te, repeats = processGFF( gffFileStr, useScaffolds, useCLM, needClean  )
	
	if outPre == None:
		outBase = os.path.basename( gffFileStr )
		rInd = outBase.rfind ('.')
		outPre = outBase[:rInd]
	outFileGenes = outPre + '_browser_genes.gff'
	print( 'Writing {:s}...'.format( outFileGenes ) )
	writeOutput( outFileGenes, genes )
	
	if includeRNA:
		outFileRNA = outPre + '_browser_rna.gff'
		print( 'Writing {:s}...'.format( outFileRNA ) )
		writeOutput( outFileRNA, rna )
		
	if includeTransposons:
		outFileTE = outPre + '_browser_transposons.gff'
		print( 'Writing {:s}...'.format( outFileTE ) )
		writeOutput( outFileTE, te )
		
	if includeRepeats:
		outFileRP = outPre + '_browser_repeats.gff'
		print( 'Writing {:s}...'.format( outFileRP ) )
		writeOutput( outFileRP, repeats )
	print( 'Done' )

def processGFF( gffFileStr, useScaffolds, useCLM, needClean  ):
	
	gffFile = open( gffFileStr, 'r' )
	genes = ''
	rna = ''
	te = ''
	repeat = ''
	
	currGeneLine = None
	writeGene = False
	writeRNA = False
	
	for line in gffFile:
		# headers and blank lines
		if line.startswith( '#' ):
			continue
		else:
			lineAr = line.rstrip().split('\t')
			# GFF: (0) chrom (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) frame/phase (8) attributes
			# don't want chromosomes or exons
			if lineAr[2] in [ 'chromosome', 'exon', 'protein', 'intron' ]:
				continue
			# format chrm name
			name = formatChrmName( lineAr[0], useScaffolds, useCLM, needClean )
			#print( lineAr[0], name )
			if name == False:	# don't include
				continue
			else:
				lineAr[0] = name
				line = '\t'.join( lineAr ) + '\n'
			# transposons
			if lineAr[2] in [ 'transposable_element', 'transposon_fragment', 'transposable_element_gene' ]:
				writeGene = False
				writeRNA = False
				te += line
			# repeats
			# transposons
			if lineAr[2] in ['similarity','RM'] or 'repeat' in lineAr[2]:
				writeGene = False
				writeRNA = False
				repeat += line
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
	return genes, rna, te, repeat

def formatChrmName( inName, useScaffolds, useCLM, needClean ):
	
	chrm = inName.lower()
	if chrm == 'c' or chrm == 'chloroplast' or chrm == 'chrc':
		if useCLM:
			return ( 'ChrC' if needClean else chrm )
		return False
	elif chrm == 'm' or chrm == 'mitochondria' or chrm == 'chrm' or chrm == 'mt':
		if useCLM:
			return ( 'ChrM' if needClean else chrm )
		return False
	elif chrm == 'l' or chrm == 'lambda' or chrm == 'chrl':
		if useCLM:
			return ( 'ChrL' if needClean else chrm )
		return False
	# digit only number
	elif chrm.isdigit():
		return ('Chr'+chrm if needClean else chrm )
	elif chrm.startswith('chr') and needClean == False:
		return chrm
	elif chrm.startswith('chr'):
		if chrm.startswith( 'chromosome0' ):
			return chrm.replace( 'chromosome0', 'Chr' )
		elif chrm.startswith( 'chromosome' ):
			return chrm.replace( 'chromosome', 'Chr' )
		elif chrm.startswith( 'chrm0' ):
			return chrm.replace( 'chrm0', 'Chr' )
		elif chrm.startswith( 'chrm' ):
			return chrm.replace( 'chrm', 'Chr' )
		elif chrm.startswith( 'chr0' ):
			return chrm.replace( 'chr0', 'Chr' )
		elif chrm.startswith( 'chr' ):
			return chrm.replace( 'chr', 'Chr' )
	elif chrm.startswith( 'scaffold' ):
		if useScaffolds:
			return chrm
		return False
	elif chrm.startswith( 'contig' ):
		if useScaffolds:
			return chrm
		return False
	else:
		return chrm

def writeOutput( outFileStr, outStr ):
	outFile = open( outFileStr, 'w' )
	outFile.write( outStr )
	outFile.close()

def parseInputs( argv ):
	includeRNA = False
	includeTransposons = False
	includeRepeats = False
	useScaffolds = True
	useCLM = True
	needClean = True
	outPre = None
	startInd = 0
	
	for i in range(min(6,len(argv))):
		if argv[i] == '-rna':
			includeRNA = True
			startInd += 1
		elif argv[i] == '-tes':
			includeTransposons = True
			startInd += 1
		elif argv[i] == '-rpt':
			includeRepeats = True
			startInd += 1
		elif argv[i] == '-no-clean':
			needClean = False
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			outPre = argv[i][3:]
			startInd += 1
		elif argv[i] == '-no-scaf' or argv[i] == '-no-scaff':
			useScaffolds = False
			startInd +=1
		elif argv[i] == '-no-clm':
			useCLM = False
			startInd +=1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp( )
			exit()
		elif argv[i].startswith('-'):
			print( 'Error: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	
	gffFileStr = argv[startInd]
	processInputs( gffFileStr, outPre, includeRNA, includeTransposons, includeRepeats, useScaffolds, useCLM, needClean  )
			
def printHelp():
	print ("Usage: python prepare_gff_for_browser.py [-rna] [-tes] [-rpt] [-o=output_prefix] [-no-scaf] [-no-clm] <gff_file>")
	print( 'Converts an existing GFF file to be optimum for JBrowse' )
	print( 'Additionally correct chromosome naming scheme to be consistent' )
	print( 'Required:' )
	print( 'gff_file\tGFF formatted file to be processed' )
	print( 'Optional:' )
	print( '-rna\tGFF file has non-coding RNA annotated; writes these to a separate file' )
	print( '-tes\tGFF file has transposons annotated; writes these to a separate file' )
	print( '-rpt\tGFF file has repeats annotated; writes these to a separate file' )
	print( '\t\trepeat is considered "similarity" or "RM"' )
	print( '-no-clean\tdo not rename chromosomes; not recommended' )
	print( '-no-scaf\tdoes not include scaffolds/contigs in the output' )
	print( '\t\tuse when majority of DNA is in chromosomes' )
	print( '-no-clm\t\tdo not inlclude chroloplast, mitochondria, lambda, and\n\t\tnon-digit chrms in output' )
	print( '-o=out_prefix\tprefix for output GFF file(s) [default: GFF file name' )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
