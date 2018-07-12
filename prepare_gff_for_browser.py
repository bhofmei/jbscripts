import sys, os
from name_formatting import *

# Usage: python3 prepare_gff_for_browser.gff [-rna] [-rpt] [-tes] [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] <gff_file>

def processInputs( gffFileStr, includeRNA, includeTransposons, includeRepeats, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	# decode option types
	chrmOptions = checkEmptyAsIs( chrmOptions, 'chrms' )
	cCap, cLong, cUnSc, cZero, cEmpty, cAsIs = decodeChrmOptions( chrmOptions )
	print( 'Chromosome formatting: {:s}'.format('None' if cCap == None else (formatChrm('1', cCap, cLong, cUnSc, cZero, cEmpty, cAsIs )) ))
	scafOptions = checkEmptyAsIs( scafOptions, 'scaffolds' )
	sCap, sShort, sUnSc, sZero, sEmpty, sAsIs = decodeScafOptions( scafOptions )
	print( 'Scaffold formatting: {:s}'.format('None' if sCap == None else (formatScaf('1', sCap, sShort, sUnSc, sZero, sEmpty, sAsIs )) ) )
	contigOptions = checkEmptyAsIs( contigOptions, 'contigs' )
	tCap, tUnSc, tZero, tEmpty, tAsIs = decodeContigOptions( contigOptions )
	print( 'Contig formatting: {:s}'.format( 'None' if tCap == None else (formatContig('1', tCap, tUnSc, tZero, tEmpty, tAsIs )) ) )
	oCap, oLower, oChrm = decodeOtherOptions( otherOptions )
	print( 'Other formatting: {:s}'.format( 'None' if oCap == None else  (formatOther('Other', oCap, oLower, oChrm ) ) ) )
	mtType, chType, lmType = decodeCLMOptions( clmOptions )
	print( 'Mitochondria format: {:s}\nChloroplast format: {:s}\nLambda format: {:s}'.format( ('None' if mtType == None else mtType ), ('None' if chType == None else chType ), ('None' if lmType == None else lmType ) ) )
	
	print( 'Reading', os.path.basename(gffFileStr) )
	
	genes, rna, tes, repeats = readGFF( gffFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )
	fileName = os.path.basename( gffFileStr )
	fileDir = os.path.dirname( gffFileStr )
	rInd = fileName.rfind( '.' )
	outName = fileName[:rInd]
	outPre = os.path.join( fileDir,  outName )
	outFileGenes = outPre + '_browser_genes.gff'
	print( 'Writing', os.path.basename(outFileGenes) )
	writeOutput( outFileGenes, genes )
	
	if includeRNA:
		outFileRNA = outPre + '_browser_rna.gff'
		print( 'Writing', os.path.basename(outFileRNA) )
		writeOutput( outFileRNA, rna )
		
	if includeTransposons:
		outFileTE = outPre + '_browser_transposons.gff'
		print( 'Writing', os.path.basename(outFileTE) )
		writeOutput( outFileTE, tes )
		
	if includeRepeats:
		outFileRP = outPre + '_browser_repeats.gff'
		print( 'Writing', os.path.basename(outFileRP) )
		writeOutput( outFileRP, repeats )
	print( 'Done' )

def readGFF( gffFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	gffFile = open( gffFileStr, 'r' )
	genes = ''
	rna = ''
	te = ''
	repeat = ''
	formatDict = {}
	
	currGeneLine = None
	writeGene = False
	writeRNA = False
	
	for line in gffFile:
		# headers and blank lines
		if line.startswith( '#' ):
			continue
		else:
			lineAr = line.rstrip().split('\t')
			if len(lineAr) < 8:
				continue
			# GFF: (0) chrom (1) source (2) feature (3) start (4) end (5) score
			# (6) strand (7) frame/phase (8) attributes
			# don't want chromosomes or exons
			if lineAr[2] in [ 'chromosome', 'protein', 'intron' ]:
				continue
			# format chrm name
			if lineAr[0] in formatDict.keys():
				name = formatDict[lineAr[0]]
			else:
				name = formatChrmName( lineAr[0], chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )
				formatDict[lineAr[0]] = name
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
			# rest of gene stuff -> don't include exons
			elif writeGene and lineAr[2] != 'exon':
				genes += line
			# possible rest of rna stuff
			elif writeRNA:
				rna += line
	# end for line
	gffFile.close()
	return genes, rna, te, repeat

def formatChrmName( name, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	
	nType = determineType( name )
	if nType == 'chr':
		cCap, cLong, cUnSc, cZero, cEmpty, cAsIs = decodeChrmOptions( chrmOptions )
		nname = formatChrm( name, cCap, cLong, cUnSc, cZero, cEmpty, cAsIs )
	elif nType == 'scaf':
		sCap, sShort, sUnSc, sZero, sEmpty, sAsIs = decodeScafOptions( scafOptions )
		nname = formatScaf( name, sCap, sShort, sUnSc, sZero, sEmpty, sAsIs )
	elif nType == 'contig':
		tCap, tUnSc, tZero, tEmpty, tAsIs = decodeContigOptions( contigOptions )
		nname = formatContig( name, tCap, tUnSc, tZero, tEmpty, tAsIs )
	elif nType == 'other':
		oCap, oLower, oChrm = decodeOtherOptions( otherOptions )
		nname = formatOther( name, oCap, oLower, oChrm )
	elif nType ==  'clm':
		mtType, chType, lmType = decodeCLMOptions( clmOptions )
		nname = formatCLM( name, mtType, chType, lmType )
	return (False if nname == None else nname )

def writeOutput( outFileStr, outStr ):
	outFile = open( outFileStr, 'w' )
	outFile.write( outStr )
	outFile.close()
	
def parseInputs( argv ):
	chrmOptions = ''
	scafOptions = ''
	contigOptions = ''
	clmOptions = ''
	otherOptions = ''
	includeRNA = False
	includeTransposons = False
	includeRepeats = False
	startInd = 0
	
	for i in range(min(9,len(argv)-1)):
		if argv[i] == '-rna':
			includeRNA = True
			startInd += 1
		elif argv[i] == '-tes':
			includeTransposons = True
			startInd += 1
		elif argv[i] == '-rpt':
			includeRepeats = True
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			chrmOptions = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-s=' ):
			scafOptions = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-t=' ):
			contigOptions = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-l=' ):
			clmOptions = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-o=' ):
			otherOptions = argv[i][3:]
			startInd += 1
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv.startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	gffFileStr = argv[startInd]
	processInputs( gffFileStr, includeRNA, includeTransposons, includeRepeats, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )

def printHelp():
	print ("Usage: python3 prepare_gff_for_browser.gff [-rna] [-rpt] [-tes] [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] <gff_file>")
	print( 'Renames chromosomes/contigs based on input parameters\nRemoves unnecessary annotations so the browser is cleaner' )
	print( 'Required:' )
	print( 'gff_file\tGFF formatted file to be processed' )
	print( 'Optional:' )
	print( '-rna\tGFF file has non-coding RNA annotated; writes these to a separate file' )
	print( '-tes\tGFF file has transposons annotated; writes these to a separate file' )
	print( '-rpt\tGFF file has repeats annotated; writes these to a separate file' )
	print( '\t\trepeat is considered "similarity" or "RM"' )
	print( 'Formatting: ')
	print( getFormattingScheme( ) )

if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
