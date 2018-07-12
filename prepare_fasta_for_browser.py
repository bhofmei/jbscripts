import sys, os
from name_formatting import *

# Usage: python3 prepare_fasta_for_browser.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] [-i=include_chrm_list] [-x=exclude_chrm_list] <input_fasta>

def processInputs( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions, includeList, excludeList ):
	fileName = os.path.basename( fastaFileStr )
	fileDir = os.path.dirname( fastaFileStr )
	rInd = fileName.rfind( '.' )
	outName = fileName[:rInd] + '_browser.fa'
	
	# Check include/exclude list
	inex = checkIncludeExclude( includeList, excludeList )
	if inex == False:
		exit()
	
	outFileStr = os.path.join( fileDir,  outName )
	print( 'FASTA file: {:s}\nOutput file: {:s}'.format( fileName, outName ) )
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
	
	print( 'Reading', os.path.basename(fastaFileStr) )
	outDict = readFasta( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions, includeList, excludeList )
	
	print( 'Writing', os.path.basename(outFileStr) )
	writeOutput( outFileStr, outDict )
	print( 'Done' )

def checkIncludeExclude( includeList, excludeList ):
	if includeList == None and excludeList == None:
		return True
	elif excludeList == None:	# include only
		print( 'Include list: {:s}'.format( ', '.join( includeList) ) )
		return True
	elif includeList == None:	# exclude only
		print( 'Exclude list: {:s}'.format( ', '.join( excludeList) ) )
		return True
	else:	# both specified
		for i in includeList:
			if i in excludeList:
				print( 'ERROR: {:s} is in include list and exclude list'.format( i) )
				return False
		print( 'Include list: {:s}'.format( ', '.join( includeList) ) )
		print( 'Exclude list: {:s}'.format( ', '.join( excludeList) ) )


def readFasta( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions, includeList, excludeList ):
	fastaFile = open( fastaFileStr, 'r' )
	# Options
	cCap, cLong, cUnSc, cZero, cEmpty, cAsIs = decodeChrmOptions( chrmOptions )
	sCap, sShort, sUnSc, sZero, sEmpty, sAsIs = decodeScafOptions( scafOptions )
	tCap, tUnSc, tZero, tEmpty, tAsIs = decodeContigOptions( contigOptions )
	mtType, chType, lmType = decodeCLMOptions( clmOptions )
	oCap, oLower, oChrm = decodeOtherOptions( otherOptions )
	
	includes = (includeList != None )
	excludes = (excludeList != None )
	outDict = {}
	curSeq = None
	curDigit = False
	
	for line in fastaFile:
		# new
		if line.startswith( '>' ):
			if curDigit:
				# check already in
				if outDict.get( curDigit ) != None:
					print( 'ERROR: {:s} listed twice'.format( curDigit ) )
					exit()
				outDict[curDigit] = curSeq
				curSeq, curDigit = None, False
			# new line
			lineAr = line.rstrip().split()
			name = lineAr[0].replace('>', '')
			# do we want to include this
			if includes and name not in includeList:
				continue
			elif excludes and name in excludeList:
				continue
			nType = determineType( name )
			if nType == 'chr':
				nname = formatChrm( name, cCap, cLong, cUnSc, cZero, cEmpty, cAsIs )
			elif nType == 'scaf':
				nname = formatScaf( name, sCap, sShort, sUnSc, sZero, sEmpty, sAsIs )
			elif nType == 'contig':
				nname = formatContig( name, tCap, tUnSc, tZero, tEmpty, tAsIs )
			elif nType == 'other':
				nname = formatOther( name, oCap, oLower, oChrm )
			elif nType ==  'clm':
				nname = formatCLM( name, mtType, chType, lmType )
			if nname != None:
				curDigit = formatForDict( nname, nType )
				curSeq = '>'+nname+' ' + ' '.join(lineAr[1:]) + '\n'
		else:
			# sequence line
			if curDigit:
				curSeq += line
	# end for line
	if curDigit:
		# check already in
		if outDict.get( curDigit ) != None:
			print( 'ERROR: {:s} listed twice'.format( curDigit ) )
			exit()
		outDict[curDigit] = curSeq
	fastaFile.close()
	return outDict
		
def writeOutput( outFileStr, outDict ):
	outFile = open( outFileStr, 'wt' )
	
	for key in sorted(outDict.keys()):
		outFile.write(outDict[key])
	outFile.close()
	
			
def parseInputs( argv ):
	chrmOptions = ''
	scafOptions = ''
	contigOptions = ''
	clmOptions = ''
	otherOptions = ''
	includeList = None
	excludeList = None
	startInd = 0
	
	for i in range(min(8,len(argv)-1)):
		if argv[i].startswith( '-c=' ):
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
		elif argv[i].startswith( '-i=' ):
			includeList = argv[i][3:].split( ',' )
			startInd += 1
		elif argv[i].startswith( '-x=' ):
			excludeList = argv[i][3:].split( ',' )
			startInd += 1
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv.startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	fastaFileStr = argv[startInd]
	processInputs( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions, includeList, excludeList )

def printHelp():
	print ("Usage: python3 prepare_fasta_for_browser.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] [-i=chrms_to_include] [-x=chrms_to_exclude] <input_fasta>")
	print( 'Renames the chromsomes/contigs based on input parameter formatting' )
	print( 'Required:' )
	print( 'fasta_file\tpath to FASTA formatted genomic DNA file' )
	print( 'Optional:' )
	print( '-i=chrm_list\tonly include seqs in this list, comma-separated' )
	print( '-x=chrm_list\texclude seqs in this list, comma-separated' )
	print( '\t\tnames must match seq names in input fasta' )
	print( 'Formatting: ')
	print( getFormattingScheme( ) )

if __name__ == "__main__":
	if len(sys.argv) < 2:
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
