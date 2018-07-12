import sys, os
from name_formatting import *

# Usage: python3 format_file_names.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] [-d=delim] [-n=col_num] <infile>

def processInputs( inFileStr, delim, colNum, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
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
	
	fileName = os.path.basename( inFileStr )
	fileDir = os.path.dirname( inFileStr )
	rInd = fileName.rfind( '.' )
	outName = fileName[:rInd] + '_fmt'+fileName[rInd:]
	outFileStr = os.path.join( fileDir,  outName )
	
	if inFileStr.endswith( '.csv' ) and delim == None:
		delim = ','
	elif inFileStr.endswith( '.tsv' ) and delim == None:
		delim = '\t'
	elif delim == 'tab':
		delim = '\t'
	elif delim == None:
		print( 'ERROR: must specify delim when file does not end in .csv or .tsv' )
		exit()
	
	print( 'Reading', os.path.basename(inFileStr) )
	print( 'Writing to', os.path.basename(outFileStr) )
	readFile( inFileStr, outFileStr, delim, colNum, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )
	print( 'Done.' )
	
def readFile( inFileStr, outFileStr, delim, colNum, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	inFile = open( inFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	formatDict = {}
	
	for line in inFile:
		# headers and blank lines
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split(delim)
		if len(lineAr) < colNum:
			continue
		# format chrm name
		name = lineAr[colNum]
		if name in formatDict.keys():
			nname = formatDict[name]
		else:
			nname = formatChrmName( name, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )
			formatDict[lineAr[0]] = nname
		if nname:
			lineAr[colNum] = nname
			line = delim.join( lineAr ) + '\n'
			outFile.write( line )
	# end for line
	inFile.close()
	outFile.close()

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

def parseInputs( argv ):
	chrmOptions = ''
	scafOptions = ''
	contigOptions = ''
	clmOptions = ''
	otherOptions = ''
	delim=None
	colNum = 0
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
		elif argv[i].startswith( '-d=' ):
			delim = argv[i][3:]
			startInd += 1
		elif argv[i].startswith( '-n=' ):
			try:
				colNum = int(argv[i][3:])
				startInd += 1
			except ValueError:
				print('ERROR: column number must be integer' )
				exit()
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	inFileStr = argv[startInd]
	processInputs( inFileStr, delim, colNum, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )

def printHelp():
	print ("Usage: python3 format_file_names.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] [-d=delim] [-n=col_num] <infile>\ncol_num is 0-indexed")
	print( 'Maintains the same underlying data for input file but changes the\nchromsome names base on parameters' )
	print( 'Required:' )
	print( 'in_file\tpath to data file\n\tinfile must be .csv or .tsv unless delim is specified' )
	print( 'Optional: ' )
	print( '-d=delim\tfield deliminator' )
	print( '-n=col_num\t0-based column number for chrm column' )
	print( 'Formatting: ')
	print( getFormattingScheme( ) )
	
if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
