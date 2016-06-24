import sys, os
from name_formatting import *

# Usage: python3 prepare_fasta_for_browser.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] <input_fasta>

def processInputs( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	fileName = os.path.basename( fastaFileStr )
	fileDir = os.path.dirname( fastaFileStr )
	rInd = fileName.rfind( '.' )
	outName = fileName[:rInd] + '_browser.fa'
	outFileStr = os.path.join( fileDir,  outName )
	print( 'FASTA file: {:s}\nOutput file: {:s}'.format( fileName, outName ) )
	# decode option types
	chrmOptions = checkEmpty( chrmOptions, 'chrms' )
	cCap, cLong, cUnSc, cZero, cEmpty = decodeChrmOptions( chrmOptions )
	print( 'Chromosome formatting: {:s}'.format( formatChrm('1', cCap, cLong, cUnSc, cZero, cEmpty ) ) )
	scafOptions = checkEmpty( scafOptions, 'scaffolds' )
	sCap, sShort, sUnSc, sZero, sEmpty = decodeScafOptions( scafOptions )
	print( 'Scaffold formatting: {:s}'.format( formatScaf('1', sCap, sShort, sUnSc, sZero, sEmpty ) ) )
	contigOptions = checkEmpty( contigOptions, 'contigs' )
	tCap, tUnSc, tZero, tEmpty = decodeContigOptions( contigOptions )
	print( 'Contig formatting: {:s}'.format( formatContig('1', tCap, tUnSc, tZero, tEmpty ) ) )
	oCap, oLower, oChrm = decodeOtherOptions( otherOptions )
	print( 'Other formatting: {:s}'.format( formatOther('Other', oCap, oLower, oChrm ) ) )
	mtType, chType, lmType = decodeCLMOptions( clmOptions )
	print( 'Mitochondria format: {:s}\nChloroplast format: {:s}\nLambda format: {:s}'.format( ('None' if mtType == None else mtType ), ('None' if chType == None else chType ), ('None' if lmType == None else lmType ) ) )
	
	print( 'Reading', os.path.basename(fastaFileStr) )
	outDict = readFasta( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )
	
	print( 'Writing', os.path.basename(outFileStr) )
	writeOutput( outFileStr, outDict )
	print( 'Done' )

def readFasta( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions ):
	fastaFile = open( fastaFileStr, 'r' )
	# Options
	cCap, cLong, cUnSc, cZero, cEmpty = decodeChrmOptions( chrmOptions )
	sCap, sShort, sUnSc, sZero, sEmpty = decodeScafOptions( scafOptions )
	tCap, tUnSc, tZero, tEmpty = decodeContigOptions( contigOptions )
	mtType, chType, lmType = decodeCLMOptions( clmOptions )
	oCap, oLower, oChrm = decodeOtherOptions( otherOptions )
	
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
			nType = determineType( name )
			if nType == 'chr':
				nname = formatChrm( name, cCap, cLong, cUnSc, cZero, cEmpty )
			elif nType == 'scaf':
				nname = formatScaf( name, sCap, sShort, sUnSc, sZero, sEmpty )
			elif nType == 'contig':
				nname = formatContig( name, tCap, tUnSc, tZero, tEmpty )
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
	startInd = 0
	
	for i in range(min(6,len(argv)-1)):
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
		elif argv[i] in ['-h','--h','--help','-help']:
			printHelp()
			exit()
		elif argv.startswith( '-' ):
			print( '{:s} is not a valid parameter'.format( argv[i] ) )
			exit()
	# end for
	fastaFileStr = argv[startInd]
	processInputs( fastaFileStr, chrmOptions, scafOptions, contigOptions, clmOptions, otherOptions )

def printHelp():
	print ("Usage: python3 prepare_fasta_for_browser.py [-c=chrm_opts] [-s=scaf_opts] [-t=contig_opts] [-l=clm_opts] [-o=other_opts] <input_fasta>")
	print( 'Required:' )
	print( 'fasta_file\tpath to FASTA formatted genomic DNA file' )
	print( 'Formatting: ')
	print( getFormattingScheme( ) )

if __name__ == "__main__":
	if len(sys.argv) < 2:
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
