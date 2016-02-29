import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3 prepare_fasta_for_browser.py [-no-scaf] [-no-clm] [-i=chrm_list] <fasta_file>
# formats the chromosome naming so it is consistent across the browser
# removes scaffolds if necessary; renames chloroplast, mitochondria, lambda

def processInputs( fastaFileStr, useScaffolds, useCLM, includeList ):
	fileName = os.path.basename( fastaFileStr )
	fileDir = os.path.dirname( fastaFileStr )
	rInd = fileName.rfind( '.' )
	outName = fileName[:rInd] + '_browser.fa'
	outFileStr = os.path.join( fileDir,  outName )
	print( 'FASTA file: {:s}\nOutput file: {:s}\nRemove scaffolds: {:s}\nRemove CLM: {:s}'.format( fileName, outName, str( (not useScaffolds) ), str( (not useCLM) ) ) )
	if includeList != None:
		print( 'Include list: {:s}'.format( ', '.join( includeList) ) )
	
	chrmDict, scafDict, clmDict = readFasta( fastaFileStr, useScaffolds, useCLM, includeList )
	print( 'Writing output to', outFileStr )
	writeOutput( outFileStr, chrmDict, scafDict, clmDict )
	print( 'Done' )

def readFasta( fastaFileStr, useScaffolds, useCLM, includeList ):
	includes = ( includeList != None )
	fastaFile = open( fastaFileStr, 'r' )
	chrmDict = {}
	scafDict = {}
	clmDict = {}
	curInfo = None
	curDigit = None
	curType = None
	getSeq = False
	
	includeSeq = False
	for line in fastaFile:
		# new chrm/scaffold
		if line.startswith( '>' ):
			# save previous info
			if getSeq:
				if curType == 'chrm':
					chrmDict = addToDict( chrmDict, curInfo, curDigit, 'chrm' )
				elif curType == 'scaf':
					scafDict = addToDict( scafDict, curInfo, curDigit, 'scaf' )
				elif curType == 'clm':
					clmDict = addToDict( clmDict, curInfo, curDigit, 'clm' )
				curInfo, curDigit, curType, getSeq = None, None, None, False
			# get new info
			lineAr = line.rstrip().split()
			chrm = lineAr[0].replace('>','').lower() # chrm name
			#print( chrm )
			if chrm == 'c' or chrm == 'chloroplast' or chrm == 'chrc':
				curInfo = '>ChrC ' + ' '.join( lineAr[1:] ) + '\n'
				curDigit = 'c'
				curType = 'clm'
			elif chrm == 'm' or chrm == 'mitochondria'or chrm == 'chrm' or chrm=='mt':
				curInfo = '>ChrM ' + ' '.join( lineAr[1:] ) + '\n'
				curDigit = 'm'
				curType = 'clm'
			elif chrm == 'l' or chrm == 'lambda' or chrm == 'chrl':
				curInfo = '>ChrL ' + ' '.join( lineAr[1:] ) + '\n'
				curDigit = 'l'
				curType = 'clm'
			# digit only number
			elif chrm.isdigit():
				curDigit = int( chrm )
				curInfo = '>Chr'+chrm+' '+' '.join( lineAr[1:] ) + '\n'
				curType = 'chrm'
			# begins chr or Chr
			elif chrm.startswith( 'chr' ):
				curType = 'chrm'
				if chrm.startswith( 'chromosome0' ):
					try:
						curDigit = int( chrm.replace('chromosome0','') )
					except ValueError:
						curDigit = chrm.replace('chromosome0', '')
						curType = 'clm'
					chrm = chrm.replace( 'chromosome0', 'Chr' )
				elif chrm.startswith( 'chromosome' ):
					try:
						curDigit = int( chrm.replace('chromosome','') )
					except ValueError:
						curDigit = chrm.replace('chromosome', '')
						curType = 'clm'
					chrm = chrm.replace( 'chromosome', 'Chr' )
				elif chrm.startswith( 'chrm0' ):
					try:
						curDigit = int( chrm.replace('chrm0','') )
					except ValueError:
						curDigit = chrm.replace('chrm0', '')
						curType = 'clm'
					chrm = chrm.replace( 'chrm0', 'Chr' )
				elif chrm.startswith( 'chrm' ):
					try:
						curDigit = int( chrm.replace('chrm','') )
					except ValueError:
						curDigit = chrm.replace('chrm', '')
						curType = 'clm'
					chrm = chrm.replace( 'chrm', 'Chr' )
				elif chrm.startswith( 'chr0' ):
					try:
						curDigit = int( chrm.replace('chr0','') )
					except ValueError:
						curDigit = chrm.replace('chr0', '')
						curType = 'clm'
					chrm = chrm.replace( 'chr0', 'Chr' )
				else:
					try:
						curDigit = int( chrm.replace('chr','') )
					except ValueError:
						curDigit = chrm.replace('chr', '')
						curType = 'clm'
					chrm = chrm.replace( 'chr', 'Chr' )
				curInfo = '>'+chrm+' '+' '.join( lineAr[1:] )+ '\n'
			elif chrm.startswith( 'scaffold' ):
				try:
					curDigit = int( chrm.replace('scaffold','').replace('_','').replace('-','') )
				except ValueError:
					curDigit = chrm.replace('scaffold', '').replace('_','').replace('-','')
				curInfo = '>'+chrm+' '+' '.join( lineAr[1:] ) + '\n'
				curType = 'scaf'
			elif chrm.startswith( 'contig' ):
				try:
					curDigit = int( chrm.replace('contig','').replace('_','').replace('-','') )
				except ValueError:
					curDigit = chrm.replace('contig', '').replace('_','').replace('-','')
				curInfo = '>'+chrm+' '+' '.join( lineAr[1:] ) + '\n'
				curType = 'scaf'
			else:
				print( 'Unknown type', chrm )
				curDigit = chrm
				curInfo = '>'+chrm+' '+' '.join( lineAr[1:] ) + '\n'
				curType = 'clm'
			
			# check types vs ones we want to include
			if includes and lineAr[0][1:] in includeList:
				getSeq = True
			# check scaffolds
			elif curType == 'scaf' and useScaffolds:
				getSeq = True
			elif curType == 'clm' and useCLM:
				getSeq = True
			elif curType == 'chrm' and not includes:
				getSeq = True
		# seq line
		else:
			if getSeq:
				curInfo += line
	# end for line
	# add last chunk
	if getSeq:
		if curType == 'chrm':
			chrmDict = addToDict( chrmDict, curInfo, curDigit, 'chrm' )
		elif curType == 'scaf':
			scafDict = addToDict( scafDict, curInfo, curDigit, 'scaf' )
		elif curType == 'clm':
			clmDict = addToDict( clmDict, curInfo, curDigit, 'clm' )
	fastaFile.close()
	return chrmDict, scafDict, clmDict

def addToDict( inDict, curInfo, curDigit, label ):
	# test
	if inDict.get(curDigit) != None:
		print( 'WARNING: already entry for {:s} in {:s} dict'.format( str(curDigit), label ) )
		return inDict
	inDict[curDigit] = curInfo
	return inDict

def writeOutput( outFileStr, chrmDict, scafDict, clmDict ):
	outFile = open( outFileStr, 'wt' )
	nChrm = len(chrmDict.keys())
	nScaf = len( scafDict.keys())
	#nClm = len( clmDict.keys() )
	#print( list(chrmDict.keys()))
	#print( sorted(list(scafDict.keys())))
	#print( sorted(list(clmDict.keys())))
	if nChrm > 0:
		# chrm
		for digit in sorted( chrmDict.keys() ):
			outFile.write( chrmDict[digit] )
		for digit in sorted( clmDict.keys() ):
			outFile.write( clmDict[digit] )
		for digit in sorted( scafDict.keys() ):
			outFile.write( scafDict[digit] )
	else:
		# scaf
		for digit in sorted( scafDict.keys() ):
			outFile.write( scafDict[digit] )
		for digit in sorted( clmDict.keys() ):
			outFile.write( clmDict[digit] )
	outFile.close()

def parseInputs( argv ):
	useScaffolds = True
	useCLM = True
	includeList = None
	startInd = 0
	
	for i in range(min(3,len(argv))):
		if argv[i] == '-no-scaf' or argv[i] == '-no-scaff':
			useScaffolds = False
			startInd +=1
		elif argv[i] == '-no-clm':
			useCLM = False
			startInd +=1
		elif argv[i].startswith( '-i=' ):
			includeList = argv[i][3:].split( ',' )
			startInd += 1
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: argument {:s} is not valid'.format( argv[i] ) )
			exit()
	# end for
	fastaFileStr = argv[startInd]
	processInputs( fastaFileStr, useScaffolds, useCLM, includeList )


def printHelp():
	print( 'Usage: python3 prepare_fasta_for_browser.py [-no-scaff] [-no-clm] [-in=chrm_to_include] <fasta_file>' )
	print( 'Cleans/formats a FASTA file to be consistent across different speicies in the browser' )
	print( 'Required:' )
	print( 'fasta_file\tpath to FASTA formatted genomic DNA file' )
	print ( 'Optional:' )
	print( '-no-scaf\tdoes not include scaffolds/contigs in the output' )
	print( '\t\tuse when majority of DNA is in chromosomes' )
	print( '-no-clm\t\tdo not inlclude chroloplast, mitochondria, lambda,\n\t\tand other non-digit chrms in output' )
	print( '-i=chrm_list\twhen specified, only includes these chromsomes in output\n\t\tcomma-separated list' )
	print( '\t\tchromosome names must exactly match those in input fasta' )
	
if __name__ == "__main__":
	if len(sys.argv) < 2 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
