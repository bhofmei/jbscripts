import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from bioFiles import *

# Usage: file_to_bigwig_pe.py [-keep] [-no-clean] [-strand] [-scale] [-union] [-p=num_proc]  <chrm_file> <bam_file | bed_file> [bam_file | bed_file]*

NUMPROC=1

def processInputs( chrmFileStr, bedFileStrAr, keepTmp, isStrand, isScale, isUnion, numProc, needClean ):
	if len(bedFileStrAr)==1:
		print( 'Input File: {:s}'.format( bedFileStrAr[0] ) )
	else:
		print( 'Input Files: {:s}'.format( ', '.join( bedFileStrAr ) ) )
	print('Chromosome sizes file: {:s}\nKeep temporary files: {:s}\nStrand-specific: {:s}\nScale by library size: {:s}'.format( chrmFileStr, str( keepTmp ), str( isStrand ), str( isScale) ) )
	if isStrand:
		print( 'Combine strand-specific output: {:s}\n'.format( str(isUnion) ) )
	print( 'Cleaning chromosome names: {:s}'.format( str(needClean)))
	# read chrm sizes file to get chrm list to include
	if needClean:
		print( 'Reading chromosome sizes file' )
		chrmList = readChrmFile( chrmFileStr )
	else:
		chrmList = None
	
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(bedFileStr, chrmFileStr, chrmList, keepTmp, isStrand, isScale, isUnion) ) for bedFileStr in bedFileStrAr ]
	suc = [ p.get() for p in results ]
	print( 'Done.' )
	
def readChrmFile( inFileStr ):
	
	inFile = open( inFileStr, 'r' )
	chrmList= []
	
	for line in inFile:
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) length ...
		chrmList += [ lineAr[0] ]
	inFile.close()
	return chrmList
	
def processFile( bedFileStr, chrmFileStr, chrmList, keepTmp, isStrand, isScale, isUnion ):
	baseDir = os.path.dirname( bedFileStr )
	ind = bedFileStr.rfind('.')
	baseName = bedFileStr[:ind]
	fileBase = os.path.basename( bedFileStr )
	rmFile = []
	
	# check for bam file -> convert to bed
	if bedFileStr.endswith( '.bam' ):
		print( 'Converting {:s} to bed'.format( fileBase ) )
		bamFileStr = bedFileStr
		bedFileStr = '{:s}.bed'.format( baseName)
		command = 'bedtools bamtobed -i {:s} > {:s}'.format( bamFileStr, bedFileStr )
		subprocess.call( command, shell=True )
		rmFile += [ bedFileStr ]
	
	if chrmList != None:
		# check chrm names and such
		bedOutStr = '{:s}_clean.bed'.format( baseName )
		print( 'Checking chromosome names in bed' )
		bedFileStr = cleanBedFile( bedFileStr, bedOutStr, chrmList )
		rmFile += [ bedFileStr ]
		
	# get scale value if necessary
	scaleVal = 1
	if isScale:
		scaleVal = getScaleValue( bedFileStr )
	
	# sort bedfile
	sortBedFile( bedFileStr )
	
	# bed to bedGraph
	if isUnion:
		bedGraphFileAr = convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal )
		rmFile += bedGraphFileAr
		bedGraphFile = uniteBedGraph( bedGraphFileAr[0], bedGraphFileAr[1], chrmFileStr, baseName )
		subAr = [ '_union' ]
	elif isStrand:
		bedGraphFile = convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal )
		subAr = [ '_plus', '_minus' ]
	else:
		bedGraphFile = convertToBedGraph( bedFileStr, chrmFileStr, baseName, scaleVal )
		subAr = [ '' ]
	
	rmFile += bedGraphFile
	
	# convert to bigwig
	bigWigFile = convertToBigWig( bedGraphFile, chrmFileStr, baseName, subAr )
	
	if keepTmp == False:
		print( 'Removing temporary files' )
		for f in rmFile:
			os.remove(f)	
	
	print( 'Output written to {:s}'.format( ' & '.join(bigWigFile) ) )

def cleanBedFile( bedFileStr, outFileStr, chrmList ):
	'''
		simply going through the BED file cleaning up chrm names and
		removing reads that aren't part of chrmList
	'''
	bedFile = open( bedFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in bedFile:
		if line.startswith( '#' ):
			outFile.write( line )
			continue
		lineAr = line.rstrip().split( '\t' )
		name = formatChrmName( lineAr[0] )
		if name in chrmList:
			lineAr[0] = name
			outFile.write( '\t'.join( lineAr ) + '\n' )
	# end for
	outFile.close()
	bedFile.close()
	return outFileStr

def formatChrmName( inName ):
	
	chrm = inName.lower()
	if chrm == 'c' or chrm == 'chloroplast' or chrm == 'chrc':
		'ChrC'
	elif chrm == 'm' or chrm == 'mitochondria' or chrm == 'chrm' or chrm == 'mt':
		return 'ChrM'
	elif chrm == 'l' or chrm == 'lambda' or chrm == 'chrl':
		return 'ChrL'
	# digit only number
	elif chrm.isdigit():
		return 'Chr'+chrm
	elif chrm.startswith( 'chromosome' ):
		return chrm.replace( 'chromosome', 'Chr' )
	elif chrm.startswith( 'chrm' ):
		return chrm.replace( 'chrm', 'Chr' )
	elif chrm.startswith( 'chr' ):
		return chrm.replace( 'chr', 'Chr' )
	elif chrm.startswith( 'scaffold' ):
		return chrm
	elif chrm.startswith( 'contig' ):
		return chrm
	else:
		return chrm


def getScaleValue( bedFileStr ):
	command = 'wc -l {:s}'.format( bedFileStr )
	readCountStr = subprocess.check_output( command, shell=True, universal_newlines = True )
	return float( readCountStr.split()[0] ) / 1000000
	
def sortBedFile( bedFileStr ):
	print( 'Sorting bed file' )
	command = 'bedSort {:s} {:s}'.format( bedFileStr, bedFileStr )
	subprocess.call( command, shell=True )
	
def convertToBedGraph( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating bedGraph of {:s}'.format( os.path.basename(bedFileStr) ) )
	bedGraphFile = '{:s}.bedGraph'.format( baseName )
	command = 'bedtools genomecov -bga -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFile )
	subprocess.call( command, shell=True )
	return [bedGraphFile]

def convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating strand-specific bedGraph of {:s}'.format( os.path.basename(bedFileStr) ) )
	bedGraphFilePlus = '{:s}_scale_plus.bedGraph'.format( baseName )
	bedGraphFileMinus = '{:s}_scale_minus.bedGraph'.format( baseName )
	
	command = 'bedtools genomecov -bga -strand + -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFilePlus )
	subprocess.call( command, shell=True )
	command = 'bedtools genomecov -bga -strand - -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFileMinus+'.tmp' )
	subprocess.call( command, shell=True )
	# negative strand -> negative scores
	negativeBedGraph( bedGraphFileMinus, bedGraphFileMinus + '.tmp' )
	os.remove( bedGraphFileMinus + '.tmp' )
	return [bedGraphFilePlus, bedGraphFileMinus]

def negativeBedGraph( bedGraphFileMinus, bedGraphTmp ):
	
	tmpFile = open( bedGraphTmp, 'r' )
	bedGraphFile = open( bedGraphFileMinus, 'w' )
	
	for line in tmpFile:
		lineAr = line.rstrip().split( '\t' )
		if len(lineAr) > 3:
			lineAr[3] = '-'+lineAr[3]
		bedGraphFile.write( '{:s}\n'.format( '\t'.join(lineAr) ) )
	tmpFile.close()
	bedGraphFile.close()

def uniteBedGraph( bedGraphPlus, bedGraphMinus, chrmFileStr, baseName ):
	print( 'Uniting bedGraphs {:s} & {:s}'.format( bedGraphPlus, bedGraphMinus ) )
	bedGraphFile = '{:s}_stranded.bedGraph'.format( baseName )
	command = 'bedtools unionbedg -empty -g {:s} -i {:s} {:s} > {:s}'.format( chrmFileStr, bedGraphPlus, bedGraphMinus, bedGraphFile + '.tmp' )
	subprocess.call( command, shell=True )
	sumBedGraph( bedGraphFile, bedGraphFile + '.tmp' )
	os.remove( bedGraphFile + '.tmp' )
	return [ bedGraphFile ]
	
def sumBedGraph( bedGraphFile, bedGraphTmp ):
	
	bedGraph = open( bedGraphFile, 'w' )
	tmpFile = open( bedGraphTmp, 'r' )
	for line in tmpFile:
		lineAr = line.rstrip().split('\t' )
		# (0) chrm (1) start (2) end (3+) samples
		lNum = [ float(x) for x in lineAr[3:] ]
		bedGraph.write( '{:s}\t{:.2f}\n'.format( '\t'.join(lineAr[0:3]), sum(lNum) ) )
	tmpFile.close()
	bedGraph.close()
	
def convertToBigWig( bedGraphFileAr, chrmFileStr, baseName, subAr ):
	
	print( 'Converting {:s} to bigWig'.format( '&'.join( bedGraphFileAr) ) )
	bwFiles = [ '{:s}{:s}.bw'.format(baseName, x) for x in subAr ]
	for i in range(len(bedGraphFileAr)):
		command = 'bedGraphToBigWig {:s} {:s} {:s}'.format( bedGraphFileAr[i], chrmFileStr, bwFiles[i] )
		subprocess.call( command, shell=True )
	return bwFiles

def parseInputs( argv ):
	keepTmp = False
	isStrand = False
	isScale = False
	isUnion = False
	needClean = True
	numProc = NUMPROC
	startInd = 0
	for i in range(min(7,len(argv)-2)):
		if argv[i] == '-keep':
			keepTmp = True
			startInd += 1
		elif argv[i] == '-no-clean':
			needClean = False
			startInd += 1
		elif argv[i] == '-strand':
			isStrand = True
			startInd += 1
		elif argv[i] == '-scale':
			isScale = True
			startInd += 1
		elif argv[i] == '-union':
			isUnion = True
			startInd += 1
		elif argv[i].startswith( '-p=' ):
			try:
				numProc = int( argv[i][3:] )
				startInd += 1
			except ValueError:
				print( 'ERROR: number of processors must be integer' )
				exit()
		elif argv[i] in [ '-h', '--help', '-help']:
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid option'.format( argv[i] ) )
			exit()
	
	# can't have union and not strand
	if isUnion and isStrand == False:
		print( 'ERROR: union should only be used with strand' )
		exit()
		
	chrmFileStr = argv[startInd]
	bedFileStrAr = []
	for j in range(startInd+1, len(argv)):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( chrmFileStr, bedFileStrAr, keepTmp, isStrand, isScale, isUnion, numProc, needClean )
	
def printHelp():
	print ("Usage: python3 file_to_bigwig_pe.py [-keep] [-no-clean] [-strand] [-scale] [-union] [-p=num_proc] <chrm_file> <bam_file | bed_file> [bam_file | bed_file]*")
	print( 'Note: bedtools and bedGraphToBigWig programs must be in the path' )
	print( 'Required:' )
	print( 'chrm_file\ttab-delimited file with chromosome names and lengths, i.e. fasta index file' )
	print( 'bam_file\tbam file that already has been indexed, i.e. file.bam.bai' )
	print( 'bed_file\tBED formatted file' )
	print( 'Optional:' )
	print( '-keep\t\tkeep intermediate files' )
	print( '-no-clean\tdoes not check chromosome names match chrm file\n\t\tnot recommended\n' )
	print( '-strand\t\tseparate reads by strand to have strand-specific bigwig files' )
	print( '-scale\t\tscale the bigwig values by total number of reads in file' )
	print( '-union\t\tused with strand, add stand-specific values for one bigwig file' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
