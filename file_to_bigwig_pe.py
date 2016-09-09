import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: file_to_bigwig_pe.py [-keep] [-scale] [-strand] [-sort] [-p=num_proc] <chrm_file> <bam_file | bed_file> [bam_file | bed_file]*

NUMPROC=1

def processInputs( chrmFileStr, bedFileStrAr, keepTmp, isScale, isStrand, isSort, numProc ):
	if len(bedFileStrAr)==1:
		print( 'Input File: {:s}'.format( bedFileStrAr[0] ) )
	else:
		print( 'Input Files: {:s}'.format( ', '.join( bedFileStrAr ) ) )
	print('Chromosome sizes file: {:s}\nKeep temporary files: {:s}\nScale by library size: {:s}\nStranded: {:s}\nSorting: {:s}'.format( chrmFileStr, str( keepTmp ), str( isScale), str(isStrand), str(isSort) ) )
	if len(bedFileStrAr) > numProc:
		numProc = len(bedFileStrAr)
	print( 'Begin processing files with {:d} processors'.format( numProc ) )
	pool = multiprocessing.Pool( processes=numProc )
	results = [ pool.apply_async( processFile, args=(bedFileStr, chrmFileStr, keepTmp, isScale, isStrand, isSort) ) for bedFileStr in bedFileStrAr ]
	suc = [ p.get() for p in results ]
	print( 'Done.' )
	
def processFile( bedFileStr, chrmFileStr, keepTmp, isScale, isStrand, isSort ):
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
		
	# get scale value if necessary
	scaleVal = 1
	if isScale:
		scaleVal = getScaleValue( bedFileStr )
	
	# bed to bedGraph
	if isStrand:
		bedGraphFileAr = convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal )
	else:	
		bedGraphFileAr = convertToBedGraph( bedFileStr, chrmFileStr, baseName, scaleVal )
	
	rmFile += bedGraphFileAr
	bigWigFileAr = []
	for bedGraphFile in bedGraphFileAr:
		# sort if necessary
		if isSort:
			sortBedFile( bedGraphFile )
		# convert to bigwig
		bigWigFile = convertToBigWig( bedGraphFile, chrmFileStr )
		bigWigFileAr += [bigWigFile]
	
	if keepTmp == False:
		print( 'Removing temporary files' )
		for f in rmFile:
			os.remove(f)	
	
	print( 'Output written to {:s}'.format( ' & '.join(bigWigFileAr) ) )

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
		if name in chrmList:
			lineAr[0] = name
			outFile.write( '\t'.join( lineAr ) + '\n' )
	# end for
	outFile.close()
	bedFile.close()
	return outFileStr

def getScaleValue( bedFileStr ):
	command = 'wc -l {:s}'.format( bedFileStr )
	readCountStr = subprocess.check_output( command, shell=True, universal_newlines = True )
	return float( readCountStr.split()[0] ) / 1000000
	
def sortBedFile( bedFileStr ):
	print( 'Sorting bedGraph', os.path.basename(bedFileStr) )
	command = 'bedSort {:s} {:s}'.format( bedFileStr, bedFileStr )
	subprocess.call( command, shell=True )
	
def convertToBedGraph( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating bedGraph of {:s}'.format( os.path.basename(bedFileStr) ) )
	bedGraphFile = '{:s}.bedGraph'.format( baseName )
	command = 'bedtools genomecov -bga -scale {:.2f} -i {:s} -g {:s} > {:s}'.format(scaleVal, bedFileStr, chrmFileStr, bedGraphFile )
	subprocess.call( command, shell=True )
	return [bedGraphFile]

def convertToBedGraphStrand( bedFileStr, chrmFileStr, baseName, scaleVal ):
	print( 'Creating strand-specific bedGraphs of {:s}'.format( os.path.basename(bedFileStr) ) )
	bedGraphFilePlus = '{:s}.bedGraph.plus'.format( baseName )
	bedGraphFileMinus = '{:s}.bedGraph.minus'.format( baseName )
	
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
	
def convertToBigWig( bedGraphStr, chrmFileStr ):
	
	print( 'Converting {:s} to bigWig'.format( os.path.basename( bedGraphStr) ) )
	bigWigStr = bedGraphStr.replace( '.bedGraph', '.bw' )
	#print( bigWigStr )
	# bedGraphToBigWig in.bedGraph chrom.sizes out.bw
	command = 'bedGraphToBigWig {:s} {:s} {:s}'.format( bedGraphStr, chrmFileStr, bigWigStr )
	subprocess.call( command, shell=True )
	return bigWigStr

def parseInputs( argv ):
	keepTmp = False
	isScale = False
	isStrand = False
	isSort = False
	numProc = NUMPROC
	startInd = 0
	for i in range(min(5,len(argv)-2)):
		if argv[i] == '-keep':
			keepTmp = True
			startInd += 1
		elif argv[i] == '-scale':
			isScale = True
			startInd += 1
		elif argv[i] == '-strand':
			isStrand = True
			startInd += 1
		elif argv[i] == '-sort':
			isSort = True
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
		
	chrmFileStr = argv[startInd]
	bedFileStrAr = []
	for j in range(startInd+1, len(argv)):
		bedFileStrAr += [ argv[j] ]
	
	processInputs( chrmFileStr, bedFileStrAr, keepTmp, isScale, isStrand, isSort, numProc )
	
def printHelp():
	print ("Usage: python3 file_to_bigwig_pe.py [-keep] [-scale] [-strand] [-sort] [-p=num_proc] <chrm_file> <bam_file | bed_file> [bam_file | bed_file]*")
	print( 'Convert BED/BAM file to bigWig format for coverage view' )
	print( 'Note: bedtools, bedGraphToBigWig, and bedSort programs must be in the path' )
	print( 'Required:' )
	print( 'chrm_file\ttab-delimited file with chromosome names and lengths\n\t\ti.e. fasta index file' )
	print( 'bam_file\tbam file that already has been indexed, i.e. file.bam.bai' )
	print( 'bed_file\tBED formatted file' )
	print( 'Optional:' )
	print( '-keep\t\tkeep intermediate files' )
	print( '-scale\t\tscale the bigwig values by total number of reads in file' )
	print( '-strand\t\toutput reads from plus and minus strand to separate files' )
	print( '-sort\t\tsort bedgraph; use when bigwig conversion fails' )
	print( '-p=num_proc\tnumber of processors to use [default 1]' )

if __name__ == "__main__":
	if len(sys.argv) < 3 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
