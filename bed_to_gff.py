import sys, math, glob, multiprocessing, subprocess, os, bisect, random

# Usage: python3.4 bed_to_gff.py [-i=peak_finder_program_name] <bed_file>

def processInputs( bedFileStr, source ):
	if bedFileStr.endswith( '.bed' ) == False:
		gffFileStr = bedFileStr + '.gff'
	else:
		rInd = bedFileStr.rfind( '.' )
		gffFileStr = bedFileStr[:rInd] + '.gff'
	print( 'Converting {:s} to {:s}'.format( bedFileStr, gffFileStr ) )
	convertFile( bedFileStr, gffFileStr, source )
	print( 'Done.' )

def convertFile( bedFileStr, gffFileStr, source ):
	bedFile = open( bedFileStr, 'r' )
	gffFile = open( gffFileStr, 'w' )
	
	for line in bedFile:
		if line.startswith( '#' ):
			continue
		lineAr = line.rstrip().split( '\t' )
		# BED: (0) chrom (1) start (2) end (3) name (4) score (5) strand ...
		# only 0-2 are needed
		chrom = lineAr[0]
		start = int( lineAr[1] ) + 1
		end = int( lineAr[2] ) + 1
		moreFields = ['.'] * 3
		name = 'NA'
		if len( lineAr ) > 3:
			name = lineAr[3]
		if len( lineAr ) > 4:
			#print( lineAr )
			for i in range(4, min(len(lineAr),6)):
				moreFields[i-4] = lineAr[i]
			
		# GFF: (0) chrom (1) source (2) feature (3) start (4) end (5) score
		# (6) strand (7) frame/phase (8) attributes
		gffLine = '{:s}\t{:s}\tpeak\t{:d}\t{:d}\t{:s}\tName={:s}\n'.format( chrom, source, start, end, '\t'.join( moreFields ), name )
		gffFile.write( gffLine )
	# end for line
	bedFile.close()
	gffFile.close()

def parseInputs( argv ):
	startInd = 0
	source = 'Unknown'
	
	for i in range( min(1, len(argv)) ):
		if argv[i].startswith( '-i=' ):
			source = argv[i][3:]
			startInd += 1
		
	bedFileStr = argv[startInd]
	
	if bedFileStr.endswith('.bed') == False:
		print( 'WARNING: bed file does not end with "bed"...check output results for correctness')
	processInputs( bedFileStr, source )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python3.4 bed_to_gff.py [-i=source] <bed_file>")
	else:
		parseInputs( sys.argv[1:] )
