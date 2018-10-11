import sys, math, glob, multiprocessing, subprocess, os

# Usage: python3 add_ortholog_gff.py [-l=label1,label2] [-c=col1,col2] [-o1=outfile_species1] [-o2=outfile_species2] <ortholog_file> <species1_gff> <species2_gff>
# note: only adding the ortholog information to genes

def processInputs( orthoFileStr, gffFileStr1, gffFileStr2, outFileStr1, outFileStr2, labels, cols ):
	
	if outFileStr1 == None:
		rInd1 = gffFileStr1.rfind( '.' )
		if labels != None:
			outFileStr1 = gffFileStr1[:rInd1] + '_' + labels[1] + '-ortholog.gff'
		else:
			outFileStr1 = gffFileStr1[:rInd1] + '_ortholog.gff'
	if outFileStr2 == None:	
		rInd2 = gffFileStr2.rfind( '.' )
		if labels != None:
			outFileStr2 = gffFileStr2[:rInd2] + '_' + labels[0] + '-ortholog.gff'
		else:
			outFileStr2 = gffFileStr2[:rInd2] + '_ortholog.gff'
	
	if labels != None:
		useLabels = True
	else:
		useLabels = False
	
	if cols == None:
		cols = [0,1]
	
	print( 'Reading ortholog file and creating dictionaries...' )
	speciesNameAr, species1Dict, species2Dict = readOrthoFile( orthoFileStr, cols, useLabels )
	if useLabels:
		speciesNameAr = labels
	#print(species1Dict)
	
	print( 'Correcting GFF for species 1...' )
	# correct species 1
	readGFF( gffFileStr1, outFileStr1, species1Dict, speciesNameAr[1] )
	# correct species 2
	#print(species2Dict)
	print( 'Correcting GFF for species 2...' )
	readGFF( gffFileStr2, outFileStr2, species2Dict, speciesNameAr[0] )
	print( 'Done' )

def readOrthoFile( orthoFileStr, cols, useLabels ):
	
	speciesName = ['','']
	# species 1 genes -> species 2 genes
	species1Dict = {}
	# species 2 genes -> species 1 genes
	species2Dict = {}
	
	orthoFile = open( orthoFileStr, 'r' )
	firstLine = not useLabels
	mCol = max(cols)
	for line in orthoFile:
		lineAr = line.rstrip().split('\t')
		if len(lineAr) <= mCol or line.startswith('#'):
			continue
		if firstLine:
			firstLine = False
			rInd1 = lineAr[0].rfind( '_' )
			speciesName[0] = lineAr[cols[0]][:rInd1]
			rInd2 = lineAr[2].rfind( '_' )
			speciesName[1] = lineAr[cols[1]][:rInd2]
			continue
		# (0) gene name (1) methylation (2) ortholog (3) methylation
		# or cols[0] = gene, name cols[1] = ortholog
		gNames = ['','']
		for i in range(2):
			gNames[i] = lineAr[ cols[i] ]
			if gNames[i].startswith( 'AT' ) or gNames[i].startswith( 'Potri' ):
				rInd = gNames[i].rfind( '.' )
				gNames[i] = gNames[i][:rInd]
			elif gNames[i].startswith( 'Th' ):
				gNames[i] += '.g'
		# assign to dictionaries
		species1Dict[gNames[0]] = gNames[1]
		species2Dict[gNames[1]] = gNames[0]
	orthoFile.close()
	return speciesName, species1Dict, species2Dict

def readGFF( gffFileStr, outFileStr, orthoDict, orthoName ):
	
	gffFile = open( gffFileStr, 'r' )
	outFile = open( outFileStr, 'w' )
	
	for line in gffFile:
		if line.startswith( '#' ):
			outFile.write( line )
			continue
		lineAr = line.rstrip().split( '\t' )
		# (0) chrm (1) source (2) type (3) start (4) end (5) ?
		# (6) strand (7) ? (8) Notes
		if lineAr[2] == 'gene':
			gName = getGeneName( lineAr[8] )
			oGene = orthoDict.get( gName )
			if oGene != None:
				lineAr[8] += ';{:s}-Ortholog={:s}'.format( orthoName, oGene )
		
		# write to new gff
		outFile.write( '{:s}\n'.format( '\t'.join( lineAr ) ) )
	# end for
	gffFile.close()
	outFile.close()
		
def getGeneName (notesStr):
	search = "ID="
	index = notesStr.find(search)
	adIndex = index + len(search)
	endIndex = notesStr[adIndex:].find(';')
	if endIndex == -1:
		return notesStr[adIndex:]
	else:
		return notesStr[adIndex:endIndex+adIndex]
def parseInputs( argv ):
	labels = None
	cols = None
	outFileStr1 = None
	outFileStr2 = None
	startInd = 0
	
	for i in range(min(4,len(argv)) ):
		if argv[i].startswith( '-l=' ):
			labels = argv[i][3:].split( ',' )
			if len(labels) != 2:
				print( 'ERROR: must be only two labels' )
				exit()
			startInd += 1
		elif argv[i].startswith( '-c=' ):
			colStr = argv[i][3:].split( ',' )
			try:
				cols = [ int(x) for x in colStr ]
				if len( cols ) != 2:
					print( 'ERROR: must be only 2 columns specified' )
					exit()
				startInd += 1
			except ValueError:
				print( 'ERROR: at least one column index not an integer' )
				exit()
		elif argv[i].startswith('-o1='):
			outFileStr1 = argv[i][4:]
			startInd += 1
		elif argv[i].startswith('-o2='):
			outFileStr2 = argv[i][4:]
			startInd += 1
		elif argv[i].startswith( '-h' ):
			printHelp()
			exit()
		elif argv[i].startswith( '-' ):
			print( 'ERROR: {:s} is not a valid parameter'.format( argv[i] ) )
			exit()
			
	# end for
	orthoFileStr = argv[startInd]
	gffFileStr1 = argv[startInd+1]
	gffFileStr2 = argv[startInd+2]
	processInputs( orthoFileStr, gffFileStr1, gffFileStr2, outFileStr1, outFileStr2, labels, cols )

def printHelp():
	print ("Usage:\tpython3 add_ortholog_gff.py [-l=label1,label2] [-c=col1,col2] \n\t\t[-o1=outfile_species1] [-o2=outfile_species2] <ortholog_file> <species1_gff>\n\t\t<species2_gff>")
	print()
	print('Required: ')
	print('ortholog_file\ttab-delimited file with pairs of orthologous genes\n\t\t\t\tif labels not specified, first line must have species name')
	print('species1_gff\tGFF file for genes in species 1; gene IDs\n\t\t\t\tshould match those in ortholog file')
	print('species2_gff\tGFF file for genes in species 2; gene IDs\n\t\t\t\tshould match those in ortholog file')
	print()
	print('Optional: ')
	print('-l=label1,label2\t\tspecies labels to use in use in output GFF of other species;\n\t\t\t\t\t\t[default None--get label from first line of ortholog file]')
	print('-c=col1,col2\t\t\t0-indexed column number which has the species gene id;\n\t\t\t\t\t\t[Default 0,1]')
	print('-o1=outfile_species1\toutput file name for species 1 [default appends "ortholog.gff" to\n\t\t\t\t\t\tinput]')
	print('-o1=outfile_species1\toutput file name for species 2 [default appends "ortholog.gff"\n\t\t\t\t\t\tto input]')

	
if __name__ == "__main__":
	if len(sys.argv) < 4 :
		printHelp()
	else:
		parseInputs( sys.argv[1:] )
