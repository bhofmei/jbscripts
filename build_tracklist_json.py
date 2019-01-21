import sys, os
from io import open

# Usage: python build_tracklist_json.py [-q] <track_info_file>

COLOR_AR = ['blue', 'red', 'orange', 'yellow', 'purple', 'green', 'black', 'gray', 'cyan', 'magenta', 'darkred', 'darkblue' ]
COMP_COLOR_AR = ['red', 'blue', 'darkviolet', 'rebeccapurple', 'limegreen', 'fuchsia', 'deepskyblue', 'skyblue', 'magenta', 'cyan', 'darkturquoise', 'mediumvioletred' ]

'''
info needed for
# header is version
All: label, key, category, height
ChIP: modification type, path/file name, ggf run, desciption
Methylation: context, path/file name, description, method
RNA-seq: path bigwig, path bam, metadata (ecotype, ggf run, source, description)'
smRNA-seq: (path bigwig), path bam, medata, (ecotype, ggf run, source, description)
genes: orthologs, genome version, source
rnas: genome version, source
tes: genome version, source
dna: genome version, source
bed: label, key, category, chip_type, meta
atac-seq: bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
vcf: vcf(bigwig), description, source, metadata
gc-content
motif density: context

Track File: tab or comma separated with headers:
trackType label key category height chip_type/orthologs chip/meth/rna-seq_bigwig rna-seq_bam genome_version/description source_label/ggf_run source_link metadata_key_value
'''

def processInputs( trackInfoStr, isQuiet ):
	baseDir = os.path.dirname( trackInfoStr )
	outFileStr = os.path.join( baseDir, 'trackList.json' )
	outFile = open( outFileStr, 'w' )
	print( 'Writing to trackList.json' )
	outFile.write( getStarting() )
	# read file
	outStr, version, includeAr = readInfoFile( trackInfoStr, isQuiet )
	outFile.write( outStr )
	outFile.write( getEnding( version, includeAr ) )
	outFile.close()
	print( 'Done' )

def getStarting():
	outStr = u'{\n'
	outStr += tab(1) + u'"tracks" : [\n'
	return outStr

def getEnding( version, includeAr ):
	outStr = u'\n'+ tab(1) + u'],\n'
	# the names part
	outStr += tab(1) + u'"names" : {\n'
	outStr += tab(2) + u'"url" : "names/",\n'
	outStr += tab(2) + u'"type" : "Hash"\n'
	outStr += tab(1) + u'},\n'
	# check includeAr
	if len(includeAr) == 1:
		outStr += tab(1) + u'"include" : "'+ includeAr[0]+'",\n'
	elif len(includeAr) > 1:
		outStr += tab(1) + u'"include" : [' + ','.join( ['"{:s}"'.format(x) for x in includeAr ] ) + '],\n'
	outStr += tab(1) + u'"formatVersion" : '+version+'\n}\n'
	return outStr

def tab( n ):
	return u' '*(3*n)

def readInfoFile( trackInfoStr, isQuiet ):

	#trackTypeDict = {'dna':['dna'], 'genes':['genes','gene'], 'rnate': ['rnas', 'te', 'tes', 'transposons', 'repeats'], 'anno': ['anno'], 'chip': ['chip', 'chipseq', 'chip-str', 'chipseq-str'], 'reads': ['reads', 'read'], 'rnaseq': ['rnaseq', 'rnaseqpe', 'rnaseq-pe', 'rnaseq-pe-str', 'rnaseq-str'], 'smrna': ['smrna','smrnaseq', 'smrna-str', 'smrnaseq-str'], 'methyl': ['methyl', 'methylwig'], 'atac': ['atac', 'atacseq', 'atac-str', 'atacseq-str'], 'peaks': ['peak', 'peaks'], 'vcf':['vcf']}
	trackTypeDict = {'dna':['dna'], 'genes':['genes','gene'], 
	'rnate': ['rnas', 'te', 'tes', 'transposons', 'repeats'], 'anno': ['anno'], 
	'chip': ['chip', 'chipseq', 'chip-str', 'chipseq-str'], 'reads': ['reads', 'read'], 
	'rnaseq': ['rnaseq', 'rnaseqpe', 'rnaseq-pe', 'rnaseq-str', 'rnaseqpe-str', 'rnaseq-pe-str'],
	'smrna': ['smrna','smrnaseq','smrna-str','smrnaseq-str'],
	'methyl': ['methyl', 'methyl3','methylwig'],
	'atac': ['atac', 'atacseq', 'atac-str', 'atacseq-str'],
	'dap': ['dap', 'dapseq', 'dap-str', 'dapseq-str'],
	'starr': ['starr', 'starrseq', 'starr-str', 'starrseq-str'],
	'peaks': ['peak', 'peaks'], 
	'vcf':['vcf'], 'gc': ['gccont', 'gcdens'], 
	'motifdens': ['nucdens','motifdens']}
	trackTypeList = []
	for x in trackTypeDict.keys():
		trackTypeList += trackTypeDict[x]
	
	trackFile = open( trackInfoStr, 'rt' )
	outStr = ''
	version = '1'
	includeAr = []
	for line in trackFile:
		if line.startswith( '#' ):
			if line.startswith( '#Version:' ):
				lineAr = line.rstrip().split(',')
				version = lineAr[0].rstrip()[9:]
			elif line.startswith( '#Includes:' ):
				# includes need to be ";" separated
				cleanLine = line.replace('#Includes:','').rstrip().lstrip()
				if '\t' in cleanLine:
					cleanLineAr = cleanLine.split('\t')
				else:
					cleanLineAr = cleanLine.split(',')
				includeAr = cleanLineAr[0].split(';')
			continue
		line = line.rstrip()
		if '\t' in line:
			lineAr = line.split( '\t' )
		else:
			lineAr = line.split(',')
		if len(lineAr) < 2:
			continue
		# (0) trackType (1) label (2) key (3) category (4) height
		# (5) chip_type/orthologs/reads_type/color
		# (6) chip/meth/rna-seq_bigwig (7) reads/rna-seq_bam
		# (8) genome_version/description
		# (9) source_label/ggf_run (10) source_link
		# (11) mapping_rate (12) percent_remaining
		# (13) metadata_key_value
		trackType = lineAr[0].lower()
		
		if trackType not in trackTypeList:
			print( 'WARNING: {:s} is not a correct track type. Skipping...'.format( lineAr[0] ) )
			continue
		# commas for previous track
		elif outStr != '':
			outStr += ',\n'
		# handle type
		if not isQuiet:
			print( '-'+lineAr[1]+' ('+trackType+')' )
		if trackType in trackTypeDict['dna']:
			# label, key, category, genome_version, source_label, source_link
			info = lineAr[1:4] + lineAr[8:11]
			outStr += generateDNAtext( info )
			
		elif trackType in trackTypeDict['genes']:
			# label, key, category, track_height, orthologs, genome_version, source_label, source_link
			info = lineAr[1:6] + lineAr[8:11]
			outStr += generateGeneText( info )
			
		elif trackType in trackTypeDict['rnate']:
			#label, key, category, track_height, color, genome_version, source_label, source_link, trackType
			info = lineAr[1:6] + lineAr[8:11] + [lineAr[0]]
			outStr += generateRnaTeText( info )
		
		elif trackType in trackTypeDict['anno']:
			#label, key, category, track_height, chip_type, genome_version, source_label, source_link, meta
			info = lineAr[1:6] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateAnnoText( info )
			
		elif trackType in trackTypeDict['chip']:
			# label, key, category, track_height, chip_type, bigwig, description, ggf_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7] + lineAr[8:14]
			# determine stranded
			info += [ '-str' in trackType ]
			outStr += generateChipText( info )
		
		elif trackType in trackTypeDict['atac']:
			# label, key, category, track_height, color, bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7]+ lineAr[8:14]
			info += [ '-str' in trackType ]
			outStr += generateAtacText( info )
		
		elif trackType in trackTypeDict['dap']:
			# label, key, category, track_height, color, bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7]+ lineAr[8:14]
			info += [ '-str' in trackType ]
			outStr += generateDapText( info )
		
		elif trackType in trackTypeDict['starr']:
			# label, key, category, track_height, color, bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7]+ lineAr[8:14]
			info += [ '-str' in trackType ]
			outStr += generateStarrText( info )
			
		elif trackType in trackTypeDict['reads']:
			# label, key, category, track_height, type, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:6] + lineAr[7:14]
			outStr += generateReadsText( info )
		
		elif trackType in trackTypeDict['rnaseq']:
			# label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta, isPE, isStranded
			info = lineAr[1:14]
			# check PE
			info +=  ['pe' in trackType ]
			# check strand
			info += [ '-str' in trackType ]
			outStr += generateRnaSeqText( info )
			
		elif trackType in trackTypeDict['smrna']:
			# label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:14]
			# check strand
			info += [ '-str' in trackType ]
			outStr += generateSmRnaSeqText( info )
		
		elif trackType in trackTypeDict['methyl']:	
			if trackType == 'methyl' or trackType == 'methyl3':	# methylation verion 2 and 3
				# label, key, category, track_height, context, bigwig, description, gff_run/source, source label, meta, isV3
				info = lineAr[1:7] + lineAr[8:11] + [ lineAr[13] ]
				info += [ '3' in trackType ]
				outStr += generateMethylationTextv2( info )
			elif trackType == 'methylwig':	# methylation version 1
				# label, key, category, track_height, (chip_type,) bigwig, description, ggf_run/source, source_link, meta
				info = lineAr[1:5] + lineAr[6:7] + lineAr[8:11] + [ lineAr[13] ]
				outStr += generateMethylationTextv1( info )
		
		elif trackType in trackTypeDict['peaks']:
			# label, key, category, track_height, chip_type, meta
			info = lineAr[1:6] + [ lineAr[13] ]
			outStr += generatePeakBed( info )

		elif trackType in trackTypeDict['vcf']:
			# label, key, category, trackHeight, bigwig, description, source, source_link, meta
			info = lineAr[1:5] + [ lineAr[6] ] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateVCFText( info )
		
		elif trackType in trackTypeDict['gc']:
			# label, key, category, trackHeight, color, description, meta
			info = lineAr[1:6] + [ lineAr[8] ] + [ lineAr[13] ]
			info += [ 'dens' in trackType ]
			outStr += generateGCContent( info )	
	
		elif trackType in trackTypeDict['motifdens']:
			# label, key, category, trackHeight, contexts (can include colors)
			info = lineAr[1:6]
			outStr += generateMotifDensContent( info )
	# end for line
	trackFile.close()
	return outStr, version, includeAr

def generateGenomeSource( genomeVersion, sourceLabel, sourceLink ):
	outStr = ',\n'
	outStr += tab(3) + '"metadata" : {\n'
	if genomeVersion == "" and sourceLabel == "" and sourceLink == "":
		return "\n"
	if genomeVersion != "":
		outStr += tab(4) + '"Genome Version" : "{:s}"'.format( genomeVersion )
		if sourceLabel == "":
			outStr += '\n'
		else:
			outStr += ',\n'
	if sourceLabel != "":
		outStr += tab(4) + '"Source" : "{:s}"\n'.format( sourceLabel )
	outStr += tab(3) + '}'
	if sourceLabel != "" and sourceLink != "":
		outStr += ',\n' + tab(3) + '"fmtMetaValue_Source" : "function(source) { return \' <a href='+sourceLink+'>'+sourceLabel+'</a>\';}"\n'
	else:
		outStr += '\n'
	return outStr

def generateMeta( description, sLabel, sLink, mapRate, perRemain, meta ):
	"""
		generate metadata information including source formating if necessary
		always ends with comma
	"""
	if description == "" and sLabel == "" and sLink == "" and meta == "" and mapRate == "" and perRemain == "":
		return '\n'
	before = ""
	outStr = tab(3) + '"metadata" : {'
	# ggf run or source
	if sLabel.startswith( 'run' ):
		outStr += '\n' +tab(4) + '"GGF Run" : "{:s}"'.format( sLabel )
	elif sLink != "":
		outStr +='\n' + tab(4) + '"Source" : "{:s}"'.format( sLabel )
		before = tab(3) + '"fmtMetaValue_Source" : "function(source) { return \' <a href='+sLink+'>'+sLabel+'</a>\';}",\n'
	elif sLabel != "":
		outStr += '\n' + tab(4) + '"Source" : "{:s}"'.format( sLabel )
	# description
	if description != "":
		if sLabel != '':
			outStr += ','
		outStr += '\n' + tab(4) + '"Description" : "{:s}"'.format( description)
	# mapping rate/download raw
	if mapRate != "":
		if sLabel != "" or description != "":
			outStr += ','
		try:
			f = float(mapRate)
			outStr += '\n' + tab(4) + '"Read Mapping Rate" : "{:s}%"'.format( mapRate )
			if perRemain != "":
				if sLabel != "" or description != "" or mapRate != "":
					outStr += ','
				outStr += '\n' + tab(4) + '"Percent Reads Remaining" : "{:s}%"'.format( perRemain )
		except ValueError:
			# we have download label and/or link
			outStr +='\n' + tab(4) + '"Download" : "{:s}"'.format( mapRate )
			if perRemain != "":
				before += tab(3) + '"fmtMetaValue_Download" : "function(source) { return \' <a href='+perRemain+'>'+mapRate+'</a>\';}",\n'
	# other metadata
	if meta != "":
		if description!= '' or sLabel != '' or mapRate != '' or perRemain != '':
			outStr += ','
		outStr += parseMetaKeys( meta )
	outStr += '\n' + tab(3) + '},\n'
	return before + outStr
	
def generateDNAtext( infoAr ):
	'''
		infoAr = [label, key, category, trackHeight, genome_version, source_label, source_link]
		infoAr[i] == '' if info not available
	'''
	label, key, category, gVersion, sLabel, sLink = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"chunkSize" : 20000,\n' + tab(3) + '"storeClass" : "JBrowse/Store/Sequence/StaticChunked",\n'
	outStr += tab(3) + '"urlTemplate" : "seq/{refseq_dirpath}/{refseq}-",\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"type" : "SequenceTrack"'
	# genome version and source
	outStr += generateGenomeSource( gVersion, sLabel, sLink )
	outStr += tab(2) + '}'
	return outStr
	
def generateGeneText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, orthologs, genome_version, source_label, source_link]
	'''
	label, key, category, tHeight, orthologs, gVersion, sLabel, sLink = infoAr
	color = getColors( label )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histogram and style
	if color != None:
		outStr += tab(3) + '"histograms" : {\n' + tab(4) + ' "color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature-genes",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	else:
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"description" : "description,note"\n' + tab(3) +'},\n'

	# check orthologs
	if orthologs != "":
		orthoAr = orthologs.split(';')
		for ortho in orthoAr:
			outStr += generateOrtholog( ortho )
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '500' if tHeight == '' else tHeight )
	outStr += tab(3) + '"maxFeatureScreenDensity" : 0.1,\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"type" : "JBrowse/View/Track/CanvasFeatures"'
	outStr += generateGenomeSource( gVersion, sLabel, sLink )
	outStr += tab(2) + '}'
	return outStr
	
def generateOrtholog( species ):
	formatSpecies = getOrthologFormat( species )
	if formatSpecies == False:
		return ''
	
	outStr = tab(3) + "\"fmtDetailValue_"+formatSpecies+"-ortholog\" : \"function(name, feature) {if(feature.get('type')=='gene'){ return ' <a href=http://epigenome.genetics.uga.edu/JBrowse/?data="+species+"&loc='+name+'>'+name+'</a>'; } }\",\n"
	return outStr

def generateRnaTeText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, genome_version, source_label, source_link, type]
	'''
	label, key, category, tHeight, sColor, gVersion, sLabel, sLink, tType = infoAr
	# no color specified -> use default
	if sColor == '':
		color = getColors( label )
		histColor = color
	# two colors specified -> feature color; hist color
	elif ';' in sColor:
		colorAr = sColor.split(';')
		color = colorAr[0]
		histColor = colorAr[1]
	# single color specified -> color
	else:
		color = sColor
		histColor = sColor
	
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histogram and style
	if color != None:
		outStr += tab(3) + '"histograms" : {\n' + tab(4) + ' "color" : "{:s}"\n'.format( histColor ) + tab(3) + '},\n'
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature-' + tType + '",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	else:
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature-' + tType + '"\n' + tab(3) +'},\n'
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '400' if tHeight == '' else tHeight )
	outStr += tab(3) + '"maxFeatureScreenDensity" : 0.1,\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"type" : "CanvasFeatures"'
	outStr += generateGenomeSource( gVersion, sLabel, sLink )
	outStr += tab(2) + '}'
	return outStr
	
def generateAnnoText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, color, genome_version, source_label, source_link, meta]
	'''
	label, key, category, tHeight, iColor, gVersion, sLabel, sLink , meta= infoAr
	
	color = getColors( (label if iColor == '' else iColor) )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# style
	outStr += tab(3) + '"histograms" : {\n' + tab(4) + ' "color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	# check for style meta info
	styleTxt, meta = parseMetaStyle( meta )
	outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature",\n' + tab(4) + '"color" : "{:s}"'.format( color )
	if styleTxt != '':
		outStr += ',\n' + styleTxt + '\n'
	outStr += tab(3) + '},\n'
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '400' if tHeight == '' else tHeight )
	outStr += tab(3) + '"maxFeatureScreenDensity" : 0.1,\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += generateMeta( gVersion, sLabel, sLink, '', '', meta )
	outStr += tab(3) + '"type" : "CanvasFeatures"\n'
	outStr += tab(2) + '}'
	return outStr
	
def generateChipText( infoAr ):
	'''
		infoAr = [label, key, category, tsHeight, chip_type, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, tsHeight, chipType, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta, stranded = infoAr
	tHeight, scaleType = getHeightScale( tsHeight )
	
	if stranded:
		colorAr =  getStrandedColors( chipType )
	else:
		color = getColors( chipType )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	if stranded and colorAr != None:
		outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( colorAr[0] )
		outStr += tab(4) + '"neg_color" : "{:s}",\n'.format( colorAr[1] )
	elif not stranded:
		outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : {:s}\n'.format( '50' if tHeight == '' else tHeight )
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	if stranded:
		outStr += tab(3) + '"storeClass" : "StrandedPlotPlugin/Store/SeqFeature/StrandedBigWig",\n'
		outStr += tab(3) + '"type" : "StrandedPlotPlugin/View/Track/Wiggle/StrandedXYPlot",\n'
	else:
		if scaleType != 'local' and scaleType != '':
			outStr += tab(3) + '"autoscale": "{:s}",\n'.format( scaleType )
		else:
			outStr += tab(3) + '"autoscale": "local",\n'
		outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
		outStr += tab(3) + '"type" : "JBrowse/View/Track/Wiggle/XYPlot",\n'
		outStr += tab(3) + '"min_score" : 0,\n' 
	
	outStr += tab(3) + '"urlTemplate" : "raw/chip/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateAtacText( infoAr ):
	# atac-seq is the same as chip except urlTemplate
	if infoAr[4] == "":
		infoAr[4] = 'atac'
	outStr = generateChipText( infoAr )
	outStr = outStr.replace( 'raw/chip/', 'raw/atac/' )
	return outStr
	
def generateDapText( infoAr ):
	# dap-seq is the same as chip except urlTemplate
	if infoAr[4] == "":
		infoAr[4] = 'dap'
	outStr = generateChipText( infoAr )
	outStr = outStr.replace( 'raw/chip/', 'raw/dap/' )
	return outStr

def generateStarrText( infoAr ):
	# starr-seq is the same as chip except urlTemplate
	if infoAr[4] == "":
		infoAr[4] = 'starr'
	outStr = generateChipText( infoAr )
	outStr = outStr.replace( 'raw/chip/', 'raw/starr/' )
	return outStr

def generateMethylationTextv1( infoAr ):
	'''
		infoAr = [label, key, category, track_height, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, tHeight, bigWig, desc, sLabel, sLink, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"height" : {:s}\n'.format( '70' if tHeight == '' else tHeight )
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"min_score" : -1,\n'
	outStr += tab(3) + '"max_score" : 1,\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/methyl/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, '','', meta )
	outStr += tab(3) + '"type" : "MethylationPlugin/View/Track/Wiggle/MethylXYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateMethylationTextv2( infoAr ):
	'''
		infoAr = [label, key, category, track_height, context, bigwig, description, ggf_run/source, source_link, meta, isV3]
	'''
	label, key, category, tHeight, mContext, bigWig, desc, sLabel, sLink, meta, isV3 = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"height" : {:s}\n'.format( '70' if tHeight == '' else tHeight )
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"min_score" : -1,\n'
	outStr += tab(3) + '"max_score" : 1,\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"storeClass" : "MethylationPlugin/Store/SeqFeature/MethylBigWig",\n'
	if mContext != '':
		mAr = mContext.split(';')
		mmAr = [ '"{:s}"'.format( x.lower() ) for x in mAr ]
		outStr += tab(3) + '"context" : [ ' + ', '.join(mmAr) + ' ],\n'
	outStr += tab(3) + '"urlTemplate" : "raw/methyl/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, '','', meta )
	outStr += tab(3) + '"type" : "MethylationPlugin/View/Track/Wiggle/MethylPlot",\n'
	# if v3, add methylation option
	if isV3:
		outStr += tab(3) + '"methylatedOption" : true,\n'
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateRnaSeqText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, meta, isPE, isStranded]
	'''
	label, key, category, tHeight, height, bigWig, bam, desc, sLabel, sLink, mapRate, perRemain, meta, isPE, isStranded = infoAr
	if height == "":
		maxheight = "7000"
		minheight = '0'
	else:
		heightAr = height.split('#')
		if len(heightAr)==1:
			minheight = '0'
			maxheight = heightAr[0]
		else:
			minheight = heightAr[0]
			maxheight = heightAr[1]
	

	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histograms
	outStr += tab(3) + '"histograms" : {\n'
	outStr += tab(4) + '"color" : "#c4c4c4",\n'
	if isStranded:
		outStr += tab(4) + '"storeClass" : "StrandedPlotPlugin/Store/SeqFeature/StrandedBigWig",\n'
	else:
		outStr += tab(4) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
		outStr += tab(4) + '"min" : {:s},\n'.format( minheight )
		outStr += tab(4) + '"max" : {:s},\n'.format( maxheight )
	outStr += tab(4) + '"urlTemplate" : "raw/rna/{:s}",\n'.format( bigWig )
	outStr += tab(4) + '"description" : "coverage depth",\n'
	outStr += tab(4) + '"height" : 100\n'
	outStr += tab(3) + '},\n'
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BAM",\n'
	#outStr += tab(3) + '"maxFeatureScreenDensity" : 2,\n'
	if isPE:
		outStr += tab(3) + '"useReverseTemplateOption": true,\n'
		outStr += tab(3) + '  "useReverseTemplate": true,\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '400' if tHeight == '' else tHeight )
	outStr += tab(3) + '"urlTemplate" : "raw/rna/{:s}",\n'.format( bam )
	outStr += tab(3) + '"category" : "{:s}",\n'.format(category)
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Alignments2",\n'
	outStr += tab(3) + '"glyph" : "JBrowse/View/FeatureGlyph/Alignment",\n'
	outStr += tab(3) + '"maxFeatureSizeForUnderlyingRefSeq" : 250000\n'
	outStr += tab(2) + '}'
	return outStr

def generateSmRnaSeqText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, meta, isStranded]
	'''
	label, key, category, tHeight, height, bigWig, bam, desc, sLabel, sLink, mapRate, perRemain, meta, isStranded = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histograms
	if bigWig != '':
		outStr += tab(3) + '"histograms" : {\n'
		outStr += tab(4) + '"color" : "#d1d1d1",\n'
		if height != "" and not isStranded:
			outStr += tab(4) + '"max" : {:s},\n'.format( height )
		if isStranded:
			outStr += tab(4) + '"storeClass" : "StrandedPlotPlugin/Store/SeqFeature/StrandedBigWig",\n'
		else:
			outStr += tab(4) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
		outStr += tab(4) + '"urlTemplate" : "raw/smrna/{:s}",\n'.format( bigWig )
		outStr += tab(4) + '"description" : "coverage depth",\n'
		outStr += tab(4) + '"height" : 100\n'
		outStr += tab(3) + '},\n'
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BAM",\n'
	outStr += tab(3) + '"maxFeatureScreenDensity" : 1.5,\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '300' if tHeight == '' else tHeight )
	outStr += tab(3) + '"urlTemplate" : "raw/smrna/{:s}",\n'.format( bam )
	outStr += tab(3) + '"category" : "{:s}",\n'.format(category)
	outStr += tab(3) + '"type" : "SmallRNAPlugin/View/Track/smAlignments",\n'
	outStr += tab(3) + '"maxFeatureSizeForUnderlyingRefSeq" : 500000\n'
	outStr += tab(2) + '}'
	return outStr

def generateReadsText( infoAr ):
	'''
	label, key, category, track_type, type, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
	'''
	label, key, category, tHeight, folder, bam, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BAM",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format('500' if tHeight == '' else tHeight)
	outStr += tab(3) + '"urlTemplate" : "raw/{:s}/{:s}",\n'.format( folder, bam )
	outStr += tab(3) + '"category" : "{:s}",\n'.format(category)
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Alignments2",\n'
	outStr += tab(3) + '"glyph" : "JBrowse/View/FeatureGlyph/Alignment",\n'
	outStr += tab(3) + '"maxFeatureSizeForUnderlyingRefSeq" : 250000\n'
	outStr += tab(2) + '}'
	return outStr

def generateGCContent( infoAr ):
	'''
		infoAr = [label, key, category, track_height, color, description, meta, isDens]
	'''
	label, key, category, tHeight, color, desc, meta, isDens = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	if tHeight != '' or color != '':
		outStr += tab(3) + '"style" : {\n'
		if tHeight != '':
			outStr += tab(4) + '"height" : {:s},\n'.format( tHeight )
		colorAr = getStrandedColors( color )
		if colorAr != None:
			outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( colorAr[0] )
			outStr += tab(4) + '"neg_color" : "{:s}",\n'.format( colorAr[1] )
	# end if tHeight
		outStr += tab(3) + '},\n'
	outStr += generateMeta( desc, '', '', '','', meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/SequenceChunks",\n'
	outStr += tab(3) + '"urlTemplate" : "seq/{refseq_dirpath}/{refseq}-",\n'
	outStr += tab(3) + '"bicolor_pivot": 0.5,\n'
	outStr += tab(3) + '"type": "GCContent/View/Track/GCContent{:s}",\n'.format( '' if isDens else 'XY' )
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateMotifDensContent( infoAr ):
	'''
		infoAr = [label, key, category, track_height, color]
	'''
	label, key, category, tHeight, contexts = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	if tHeight != '':
		outStr += tab(3) + '"style" : {\n'
		outStr += tab(4) + '"height" : {:s},\n'.format( tHeight )
		outStr += tab(3) + '},\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/SequenceChunks",\n'
	outStr += tab(3) + '"urlTemplate" : "seq/{refseq_dirpath}/{refseq}-",\n'
	ctxStr, clrStr = parseContexts( contexts )
	outStr += tab(3) + '"motifs" : {:s},\n'.format( ctxStr )
	if clrStr != '':
		outStr += tab(3) + '"colors" : {:s},\n'.format( clrStr )
	outStr += tab(3) + '"type": "MotifDensityPlugin/View/Track/MotifDensity",\n'
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr
	
def parseContexts( contextsStr ):
	if contextsStr == '':
		return '["CG"]', ''
	outColor = '{'
	outContexts = '['
	contextsAr = contextsStr.split(';')
	for context in contextsAr:
		if ':' in context: # color specified
			ctx, clr = context.split( ':' )
			outColor += '"{:s}" : "{:s}",'.format( ctx, clr )
			outContexts += '"{:s}",'.format( ctx ) 
		else:
			outContexts += '"{:s}",'.format( context.upper() ) 
	outColor = outColor[:-1]
	outContexts = outContexts[:-1]
	if outColor != '':
		outColor += '}'
	outContexts += ']'
	return outContexts, outColor
	
def parseMetaStyle( metaStr ):
	if metaStr == "":
		return "", ""
	styleStr = ""
	outStr = ""
	isFirst = True
	pairsAr = metaStr.split(';')
	for pair in pairsAr:
		key, value = pair.split(':')
		if key.startswith( 'style.' ):
			key = key.replace('style.', '' )
			if isFirst:
				isFirst  = False
			else:
				styleStr += ',\n'
			styleStr += tab(4) + '"'+key+'" : "'+value+'"'
		else:
			outStr += pair + ';'
	# end for
	if outStr != '':
		outStr = outStr[:-1]
	return styleStr, outStr
	
def parseMetaKeys( metaStr ):
	if metaStr == "":
		return '\n'
	outStr = ''
	isFirst = True
	pairsAr = metaStr.split(';')
	for pair in pairsAr:
		key, value = pair.split(':')
		if isFirst:
			isFirst = False
			outStr += '\n'
		else:
			outStr += ',\n'
		outStr += tab(4) + '"'+key+'" : "'+value+'"'
	return outStr

def generatePeakBed( infoAr ):
	'''
		infoAr = [label, key, category, track_height, chipType, meta]
	'''
	label, key, category, tHeight, chipType, meta = infoAr
	color = getColors( chipType )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"className" : "feature6",\n'
	outStr += tab(4) + '"showLabels" : false,\n'
	outStr += tab(4) + '"arrowheadClass" : null,\n'
	outStr += tab(4) + '"featureCss" : "background-color:{:s};border-color:black",\n'.format( color )
	outStr += tab(4) + '"histCss" : "background-color:{:s};border-color:black"\n'.format( color )
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += generateMeta( '', '', '', '','',meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "HTMLFeatures",\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"type" : "JBrowse/View/Track/HTMLFeatures",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format('500' if tHeight == '' else tHeight)
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

	
def generateVCFText( infoAr ):
	label, key, category, tHeight, vcf, desc, sLabel, sLink, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"urlTemplate" : "raw/vcf/{:s}",\n'.format( vcf )
	outStr += generateMeta( desc, sLabel, sLink, '', '', meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/VCFTabix",\n'
	outStr += tab(3) + '"type" : "JBrowse/View/Track/HTMLVariants",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format('200' if tHeight == '' else tHeight)
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def getHeightScale( heightStr ):
	height = ''
	scale = ''
	scaleTypes = ['local', 'global', 'clipped_global']
	
	ar = heightStr.split(';')
	for x in ar:
		if x in scaleTypes:
			scale = x
		else:
			try:
				n = int( x.replace('px','') )
				height = str(n)
			except ValueError:
				height =''
	# end for
	return height, scale
		

def getStrandedColors( typeStr ):
	# no color specifed, return None
	if typeStr == '':
		return None
	# otherwise split by ';'
	typeAr = typeStr.split(';')
	# two colors listed
	if len(typeAr) == 2:
		tmpAr = [ getColors(x) for x in typeAr ]
		if tmpAr[0] == 'black' and tmpAr[1] == 'black':
			typeAr = tmpAr # get the complimentary color
		else:
			return tmpAr
	# one color listed -> get the color then complement it
	color1 = getColors( typeAr[0] )
	color2 = getComplimentaryColor( color1 )
	return [color1, color2]

def getColors( typeStr ):
	if typeStr.startswith( 'gene' ):
		typeStr = 'gene'
	elif typeStr.startswith('transcript'):
		typeStr = 'gene'
	elif typeStr.startswith( 'transposon' ):
		typeStr = 'tes'
	elif typeStr.startswith( 'rnas' ):
		typeStr = 'rnas' 
	typeDict = { 'genes':'#daa520', 'gene':'#daa520', 
		'rnas':'#18A071', 'rna':'#18A071',
		'tes':'#77158D', 'repeats':'#77158D', 'transposons':'#77158D',
		'h2ax-elements':'#daa520', 'regulatory':'#da2055',
		'mcg':'#A36085', 'mchg':'#0072B2', 'mchh':'#CF8F00',
		'cg':'#A36085', 'chg':'#0072B2','chh':'#CF8F00',
		'h2az':'#EE7600', 'h3':'#937E6B', 
		'h3k4m1':'#701941', 'h3k4m2':'#711B71', 'h3k4m3':'#862BB3',
		'h3k4me1':'#701941', 'h3k4me2':'#711B71', 'h3k4me3':'#862BB3',
		'h3k9ac':'#51d657', 'h3k9m2':'#228B22', 'h3k9m3':'#165C16',
		'h3k9me2':'#228B22', 'h3k9me3':'#165C16',
		'h3k23ac':'#ED735D', 'h3k27ac':'#11B1DC',
		'h3k27m3':'#225EA8', 'h3k27me3':'#225EA8',
		'h3k36m1':'#701D19', 'h3k36m2':'#A61B03', 'h3k36m3':'#D42727',
		'h3k36me1':'#701D19', 'h3k36me2':'#A61B03', 'h3k36me3':'#D42727', 
		'h3k56ac':'#DB3B7E', 'input':'#E3AC22', 'polii':'#006f6f',
		'h3t32':'#008b8b', 'h2bs112':'#6d32c3','h3t32':'#32a2a2', 
		'sdg7':'#2e8b57', 'basej':'#228b22', 
		'h3t32g':'#00688b', 'methyl':'#a1a1a1',
		'atac': '#768C9C', 'dap':'#343A40', 'starr':'#685D79'
		 }
	outStr = typeDict.get( typeStr.lower() )
	if outStr == None:
		if typeStr.lower() in COLOR_AR:
			return typeStr
		elif typeStr.startswith( '#' ):
			return typeStr
		return 'black'
	return outStr

def getComplimentaryColor( colorStr ):
	if colorStr in COLOR_AR:
		fInd = COLOR_AR.index(colorStr)
		return COMP_COLOR_AR[fInd]
	else:
		colorStr=colorStr.replace('#','')
		r = hex( 255-int(colorStr[0:2],base=16) )[2:]
		if len(r) == 1:
			r = '0'+r
		g = hex( 255-int(colorStr[2:4],base=16) )[2:]
		if len(g) == 1:
			g = '0'+g
		b = hex( 255-int(colorStr[4:6],base=16) )[2:]
		if len(b) == 1:
			b = '0'+b
		return '#' + r + g + b
		

def getOrthologFormat( orthoStr ):
	strDict = { 'poplar': 'Ptrichocarpa', 'arabidopsis':'Athaliana',
		'eutrema':'Esalsugineum', 'maize':'Zmays', 'bdistachyon': 'Bdistachyon',
		'maize_v4': 'B73', 'maize_ph207':'Ph207'}
	outStr = strDict.get( orthoStr )
	if outStr == None:
		print( 'WARNING: ortholog name {:s} not recognized'.format( orthoStr ) )
		return False
	return outStr
	
		
def parseInputs( argv ):
	if argv[0] == '-q':
		trackInfoStr = argv[1]
		isQuiet = True
	else:
		trackInfoStr = argv[0]
		isQuiet = False
	processInputs( trackInfoStr, isQuiet )


if __name__ == "__main__":
	if len(sys.argv) < 2 :
		print ("Usage: python2 build_tracklist_json.py [-q] <track_info_file>\n converts a specifically formatted CSV file to a trackList.json file\n use -q for quiet option")
	else:
		parseInputs( sys.argv[1:] )
