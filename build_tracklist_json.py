import sys, math, glob, multiprocessing, subprocess, os, bisect, random
from io import open

# Usage: python build_tracklist_json.py <track_info_file>

'''
info needed for
# header is version
All: label, key, category
ChIP: modification type, path/file name, ggf run, desciption
Methylation: context, path/file name, description, method
RNA-seq: path bigwig, path bam, metadata (ecotype, ggf run, source, description)
genes: orthologs, genome version, source
rnas: genome version, source
tes: genome version, source
dna: genome version, source
bed: label, key, category, chip_type, meta
atac-seq: bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta

Track File: tab separated with headers:
trackType label key category chip_type/meth_con/orthologs chip/meth/rna-seq_bigwig rna-seq_bam genome_version/description source_label/ggf_run source_link metadata_key_value
'''

def processInputs( trackInfoStr ):
	baseDir = os.path.dirname( trackInfoStr )
	outFileStr = os.path.join( baseDir, 'trackList.json' )
	outFile = open( outFileStr, 'w' )
	print( 'Writing to trackList.json' )
	outFile.write( getStarting() )
	# read file
	outStr, version = readInfoFile( trackInfoStr )
	outFile.write( outStr )
	outFile.write( getEnding( version ) )
	outFile.close()
	print( 'Done' )
	

def getStarting():
	outStr = u'{\n'
	outStr += tab(1) + u'"tracks" : [\n'
	return outStr

def getEnding( version ):
	outStr = u'\n'+ tab(1) + u'],\n'
	# the names part
	outStr += tab(1) + u'"names" : {\n'
	outStr += tab(2) + u'"url" : "names/",\n'
	outStr += tab(2) + u'"type" : "Hash"\n'
	outStr += tab(1) + u'},\n'
	outStr += tab(1) + u'"formatVersion" : '+version+'\n}\n'
	return outStr

def tab( n ):
	return u' '*(3*n)

def readInfoFile( trackInfoStr ):
	
	trackFile = open( trackInfoStr, 'rt' )
	outStr = ''
	version = '1'
	for line in trackFile:
		if line.startswith( '#' ):
			if line.startswith( '#Version:' ):
				lineAr = line.rstrip().split(',')
				version = lineAr[0].rstrip()[9:]
			continue
		lineAr = line.rstrip().split(',')
		if len(lineAr) < 2:
			continue
		# (0) trackType (1) label (2) key (3) category
		# (4) chip_type/orthologs
		# (5) chip/meth/rna-seq_bigwig (6) rna-seq_bam
		# (7) genome_version/description
		# (8) source_label/ggf_run (9) source_link
		# (10) mapping_rate (11) percent_remaining
		# (12) metadata_key_value
		# acceptable track types: dna, genes, rnas, te, repeats, chip, rnaseq, methyl, peaks
		trackType = lineAr[0].lower()
		if trackType not in ['dna', 'genes', 'rnas', 'te', 'repeats', 'chip', 'rnaseq', 'methyl', 'peaks', 'atac','atacseq','rnastrand']:
			print( 'WARNING: {:s} is not a correct track type. Skipping...'.format( lineAr[0] ) )
			continue
		# commas for previous track
		elif outStr != '':
			outStr += ',\n'
		print( '-'+lineAr[1] )
		# handle type
		if trackType == 'dna':
			# label, key, category, genome_version, source_label, source_link
			info = lineAr[1:4] + lineAr[7:10]
			outStr += generateDNAtext( info )
			
		elif trackType == 'genes':
			# label, key, category, orthologs, genome_version, source_label, source_link
			info = lineAr[1:5] + lineAr[7:10]
			outStr += generateGeneText( info )
			
		elif trackType in ['rnas', 'te', 'repeats']:
			#label, key, category, genome_version, source_label, source_link
			info = lineAr[1:4] + lineAr[7:10]
			outStr += generateRnaTeText( info )
			
		elif trackType == 'chip':
			# label, key, category, chip_type, bigwig, description, ggf_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:6] + lineAr[7:13]
			outStr += generateChipText( info )
		elif trackType == 'rnaseq':
			# label, key, category, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:13]
			outStr += generateRnaSeqText( info )
		elif trackType == 'methyl':
			# label, key, category, (chip_type,) bigwig, description, ggf_run/source, source_link, meta
			''' info = lineAr[1:6] + lineAr[7:10] + [ lineAr[12] ] '''
			info = lineAr[1:4] + lineAr[5:6] + lineAr[7:10] + [ lineAr[12] ]
			outStr += generateMethylationText( info )
		elif trackType == 'peaks':
			# label, key, category, chip_type, meta
			info = lineAr[1:5] + [ lineAr[12] ]
			outStr += generatePeakBed( info )
		elif trackType in ['atac','atacseq']:
			# label, key, category, bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:4] + [ lineAr[5] ] + lineAr[7:13]
			outStr += generateAtacText( info )
		elif trackType == 'rnastrand':
			# label, key, category, chip_type, bigwig, description, ggf_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:6] + lineAr[7:13]
			outStr += generateRNAStrandText( info )
	# end for line
	trackFile.close()
	return outStr, version

def generateDNAtext( infoAr ):
	'''
		infoAr = [label, key, category, genome_version, source_label, source_link]
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
	
def generateGeneText( infoAr ):
	'''
		infoAr = [label, key, category, orthologs, genome_version, source_label, source_link]
	'''
	label, key, category, orthologs, gVersion, sLabel, sLink = infoAr
	color = getColors( label )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histogram and style
	if color != None:
		outStr += tab(3) + '"histograms" : {\n' + tab(4) + ' "color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	else:
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature"\n' + tab(3) +'},\n'
	
	# check orthologs
	if infoAr[3] != "":
		orthoAr = orthologs.split(';')
		for ortho in orthoAr:
			outStr += generateOrtholog( ortho )
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : 500,\n'
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
		infoAr = [label, key, category, genome_version, source_label, source_link]
	'''
	label, key, category, gVersion, sLabel, sLink = infoAr
	color = getColors( label )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histogram and style
	if color != None:
		outStr += tab(3) + '"histograms" : {\n' + tab(4) + ' "color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	else:
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature"\n' + tab(3) +'},\n'
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : 400,\n'
	outStr += tab(3) + '"maxFeatureScreenDensity" : 0.1,\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"type" : "CanvasFeatures"'
	outStr += generateGenomeSource( gVersion, sLabel, sLink )
	outStr += tab(2) + '}'
	return outStr
	
def generateChipText( infoAr ):
	'''
		infoAr = [label, key, category, chip_type, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, chipType, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	color = getColors( chipType)
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : 50\n'
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"autoscale": "clipped_global",\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/chip/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Wiggle/XYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"min_score" : 0\n' + tab(2) + '}'
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
	# mapping rate
	if mapRate != "":
		if sLabel != "" or description != "":
			outStr += ','
		outStr += '\n' + tab(4) + '"Read Mapping Rate" : "{:s}%"'.format( mapRate )
	if perRemain != "":
		if sLabel != "" or description != "" or mapRate != "":
			outStr += ','
		outStr += '\n' + tab(4) + '"Percent Reads Remaining" : "{:s}%"'.format( perRemain )
	# other metadata
	if meta != "":
		if description!= '' or sLabel != '' or mapRate != '' or perRemain != '':
			outStr += ','
		outStr += parseMetaKeys( meta )
	outStr += '\n' + tab(3) + '},\n'
	return before + outStr

def generateMethylationText( infoAr ):
	'''
		infoAr = [label, key, category, chip_type, bigwig, description, ggf_run/source, source_link, meta]
	'''
	'''label, key, category, context, bigWig, desc, sLabel, sLink, meta = infoAr
	color = getColors( context)
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"neg_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : 50\n'
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/methyl/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, '','', meta )
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Wiggle/XYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"min_score" : -1,\n'
	outStr += tab(3) + '"max_score" : 1\n'
	outStr += tab(2) + '}' '''
	label, key, category, bigWig, desc, sLabel, sLink, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"height" : 70\n'
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"min_score" : -1,\n'
	outStr += tab(3) + '"max_score" : 1,\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/methyl/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, '','', meta )
	outStr += tab(3) + '"type" : "MethylationPlugin/View/Track/Wiggle/MethylXYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateRnaSeqText( infoAr ):
	'''
		infoAr = [label, key, category, height, bigwig, bam, description, source/ggf_run, source_link, meta]
	'''
	label, key, category, height, bigWig, bam, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	if height == "":
		maxheight = "1000"
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
	outStr += tab(4) + '"color" : "gray",\n'
	outStr += tab(4) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(4) + '"min" : {:s},\n'.format( minheight )
	outStr += tab(4) + '"max" : {:s},\n'.format( maxheight )
	outStr += tab(4) + '"urlTemplate" : "raw/rna/{:s}",\n'.format( bigWig )
	outStr += tab(4) + '"description" : "coverage depth",\n'
	outStr += tab(4) + '"height" : 100\n'
	outStr += tab(3) + '},\n'
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BAM",\n'
	outStr += tab(3) + '"maxFeatureScreenDensity" : 2,\n'
	outStr += tab(3) + '"maxHeight" : 500,\n'
	outStr += tab(3) + '"urlTemplate" : "raw/rna/{:s}",\n'.format( bam )
	outStr += tab(3) + '"category" : "{:s}",\n'.format(category)
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Alignments2",\n'
	outStr += tab(3) + '"glyph" : "JBrowse/View/FeatureGlyph/Alignment",\n'
	outStr += tab(3) + '"maxFeatureSizeForUnderlyingRefSeq" : 250000\n'
	outStr += tab(2) + '}'
	return outStr

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
		infoAr = [label, key, category, chipType, meta]
	'''
	label, key, category, chipType, meta = infoAr
	color = getColors( chipType )
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
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
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateAtacText( infoAr ):
	label, key, category, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( 'gray' )
	outStr += tab(4) + '"height" : 50\n'
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"autoscale": "clipped_global",\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/atac/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Wiggle/XYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += tab(3) + '"min_score" : 0\n' + tab(2) + '}'
	return outStr

def generateRNAStrandText( infoAr ):
	'''
		infoAr = [label, key, category, color, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, color, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"neg_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : 50\n'
	outStr += tab(3) + '},\n'
	outStr += tab(3) + '"variance_band" : false,\n'
	outStr += tab(3) + '"autoscale": "clipped_global",\n'
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/BigWig",\n'
	outStr += tab(3) + '"urlTemplate" : "raw/rna/{:s}",\n'.format( bigWig )
	outStr += generateMeta( desc, sLabel, sLink, mapRate, perRemain, meta )
	outStr += tab(3) + '"type" : "JBrowse/View/Track/Wiggle/XYPlot",\n'
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr +=  tab(2) + '}'
	return outStr

def getColors( typeStr ):
	colors = ['blue', 'red', 'orange', 'yellow', 'purple', 'green', 'black', 'gray', 'cyan', 'magenta' ]
	typeDict = { 'h2az':'#ee7600', 'h3':'#8b7765', 'h3k4m3':'#9a32cd',
		'h3k9m2':'#228b22', 'h3k56ac':'#ee1289', 'input':'#708090',
		'h3k36m3':'#ee2c2c', 'h3k27m3':'#3a5fcd', 'h3t32':'#008b8b',
		'h2bs112':'#6d32c3', 'genes':'#daa520', 'gene':'#daa520', 
		'rnas':'#18A071', 'rna':'#18A071', 'tes':'#77158D', 'repeats':'#77158D',
		'transposons':'#77158D', 'mcg':'#b03060', 'mchg':'#2e8b57',
		'mchh':'#1e90ff', 'cg':'#b03060', 'chg':'#2e8b57',
		'chh':'#1e90ff','h3k27m3':'#617ed7','h3t32':'#32a2a2',
		'h3k36m1':'#AE2020','h3k36m2':'#D42727', 'h3k4m1':'#6A228D',
		'h3k4m2':'#872CB3','sdg7':'#2e8b57','basej':'#228b22', 
		'h3t32g' :'#00688b', 'methyl':'#a1a1a1' }
	outStr = typeDict.get( typeStr.lower() )
	if outStr == None:
		if typeStr in colors:
			return typeStr
		return 'black'
	return outStr

def getOrthologFormat( orthoStr ):
	strDict = { 'poplar': 'Ptrichocarpa', 'arabidopsis':'Athaliana',
		'eutrema':'Esalsugineum', 'maize':'Zmays' }
	outStr = strDict.get( orthoStr )
	if outStr == None:
		print( 'WARNING: ortholog name {:s} not recognized'.format( orthoStr ) )
		return False
	return outStr
	
		
def parseInputs( argv ):
	trackInfoStr = argv[0]
	processInputs( trackInfoStr )


if __name__ == "__main__":
	if len(sys.argv) != 2 :
		print (" Usage: python2 build_tracklist_json.py <track_info_file>")
	else:
		parseInputs( sys.argv[1:] )
