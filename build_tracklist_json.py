import sys, os
from io import open

# Usage: python build_tracklist_json.py [-q] <track_info_file>

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

Track File: tab or comma separated with headers:
trackType label key category height chip_type/meth_con/orthologs chip/meth/rna-seq_bigwig rna-seq_bam genome_version/description source_label/ggf_run source_link metadata_key_value
'''

def processInputs( trackInfoStr, isQuiet ):
	baseDir = os.path.dirname( trackInfoStr )
	outFileStr = os.path.join( baseDir, 'trackList.json' )
	outFile = open( outFileStr, 'w' )
	print( 'Writing to trackList.json' )
	outFile.write( getStarting() )
	# read file
	outStr, version = readInfoFile( trackInfoStr, isQuiet )
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

def readInfoFile( trackInfoStr, isQuiet ):
	
	trackFile = open( trackInfoStr, 'rt' )
	outStr = ''
	version = '1'
	for line in trackFile:
		if line.startswith( '#' ):
			if line.startswith( '#Version:' ):
				lineAr = line.rstrip().split(',')
				version = lineAr[0].rstrip()[9:]
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
		if trackType not in ['dna', 'genes', 'rnas', 'te', 'tes', 'transposons', 'repeats', 'chip', 'rnaseq', 'rnaseqpe','smrna','smrnaseq','methyl','methylwig', 'peaks', 'atac', 'atacseq', 'reads', 'read', 'rnastrand', 'vcf','anno']:
			print( 'WARNING: {:s} is not a correct track type. Skipping...'.format( lineAr[0] ) )
			continue
		# commas for previous track
		elif outStr != '':
			outStr += ',\n'
		# handle type
		if not isQuiet:
			print( '-'+lineAr[1]+' ('+trackType+')' )
		if trackType == 'dna':
			# label, key, category, genome_version, source_label, source_link
			info = lineAr[1:4] + lineAr[8:11]
			outStr += generateDNAtext( info )
			
		elif trackType == 'genes':
			# label, key, category, track_height, orthologs, genome_version, source_label, source_link
			info = lineAr[1:6] + lineAr[8:11]
			outStr += generateGeneText( info )
			
		elif trackType in ['rnas', 'te', 'tes', 'transposons', 'repeats']:
			#label, key, category, track_height, genome_version, source_label, source_link
			info = lineAr[1:5] + lineAr[8:11]
			outStr += generateRnaTeText( info )
		
		elif trackType == 'anno':
			#label, key, category, track_height, chip_type, genome_version, source_label, source_link, meta
			info = lineAr[1:6] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateAnnoText( info )
			
		elif trackType == 'chip':
			# label, key, category, track_height, chip_type, bigwig, description, ggf_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7] + lineAr[8:14]
			outStr += generateChipText( info )
		elif trackType in ['reads','read']:
			# label, key, category, track_height, type, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:6] + lineAr[7:14]
			outStr += generateReadsText( info )
		elif trackType == 'rnaseq':
			# label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta, isPE
			info = lineAr[1:14] + [ False ]
			outStr += generateRnaSeqText( info )
		elif trackType == 'rnaseqpe':
			# label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:14] + [ True ]
			outStr += generateRnaSeqText( info )
		elif trackType in ['smrna','smrnaseq']:
			# label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:14]
			outStr += generateSmRnaSeqText( info )
		elif trackType == 'methyl':	# methylation verion 2
			# label, key, category, track_height, context, bigwig, description, gff_run/source, source label, meta
			info = lineAr[1:7] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateMethylationTextv2( info )
		elif trackType == 'methylwig':	# methylation version 1
			# label, key, category, track_height, (chip_type,) bigwig, description, ggf_run/source, source_link, meta
			info = lineAr[1:5] + lineAr[6:7] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateMethylationTextv1( info )
		elif trackType == 'peaks':
			# label, key, category, track_height, chip_type, meta
			info = lineAr[1:6] + [ lineAr[13] ]
			outStr += generatePeakBed( info )
		elif trackType in ['atac','atacseq']:
			# label, key, category, track_height, bigwig, description, gff_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:5] + [ lineAr[6] ] + lineAr[8:14]
			outStr += generateAtacText( info )
		elif trackType == 'rnastrand':
			# label, key, category, trackHeight, chip_type, bigwig, description, ggf_run/source, source_link, mapping_rate, percent_remaining, meta
			info = lineAr[1:7] + lineAr[8:14]
			outStr += generateRNAStrandText( info )
		elif trackType == 'vcf':
			# label, key, category, trackHeight, bigwig, description, source, source_link, meta
			info = lineAr[1:5] + [ lineAr[6] ] + lineAr[8:11] + [ lineAr[13] ]
			outStr += generateVCFText( info )
	# end for line
	trackFile.close()
	return outStr, version

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
		outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
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
		infoAr = [label, key, category, track_height, genome_version, source_label, source_link]
	'''
	label, key, category, tHeight, gVersion, sLabel, sLink = infoAr
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
	outStr += tab(3) + '"style" : {\n' + tab(4) + '"className" : "feature",\n' + tab(4) + '"color" : "{:s}"\n'.format( color ) + tab(3) + '},\n'
	# basics
	outStr += tab(3) + '"storeClass" : "JBrowse/Store/SeqFeature/NCList",\n'
	outStr += tab(3) + '"trackType" : "CanvasFeatures",\n'
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format( '400' if tHeight == '' else tHeight )
	outStr += tab(3) + '"maxFeatureScreenDensity" : 0.1,\n'
	outStr += tab(3) + '"urlTemplate" : "tracks/'+label+'/{refseq}/trackData.json",\n'
	outStr += tab(3) + '"compress" : 0,\n'
	outStr += tab(3) + '"category" : "{:s}",\n'.format( category )
	outStr += generateMeta( gVersion, sLabel, sLink, '', '', meta )
	outStr += tab(3) + '"type" : "CanvasFeatures"'
	outStr += tab(2) + '}'
	return outStr
	
def generateChipText( infoAr ):
	'''
		infoAr = [label, key, category, chip_type, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, tHeight, chipType, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	color = getColors( chipType)
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : {:s}\n'.format( '50' if tHeight == '' else tHeight )
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
		infoAr = [label, key, category, track_height, context, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, tHeight, mContext, bigWig, desc, sLabel, sLink, meta = infoAr
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
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateRnaSeqText( infoAr ):
	'''
		infoAr = [label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, meta, isPE]
	'''
	label, key, category, tHeight, height, bigWig, bam, desc, sLabel, sLink, mapRate, perRemain, meta, isPE = infoAr
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
	#outStr += tab(3) + '"maxFeatureScreenDensity" : 2,\n'
	if isPE:
		outStr += tab(3) + '"useReverseTemplateOption": true,\n'
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
		infoAr = [label, key, category, track_height, height, bigwig, bam, description, source/ggf_run, source_link, meta]
	'''
	label, key, category, tHeight, height, bigWig, bam, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	# histograms
	if bigWig != '':
		outStr += tab(3) + '"histograms" : {\n'
		outStr += tab(4) + '"color" : "#d1d1d1",\n'
		if height != "":
			outStr += tab(4) + '"max" : {:s},\n'.format( height )
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
	outStr += tab(3) + '"maxHeight" : {:s},\n'.format('500' if tHeight == '' else tHeight)
	outStr += tab(3) + '"category" : "{:s}"\n'.format( category )
	outStr += tab(2) + '}'
	return outStr

def generateAtacText( infoAr ):
	label, key, category, tHeight, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( 'gray24' )
	outStr += tab(4) + '"height" : {:s}\n'.format('50' if tHeight == '' else tHeight)
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
		infoAr = [label, key, category, track_height, color, bigwig, description, ggf_run/source, source_link, meta]
	'''
	label, key, category, tHeight, color, bigWig, desc, sLabel, sLink, mapRate, perRemain, meta = infoAr
	outStr = tab(2) + '{\n'
	outStr += tab(3) + '"key" : "{:s}",\n'.format( key )
	outStr += tab(3) + '"label" : "{:s}",\n'.format( label )
	outStr += tab(3) + '"style" : {\n'
	outStr += tab(4) + '"clip_marker_color" : "black",\n'
	outStr += tab(4) + '"pos_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"neg_color" : "{:s}",\n'.format( color )
	outStr += tab(4) + '"height" : {:s}\n'.format('50' if tHeight == '' else tHeight)
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

def getColors( typeStr ):
	colors = ['blue', 'red', 'orange', 'yellow', 'purple', 'green', 'black', 'gray', 'cyan', 'magenta', 'darkred', 'darkblue' ]
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
		'h3t32g' :'#00688b', 'methyl':'#a1a1a1','h2ax-elements':'#daa520' }
	outStr = typeDict.get( typeStr.lower() )
	if outStr == None:
		if typeStr.lower() in colors:
			return typeStr
		elif typeStr.startswith( '#' ):
			return typeStr
		return 'black'
	return outStr

def getOrthologFormat( orthoStr ):
	strDict = { 'poplar': 'Ptrichocarpa', 'arabidopsis':'Athaliana',
		'eutrema':'Esalsugineum', 'maize':'Zmays', 'bdistachyon': 'Bdistachyon' }
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
