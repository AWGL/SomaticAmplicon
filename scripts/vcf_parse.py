#!/usr/bin/env python

"""
Convert annotated vcf from 1_SomaticAmplicon.sh into readable variant report file
"""
import re
import os
import csv
import logging
import argparse
import textwrap
from pysam import VariantFile
from pyvariantfilter.utils import parse_csq_field

# Main functions to create variant report file
#	"""
#	Logger deals with all status messages
#	"""

logger = logging.getLogger('vcf_parse.vcf')

def parse_transcripts(pref_trans):
	"""
	Get list of preferred transcripts
	"""

	#Empty variable
	pref_transcript = []

	#Make list of transcripts
	trans = open(pref_trans, 'r')

	for line in trans:

		#Skip header
		if line.startswith('Transcript'):
			next
		else:
			transcript = line.strip().split('\t')[0]
			pref_transcript.append(transcript)
	trans.close()

	return(pref_transcript)

	logger.info('Loaded preferred transcripts')

def parse_ntc(ntc_file):
	"""
	Create dictionary of all NTC variants
	"""

	# Empty ntc list
	ntc_list = []

	# Loop over the ntc file
	ntc = VariantFile(ntc_file)

	for line in ntc.fetch():

		chrom = line.chrom
		pos = line.pos
		ref = line.ref
		alt = line.alts[0]
		var = chrom+":"+str(pos)+ref+">"+alt
		vaf = line.info['AF']
		# Convert vaf tuple to string and remove parentheses
		vaf = str(vaf).replace('(','').replace(')','').replace(',','')
		dp = line.samples[0]['DP']

		#For alt read count, get last value from AD
		alt = line.samples[0]['AD'][1]

		ntc_var = [var, vaf, dp, alt]
		ntc_list.append(ntc_var)

	return(ntc_list)

	logger.info('created list of NTC variants')

def parse_sample(inp, outfile, ntc_list, pref_transcript):
	"""
	Loads in data from vcf
	Returns list of variants in correct format
	Only gives ntc_vaf and ntc_depth if in_ntc = True
	"""

	# Read filtered meta annotated vcf produced in 1_SomaticAmplicon.sh
	logger.info(
		'loading VCF file from {}'.format(os.path.abspath(inp)))

	# Empty list to hold variants
	variant_list = []

	# Read vcf file
	inp = VariantFile(inp)
	# CSQ fields from header to be parsed with pyvariantfilter
	csq_fields = str(inp.header.info['CSQ'].record)
	csq_fields = csq_fields.strip()
	index = csq_fields.index('Format') + 8
	csq_fields = csq_fields[index:len(csq_fields)-2].split('|')

	for record in inp.fetch():

		if 'CSQ' in record.info:

			# Load chrom, pos, ref & alt fields
			chrom = record.chrom
			pos = record.pos
			ref = record.ref
			alt = record.alts[0]

			# Load vaf and depth
			vaf = record.info['AF']
	                # Convert vaf tuple to string and remove parentheses
			vaf = str(vaf).replace('(','').replace(')','').replace(',','')	
			depth = record.info['DP']
			# Alt-reads - take second value from AD
			alt_reads = record.samples[0]['AD'][1]

			# Sample for output filepath

			# Get annotations and go through these
			csq = record.info['CSQ']
			annotations = parse_csq_field(csq, csq_fields)

			# Empty list of annotations
			transcripts = []

			# If any of the 'Feature' values match the preferred transcript
			for entry in annotations:

				# Loop over annotations and if preferred transcript add to transcript list
				if entry['Feature'] in pref_transcript:
					gene = entry['SYMBOL']

			# If transcripts for variant still empty, add MANE select
			if len(transcripts) == 0:

				for entry in annotations:

					if entry['MANE']:

						gene = entry['SYMBOL']
						hgvs_p = entry['HGVSp']
						hgvs_c = entry['HGVSc']
						consequence = entry['Consequence']
						exon = entry['EXON']

						transcript_info = (gene
								, hgvs_p
								, hgvs_c
								, consequence
								, exon)

						transcripts.append(transcript.info)

			# If transcripts still empty, use PICK annotation
			if len(transcripts) == 0:

				for entry in annotations:

					if 'PICK':

						gene = entry['SYMBOL']
						hgvs_p = entry['HGVSp']
						hgvs_c = entry['HGVSc']
						consequence = entry['Consequence']
						exon = entry['EXON']

						transcript_info = (gene
								, hgvs_p
								, hgvs_c
								, consequence
								, exon)

						transcripts.append(transcript_info)

			# If still empty leave empty
			if len(transcripts) == 0:

				gene = ""
				hgvs_p = ""
				hgvs_c = ""
				consequence = ""
				exon = ""

				transcript_info = (gene
						, hgvs_p
						, hgvs_c
						, consequence
						, exon)

				transcripts.append(transcript.info)

			#Check if NTC variant list empty
			if len(ntc_list) == 0:
				in_ntc = ""
				ntc_vaf = ""
				ntc_depth = ""
				ntc_alt_reads = ""
			else:
				#Check if variant in NTC
				for entry in ntc_list:
					variant_full = chrom+":"+str(pos)+ref+">"+alt
					if variant_full == entry[0]:
						in_ntc = "True"
						ntc_vaf = entry[1]
						ntc_depth = entry[2]
						ntc_alt_reads = entry[3]
					else:
						in_ntc = "False"
						ntc_vaf = ""
						ntc_depth = ""
						ntc_alt_reads = ""

			#Go through all annotations we're keeping for this variant and add them
			for keep in transcripts:
				var = f'{keep[0]}\t{chrom}\t{pos}\t{ref}\t{alt}\t{vaf}\t{depth}\t{keep[1]}\t{keep[2]}\t{keep[3]}\t{keep[4]}\t{alt_reads}\t{in_ntc}\t{ntc_vaf}\t{ntc_depth}\t{ntc_alt_reads}'
				variant_list.append(var)


	return(variant_list)


## Input arguments
def get_args():
	"""
	Uses argparse to take arguments from command line
	"""

	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		description=textwrap.dedent(
		'''
		Summary:
		Takes VCF file created by 1_SomaticAmplicon.sh and parses
		the variants to produce a tab delimited variant report
                '''
	))

	# Required: VCF input file
	parser.add_argument('-i', '--inp', action='store',
		help ='Filepath to input VCF file. REQUIRED.'
	)

	# OPTIONAL: Output folder, defaults to current directory if empty
	parser.add_argument(
		'-o', '--outfile', action='store', 
		help=textwrap.dedent(
		'''
		Filepath to folder where output reports will be saved. 
		If missing, defaults to current directory.
		\n'''
	))

	# Required: Preferred transcript file
	parser.add_argument(
		'-t', '--pref_trans', action='store',
		help='Filepath to preferred transcrips'
	)

	# Required: Psth to ntc vcf file
	parser.add_argument(
		'-n', '--ntc_vcf', action='store',
		help='Filepath to NTC VCF'
	)

	return parser.parse_args()

# Function to run main functions
def main(args):

	# Setup logger
	logger = logging.getLogger('vcf_parse')
	logger.setLevel(logging.DEBUG)

	handler = logging.StreamHandler()
	handler.setLevel(logging.DEBUG)

	formatter = logging.Formatter(
		'%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s'
	)

	handler.setFormatter(formatter)

	logger.addHandler(handler)
	logger.info('running vcf_parse.py...')

	#Load arguments, make vcf report object and load data

	# Parse preferred transcripts
	if args.pref_trans:
		pref_transc = parse_transcripts(args.pref_trans)
	else:
		logger.info('no preffered transcript file provided - please provide')

	# Parse NTC
	if args.ntc_vcf:
		ntc_list = parse_ntc(args.ntc_vcf)
	else:
		logger.info('no ntc vcf provided - please provide')

	# Parse sample
	variant_list = parse_sample(args.inp, args.outfile, ntc_list, pref_transc)

	# Finish
	logger.info('steps completed\n{}'.format('---'*30))

	# Export data
	# Header lines
	logger.info('creating output file')

	outfile_headers = ['gene', 'chr', 'pos', 'ref', 'alt', 'vaf', 'depth', 'hgvs_p', 'hgvs_c', 'consequence', 'exon', 'alt_reads', 'in_ntc', 'ntc_vaf', 'ntc_depth', 'ntc_alt_reads']
	outfile_headers_string = '\t'.join(outfile_headers)

	#Info
	with open(args.outfile,'w',newline='') as file:
		file.write(outfile_headers_string + '\n')

		for line in variant_list:
			file.write(str(line) + '\n')

	# Final file
	logger.info('vcf_parse.py completed\n{}'.format('---'*30))

# Call functions
if __name__ == '__main__':
	args = get_args()
	main(args)
