#!/usr/bin/env python

import sys
import csv
import re
import argparse
from pysam import VariantFile


def load_ntc_dict(ntc_tsv, ntc_vcf):
	""" make a dictionary of all variants in the NTC """
	ntc_dict = {}
	with open(ntc_tsv) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for line in reader:
			var = line['Variant']
			ntc_dict[var] = {
				'ntc_vaf': line['Frequency'],
				'ntc_depth': line['Depth'],
				'ntc_alt_reads': add_from_vcf(var, ntc_vcf)['alt_reads'],
			}
	return ntc_dict


def split_genomic_coords(var):
	""" split genomic coords into its parts """
	chr = var.split(':')[0]
	pos = "".join(re.findall('\d+', var.split(':')[1]))
	ref, alt = re.sub(r'[0-9]', '', var.split(':')[1]).split('>')
	return chr, pos, ref, alt


def refomat_data(var):
	""" reformat SomAmp TSV into dictionary ready for SVD format TSV """
	chr, pos, ref, alt = split_genomic_coords(var['Variant'])
	reformatted_data = {
		'gene': var['Gene'],
		'chr': 'chr' + chr,
		'pos': pos,
		'ref': ref,
		'alt': alt,
		'vaf': var['Frequency'],
		'depth': var['Depth'],
		'hgvs_p': var['HGVSp'],
		'hgvs_c': var['HGVSc'],
		'consequence': var['Consequence'],
		'exon': var['EXON'].replace('|', '/'),
	}
	return reformatted_data


def check_ntc(var, ntc_dict):
	""" check if a variant is in the NTC and return values for final TSV """
	if var['Variant'] in ntc_dict.keys():
		variant = var['Variant']
		out_dict = {
			'in_ntc': True,
			'ntc_vaf': ntc_dict[variant]['ntc_vaf'],
			'ntc_depth': ntc_dict[variant]['ntc_depth'],
			'ntc_alt_reads': ntc_dict[variant]['ntc_alt_reads']
		}
	else:
		out_dict = {
			'in_ntc': False,
			'ntc_vaf': '',
			'ntc_depth': '',
			'ntc_alt_reads': ''
		}
	return out_dict


def add_from_vcf(var, vcf):
	""" add alt reads and gnomad pop max from the VCF for a specific variant """
	# split genomic coords
	chr, pos, ref, alt = split_genomic_coords(var)

	# find variant in VCF and pull out variables
	for v in vcf.fetch(chr, int(pos)-1, int(pos)):
		if v.ref == ref and v.alts[0] == alt:
			gnomad = v.info['gnomad_popmax_af']
			alt_reads = v.samples[0]['AD'][1]

	# return as dict
	out_dict = {
		'alt_reads': alt_reads,
		'gnomad_popmax_AF': gnomad,
	}
	return out_dict


def get_args():
	""" uses argparse to take arguments from command line """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--sample_tsv', action='store',
		help ='Filepath to sample *VariantReport.tsv file. REQUIRED.'
	)
	parser.add_argument(
		'--sample_vcf', action='store',
		help='Filepath to sample VCF. REQUIRED.'
	)
	parser.add_argument(
		'--ntc_tsv', action='store',
		help='Filepath to NTC *VariantReport.tsv file. REQUIRED.'
	)
	parser.add_argument(
		'--ntc_vcf', action='store',
		help='Filepath to NTC VCF. REQUIRED.'
	)
	parser.add_argument(
		'--output', action='store',
		help='Path to output file. REQUIRED.'
	)

	return parser.parse_args()


if __name__ == '__main__':
	# get input args
	args = get_args()
	sample_tsv = args.sample_tsv
	sample_vcf = VariantFile(args.sample_vcf)
	ntc_tsv = args.ntc_tsv
	ntc_vcf = VariantFile(args.ntc_vcf)

	out_list = []
	output_location = args.output

	# load in NTC variants
	ntc_dict = load_ntc_dict(ntc_tsv, ntc_vcf)

	# loop through sample variants TSV file
	with open(sample_tsv) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for var in reader:
			# only keep preferred transcripts
			if var['Preferred'] != 'True':
				continue

			# reformat data into SVD style
			reformatted_data = refomat_data(var)

			# check if in NTC and add results to reformatted dict
			temp_dict = check_ntc(var, ntc_dict)
			reformatted_data.update(temp_dict)

			# add gnomad and alt reads from VCF
			temp_dict = add_from_vcf(var['Variant'], sample_vcf)
			reformatted_data.update(temp_dict)

			out_list.append(reformatted_data)

	# output
	headers = ['gene', 'chr', 'pos', 'ref', 'alt', 'vaf', 'depth', 'hgvs_p', 'hgvs_c', 'consequence', 'exon', 'alt_reads', 'in_ntc', 'ntc_vaf', 'ntc_depth', 'ntc_alt_reads', 'gnomad_popmax_AF']
	with open(output_location, 'w') as f:
		csv_writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
		csv_writer.writeheader()
		csv_writer.writerows(out_list)
