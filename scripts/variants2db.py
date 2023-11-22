import sys
import csv
import re
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

if __name__ == '__main__':
	# TODO make argparse
	sample_tsv = sys.argv[1]
	sample_vcf = VariantFile(sys.argv[2])
	ntc_tsv = sys.argv[3]
	ntc_vcf = VariantFile(sys.argv[4])

	ntc_dict = load_ntc_dict(ntc_tsv, ntc_vcf)
	out_list = []

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

	# TODO - format and output
	headers = ['gene', 'chr', 'pos', 'ref', 'alt', 'vaf', 'depth', 'hgvs_p', 'hgvs_c', 'consequence', 'exon', 'alt_reads', 'in_ntc', 'ntc_vaf', 'ntc_depth', 'ntc_alt_reads', 'gnomad_popmax_AF']

	print('\t'.join(headers))
	for v in out_list:
		l = []
		for field in headers:
			l.append(str(v[field]))
		print('\t'.join(l))
