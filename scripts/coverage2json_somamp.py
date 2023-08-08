#!/usr/bin/env python

import json
import os
import argparse
import decimal
from decimal import Decimal
import logging
import pandas as pd
import numpy as np

"""
- script to format the output of CoverageCalculatorPy (.totalcoverage, .coverage) into a JSON format for importing into the new database
- for use as part of the CRM/BRCA  pipeline
- file requirements:
	- for main gene regions:
		- NGHS-101X/hotspot_coverage/<referral>.groups
		- NGHS-101X/NTC_coverage_results/<NTC_id>_<referral_type>.totalCoverage
		- NGHS-101X/sample_coverage_results/<sample_id>_<referral_type>.totalCoverage
	- for 'hotspots':
		- NGHS-101X/hotspot_coverage/<referral>.bed
		- NGHS-101X/NTC_coverage_results/<NTC_id>_<referral_type>.coverage
		- NGHS-101X/sample_coverage_results_270/<sample_id>_<referral_type>.coverage
	- for gaps (hotspot only):
		- NGHS-101X/sample_coverage_results_500/<sample_id>_<referral_type>.gaps
- usage: python coverage2json.py -r <referral type> -g <hotspots_coverage folder path> -s <sample Coverage_results folder> -n <ntc Coverage_results folder> -o <output file name>
- output: {sample_id}_{referral}_coverage.json file
{ 
	"gene1"	:	{
		"average_depth"		:	<integer>,
		"percent_500"		:	<integer>,
		"average_ntc"		:	<integer>,
		"percent_ntc"		:	<integer>,
		"hotspot_regions"	:	[
			[<chr>, <start>, <end>, "gene1"<gene info>, <avdepth>, <%@500>, <ntc depth>, <%ntc>]
		],
		"gaps_500"	:	[
			[<chr>, <start>, <end>, "gene1"<gene info>, <cosmic>],
			[<chr>, <start>, <end>, "gene1"<gene info>, <cosmic>]
		]
	}
	"gene2"	:	{
		"average_depth"		:	<integer>,
		"percent_500"		:	<integer>,
		"average_ntc"		:	<integer>,
		"percent_ntc"		:	<integer>,
		"hotspot_regions"	:	[
			[<chr>, <start>, <end>, "gene2"<gene info>, <avdepth>, <%@500>, <ntc depth>, <%ntc>]
		],
		"gaps_500"	:	[
			[<chr>, <start>, <end>, "gene2"<gene info>, <cosmic>],
			[<chr>, <start>, <end>, "gene2"<gene info>, <cosmic>]
		]
	}

}
"""

#####################################
#####         Functions         #####
#####################################

# Adding logging info to detect errors and workflow progress
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def np_encoder(object):
	"""
	- function needed to allow json.dump to parse np values correctly
	"""
	if isinstance(object, np.generic):
		return object.item()


def parse_referral_type_files(referral_type, groups_folder):
	"""
	- Function to parse referral type bed and group files
	- input: referral type (eg: 'Thyroid'), pathway to hotspot_coverage folder containing the groups/bed files
	- output:
	1) list of distinct gene names from <referral>.groups file
	3) df of <referral>.bed if present, blank if not present for referral type
	4) True or False value for if <referral>.bed is present
	5) True or False value for if <referral>.bed is present
	"""

	## create zero counter
	counter = 0

	dir_list = os.listdir(groups_folder)

	for file in dir_list:
		filepath = os.path.join(groups_folder,file)

		## parse groups file
		if file == f'{referral_type}.groups':
			groups_df = pd.read_csv(filepath, sep = '\t', index_col = False)
			gene_df = groups_df.drop_duplicates()
			gene_list = gene_df['GENE'].values.tolist()

		## parse bed file
		elif file == f'{referral_type}.bed':
			counter += 1
			df = pd.read_csv(filepath, sep = '\t', names = ['CHR', 'START', 'END', 'INFO'], index_col = False)
			df.drop(columns = ['INFO'], inplace = True)

	## set variable as not all referral types have bed/groups
	if counter > 0:
		present = True
		logger.info('bed/groups file available for {referral_type}')
	else:
		present = False
		logger.info('No bed/groups file for {referral_type}, creating blank dataframe')

		## create blank dataframe to return when no <referral>.bed file present
		df = pd.DataFrame()

	return gene_list, df, present


def parse_NTC_data(NTC_coverage_folder, referral_type, present):
	"""
	- function to parse NTC data from .totalCoverage and .coverage files from CoverageCalculatorPy
	- input: path to NTC Coverage_results folder and referral type
	- output:
	1) df of .totalcoverage gene level (gene id, ntc avg depth)
	2) df of .coverage region level (chr, start, end, meta, ntc avg depth) or blank df
	"""

	## parse .totalCoverage and .coverage for NTC_<referral>
	cov_filepath = NTC_coverage_folder
	dir_list = os.listdir(cov_filepath)

	## rip sample_id from filename
	sample_id = dir_list[0].split('_')[4]

	for file in dir_list:
		filepath = os.path.join(cov_filepath, file)
		
		## parse .totalCoverage
		if (referral_type in file) and ("NTC-" in file) and (".totalCoverage" in file):
			ntc_gene_df = pd.read_csv(filepath, sep = '\t', index_col = False)

			logger.info('NTC coverage file identified')

			## drop final row in df as this is the total average across total panel
			ntc_gene_df.drop(ntc_gene_df.tail(1).index, inplace = True)

			## get gene name positioning in 'feature' column (can change dependent on referral type)
			feature_line = ntc_gene_df['FEATURE'][0].split('_')
			gene_pos = len(feature_line) - 2

			## create 'gene' list and add to DF as a column
			gene_list = []
		
			for line in ntc_gene_df['FEATURE']:
				splitline = line.split('_')
				## if-loop to get around final row in the DF not containing gene name
				if len(splitline) > gene_pos:
					gene_list.append(splitline[gene_pos])
				else:
					gene_list.append()

			ntc_gene_df.insert(0,'GENE', gene_list)

			## drop irrelevant perc_coverage and feature column
			ntc_gene_df.drop(columns = ['PERC_COVERAGE@500', 'FEATURE'], inplace = True)

			ntc_gene_df.rename(columns={'AVG_DEPTH':'NTC_AVG_DEPTH'}, inplace = True)

		## parse .coverage for hotspots
		elif (referral_type in file) and ("NTC-" in file) and (".coverage" in file):
			ntc_region_df = pd.read_csv(filepath, sep = '\t', index_col = False)
			ntc_region_df.drop(columns = ['PERC_COVERAGE@500'], inplace = True)
			ntc_region_df.rename(columns={'AVG_DEPTH':'NTC_AVG_DEPTH'}, inplace = True)

	ntc_gene_df['NTC_AVG_DEPTH'] = ntc_gene_df['NTC_AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	
	## if hotspots present, format avg depth, or create empty df if not present
	if present:
		ntc_region_df['NTC_AVG_DEPTH'] = ntc_region_df['NTC_AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	else:
		ntc_region_df = pd.DataFrame()

	return ntc_gene_df, ntc_region_df


def parse_sample_data(sample_coverage_folder, referral_type, sample_id, present):
	"""
	- function to parse sample data from .totalCoverage and .coverage files from CoverageCalculatorPy
	- input: path to sample Coverage_results folder and referral type
	- output:
	1) sample_id
	2) df of .totalcoverage gene level (gene id, avg depth, % cov 500)
	3) df of .coverag region level (chr, start, end, meta, avg depth, % cov 500) or blank df
	4) df of 500 .gaps (chr, start, end, meta, cosmic)
	"""

	## parse 500x
	dir_list = os.listdir(sample_coverage_folder)

	for file in dir_list:
		filepath = os.path.join(sample_coverage_folder, file)

		## parse .totalCoverage for gene level info
		if (sample_id in file) and (referral_type in file) and (".totalCoverage" in file):
			sample_500_gene_df = pd.read_csv(filepath, sep = '\t', index_col = False)


			## drop final row in df as this is the total average across total panel
			sample_500_gene_df.drop(sample_500_gene_df.tail(1).index, inplace = True)

			## get gene name positioning in 'feature' column (can change dependent on referral type)
			feature_line = sample_500_gene_df['FEATURE'][0].split('_')
			gene_pos = len(feature_line) - 2

			## create 'gene' list and add to DF as a column
			gene_list = []
			for line in sample_500_gene_df['FEATURE']:
				splitline = line.split('_')

				## if loop to get around final row in the DF not containing gene name
				if len(splitline) > gene_pos:
					gene_list.append(splitline[gene_pos])
				else:
					gene_list.append()

			sample_500_gene_df.insert(0,'GENE', gene_list)

			## drop irrelevant feature column
			sample_500_gene_df.drop(columns = ['FEATURE'], inplace = True)

		## parse .coverage
		elif (sample_id in file) and (referral_type in file) and (".coverage" in file):
			sample_500_region_df = pd.read_csv(filepath, sep = '\t', index_col = False)
		
	return sample_id, sample_500_gene_df, sample_500_region_df

def parse_cosmic_data(sample_id, referral_type, cosmic_file_loc):
	'''
	- function to take output from cosmic annotations and put into dataframe - this includes the gaps
	'''
	## list files
	dir_list = os.listdir(cosmic_file_loc)
	
	for file in dir_list:

		filepath = os.path.join(cosmic_file_loc, file)
				
		## parse cosmic file - at 500X
		if (sample_id in file) and (referral_type in file) and ("cosmic.csv" in file):

			cosmic_500_df = pd.read_csv(filepath, sep = ',', index_col = False)

	return cosmic_500_df
	

def create_output_dict(gene_list, main_gene_df, region_df, present, cosmic_500_df):
	'''
	- function to collate all dataframes into one nested dictionary to export to JSON format
	- input: gene list, gene level df, region df, 500 gaps df, cosmic 500 df
	- output: dictionary formatted as JSON specified at top of script
	'''

	output_dict = {}
	for gene in gene_list:
		output_dict[gene] = {}

	for key in output_dict.keys():

		output_dict[key]['average_depth'] = main_gene_df.at[key, 'AVG_DEPTH']
		output_dict[key]['percent_500'] = main_gene_df.at[key, 'PERC_COVERAGE@500']
		output_dict[key]['average_ntc'] = main_gene_df.at[key, 'NTC_AVG_DEPTH']
		output_dict[key]['percent_ntc'] = main_gene_df.at[key, 'PERC_NTC_DEPTH']

		## bed
		if present:
			filtered_region_list = []
			region_list = region_df.to_dict(orient='records')

			for item in region_list:
				#Renaming the dictionary keys
				item['chr'] = item.pop('CHR')
				item['pos_start'] = item.pop('START')
				item['pos_end'] = item.pop('END')
				item['hgvs_c'] = item.pop('META')
				item['average_coverage'] = item.pop('AVG_DEPTH')
				item['percent_500'] = item.pop('PERC_COVERAGE@500')
				item['ntc_coverage'] = item.pop('NTC_AVG_DEPTH')
				item['percent_ntc'] = item.pop('PERC_NTC_DEPTH')

				# check if gene is the gene in the 4th column then add
				if key == item['hgvs_c'].split('(')[0]:
					filtered_region_list.append(item)

			output_dict[key]['hotspot_regions'] = filtered_region_list
												
		else:
			output_dict[key]['hotspot_regions'] = []

		## gaps and cosmic - converting into a dictionary instead of a list
		#cosmic_500_df.columns = ['chr', 'pos_start', 'pos_end', 'info', 'gene', 'count_cosmic','percent_cosmic']
		if not cosmic_500_df.empty:

			cosmic_500_list = cosmic_500_df.to_dict(orient='records')
			print(cosmic_500_list)
			cosmic_500_final_list = []
		
			for item in cosmic_500_list:
				print(item)
				#Renaming the dictionary keys
				item['chr'] = item.pop('Chr')
				item['pos_start'] = item.pop('Start')
				item['pos_end'] = item.pop('End')
				item['hgvs_c'] = item.pop('Info')
				item['gene'] = item.pop('Gene')
				item['counts_cosmic'] = item.pop('Counts')
				item['percent_cosmic'] = item.pop('Percentage')
				
				if key == item['gene']:
					cosmic_500_final_list.append(item)

			#For referrals with hotspot that don't have cosmic annotation, gene is not in fourth column so get that info
				elif isinstance(item['gene'],float):
					if np.isnan(item['gene']):
						gene = item['hgvs_c'].split("(")[0]
						item['gene'] = gene
						if key == gene:
							cosmic_500_final_list.append(item)

			output_dict[key]['gaps_500'] = cosmic_500_final_list 

		else:
			output_dict[key]['gaps_500'] = []

	return output_dict



#########################################
#####           Programme           #####
#########################################


if __name__ == '__main__':


	## args
	parser = argparse.ArgumentParser()
	parser.add_argument('--referral','-r', help = 'Referral type. eg: Thyroid')
	parser.add_argument('--groups_folder','-g', help = 'pathway to hotspots_coverage folder')
	parser.add_argument('--ntc_coverage','-n', help = 'pathway to NTC Coverage_results folder')
	parser.add_argument('--sample_id','-id', help = 'prefix for sample ID')
	parser.add_argument('--sample_coverage','-s', help = 'pathway to sample Coverage_results folder')
	parser.add_argument('--cosmic_file', '-c', help = 'pathway to sample Cosmic annotations folder')
	parser.add_argument('--outfile', '-o', help = 'pathway/file of output')
	args = parser.parse_args()


	## set decimal settings
	decimal.getcontext().rounding = decimal.ROUND_DOWN


	### parse referral_type group/bed files
	gene_list, df, present = parse_referral_type_files(args.referral, args.groups_folder)
	
	### parse NTC sample 500x for average depth per gene (.totalCoverage) and per region (.coverage)
	ntc_gene_df, ntc_region_df = parse_NTC_data(args.ntc_coverage, args.referral, present)

	### parse sample 135x and 500x files
	sample_id, sample_500_gene_df, sample_500_region_df = parse_sample_data(args.sample_coverage, args.referral, args.sample_id, present)

	### parse cosmic annotation output
	cosmic_500_df = parse_cosmic_data(args.sample_id, args.referral, args.cosmic_file)

	### create json pieces
	## join NTC data to df if present, else create blank df
	if present:

		logger.info('applying NTC data to json')
		main_region_df = pd.merge(sample_500_region_df, ntc_region_df, how = 'outer', on = ['CHR', 'START', 'END', 'META'])
		#Need to catch instances where AVG DEPTH = 0, and therefore PERC NTC DEPTH would calculate to be inf
		#First work out perc ntc depth
		PERC_NTC_DEPTH = (main_region_df['NTC_AVG_DEPTH'] / main_region_df['AVG_DEPTH']) * 100
		#If average depth is greater than 0, use the value above
		main_region_df.loc[main_region_df['AVG_DEPTH'] > 0, 'PERC_NTC_DEPTH'] = PERC_NTC_DEPTH
		#Otherwise use 0
		main_region_df.loc[main_region_df['AVG_DEPTH'] == 0, 'PERC_NTC_DEPTH'] = 0

		#Then reformat
		main_region_df['PERC_NTC_DEPTH'] = main_region_df['PERC_NTC_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))

	else:
		main_region_df = pd.DataFrame()
		logger.info('no NTC data to append')

	## join ntc information to gene level df
	main_gene_df = pd.merge(sample_500_gene_df, ntc_gene_df, how = 'outer', on = ['GENE'])


	## change index of gene df to be the gene name for iloc later
	main_gene_df.set_index('GENE', inplace = True)


	## create and format percent NTC column - need to catch instances where average depth = 0 as above
	main_gene_df.loc[main_gene_df['AVG_DEPTH'] > 0, 'PERC_NTC_DEPTH'] = (main_gene_df['NTC_AVG_DEPTH'] / main_gene_df['AVG_DEPTH']) * 100
	main_gene_df.loc[main_gene_df['AVG_DEPTH'] == 0, 'PERC_NTC_DEPTH'] = 0
	main_gene_df['PERC_NTC_DEPTH'] = main_gene_df['PERC_NTC_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
		

	### create output dict
	output_dict = create_output_dict(gene_list, main_gene_df, main_region_df, present, cosmic_500_df)

	## export dict to JSON
	with open(args.outfile,'w') as f:
		json.dump(output_dict, f, indent = 4, default = np_encoder)

		logger.info('json file complete')
