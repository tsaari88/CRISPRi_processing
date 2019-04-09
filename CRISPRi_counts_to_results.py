#!/usr/bin/env python3

import argparse
import os
import re
import json
import numpy as np #Tested w/ v 1.16.2
import pandas as pd #Tested w/ v 0.24
from scipy.stats.mstats import gmean


def process_all(file_df, thres, norm_meth, df_to_merge, outdir):
    '''
    Process all replicates within an experiment, merge these results, and write to file
    '''
	expt_groups = file_df.groupby('expt_name')
	#Process each experiment
	for expt_name in expt_groups.groups.keys():
		#Get rows from file_df for a given experiment
		expt = file_df.iloc[expt_groups.groups[expt_name]].reset_index()
		#Now process each replicate in the expt
		results_reps = []
		rep_groups = expt.groupby('replicate')
		#Gather results from processed replicates in a list
		for rep in rep_groups.groups.keys():
			rep_inpt_df = expt.iloc[rep_groups.groups[rep]].reset_index()
			results_reps.append(process_replicate(rep_inpt_df, rep, thres, norm_meth))
		#Start with merging first rep with df_to_merge, then merge the other experimental reps iteratively
		expt_merge = pd.merge(df_to_merge, results_reps[0], how='left', on='OligoID')
		for i in range(1,len(results_reps)):
			expt_merge = expt_merge.merge(results_reps[i], how='outer', on='OligoID')
		#Calculate mean log2FC from replicates
		log2FC_cols = [col for col in expt_merge.columns if re.search("log2FC", col)]
		expt_merge['mean_log2FC'] = expt_merge[log2FC_cols].mean(axis=1)
		#Write to file
		fpath = os.path.join(outdir, expt_name + "_results.txt")
		expt_merge.to_csv(fpath, sep="\t", index=False)


def process_replicate(file_df, rep, thres, norm_meth):
    '''
    Process a single replicate of an experiment (calclate log2FC)
    '''
	# "".join(list()) is used to get a string from a single-valued Series
	start_countsfile = "".join(list(file_df.query('timepoint == "START"')['countsfilepath']))
	end_countsfile = "".join(list(file_df.query('timepoint == "END"')['countsfilepath']))
	start_count_col = "rep" + rep + "_START"
	start_df = pd.read_csv(start_countsfile, sep="\t", header=None, names=['OligoID', start_count_col], dtype={'OligoID': str, start_count_col: np.int32})
	#Remove those with T0 counts lower than threshold
	start_df.drop(start_df.loc[start_df[start_count_col] < thres].index, inplace=True)
	end_count_col = "rep" + rep + "_END"
	end_df = pd.read_csv(end_countsfile, sep="\t", header=None, names=['OligoID', end_count_col], dtype={'OligoID': str, end_count_col: np.int32})
	#Merge and normalize
	rep_df = pd.merge(start_df, end_df, how='inner', on='OligoID')
	rep_df = normalize_counts(rep_df, norm_meth)
	#Setup for calculation
	colname_log2FC = "rep" + rep + "_log2FC"
	#get appropriate start and end col via regex
	start_normcol_l = [col for col in rep_df.columns if re.search("START.*Norm", col)]
	start_normcol = "".join(start_normcol_l) if len(start_normcol_l) == 1 else "" #Should only return 1 col name, otherwise init empty string to throw error
	end_normcol_l = [col for col in rep_df.columns if re.search("END.*Norm", col)]
	end_normcol = "".join(end_normcol_l) if len(end_normcol_l) == 1 else ""
	#Calculate log2FC from normalized start/end values
	rep_df[colname_log2FC] = np.log2(rep_df[end_normcol] / rep_df[start_normcol])
	return(rep_df)


def normalize_counts(df, method):
    '''
    Normalize counts with either total readcount or median-ratio method
    '''
	count_slice = df.iloc[:,1:3] #df should only have 3 columns: OligoID, start counts, and end counts. Grab the latter two
	if method == "readcount":
		normed = count_slice / count_slice.apply(np.sum, axis=0) * 1000000
		normed.rename(lambda x: x + "_rcNorm", axis=1, inplace=True)
		return(pd.concat((df, normed), axis=1))
	elif method == "median_ratio":
		#Median ratio normalization: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106
		gm_normed = count_slice / count_slice.apply(gmean, axis=0)
		median_ratio_normed = count_slice / gm_normed.apply(np.median, axis=0)
		median_ratio_normed.rename(lambda x: x + "_mrNorm", axis=1, inplace=True)
		return(pd.concat((df, median_ratio_normed), axis=1))




if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='python CRISPRi_counts_to_results.py', description = "Combines counts files in an experiment, performs initial filtering, read count normalization, and guide-level CiScore calculation")
	parser.add_argument('-m', '--T0_mincount', default=50, help="Minimum count in T0 datapoint required to use for processing")
	parser.add_argument('-o', '--output_dir', required=True, help="Directory to place the results. Each experiment will have its own results file in this directory with the naming convention expt_name_results.txt")
	merge_with = parser.add_mutually_exclusive_group(required=True)
	merge_with.add_argument('-i', '--IDs_list', help="Filename of a list of unique OligoIDs to merge results upon")
	merge_with.add_argument('-d', '--data_file', help="Filename of a data file, with chr, start, end, OligoID, GuideSequence, and target columns (for compatibility with J. Engreitz CalculateTilingStatistic.R) to merge results upon. Each row must feature a unique OligoID.")
	input_data = parser.add_mutually_exclusive_group(required=True)
	input_data.add_argument('-c', '--input_csv', help="Filename of a comma-separated value table containing experiment name, count filenames, replicates, timepoints, etc. Multiple experiments can be defined in this table, so the same table can be used for multiple experiments")
	input_data.add_argument('-j', '--input_json', help="JSON string containing experiment name, sample filenames, replicates, timepoints, etc. Multiple experiments can be defined in this table, so the same table can be used for multiple experiments")
	norm_method = parser.add_mutually_exclusive_group(required=True)
	norm_method.add_argument('-r', '--median_ratio_normalize', action="store_true", help="Use median ratio normalization, the type of normalization used by the DEseq2, etc. Default.")
	norm_method.add_argument('-R', '--readcount_normalize', action="store_false", dest="median_ratio_normalize", help="Use readcount normalization - counts per million reads. Overrides -r")

	args = parser.parse_args()

	if args.median_ratio_normalize:
		norm_method = "median_ratio"
	else:
		norm_method = "readcount"

	if args.input_csv:
		dat = pd.read_csv(args.input_csv, dtype=str)
	elif args.input_json:
		dat = pd.read_json(args.input_json, dtype=str)

	dat_required_cols = ['expt_name', 'countsfilepath', 'replicate', 'timepoint']
	if not all(x in dat.columns.tolist() for x in dat_required_cols):
		raise RuntimeError("Input DataFrame must contain the following columns: " + ", ".join(dat_required_cols))
	#Set up df to merge upon
	if args.IDs_list:
		#Minimal - OligoIDs only
		df_to_merge = pd.read_csv(args.IDs_list, header=None, names=["OligoID"])
	elif args.data_file:
		#Set up df Compatible with J. Engreitz CalulateTilingStatistic.R and/or downstream bedfile creation
		tomerge_required_cols = ['chr', 'start', 'end', 'OligoID', 'GuideSequence', 'target']
		df_to_merge = pd.read_csv(args.data_file, sep='\t', dtype=object)
		if not all(x in df_to_merge.columns.tolist() for x in tomerge_required_cols):
			raise RuntimeError("input -d/--data_file must contain the following columns: " + ", ".join(tomerge_required_cols))

	process_all(dat, int(args.T0_mincount), norm_method, df_to_merge, args.output_dir)
