#!/usr/bin/env python
# -*- coding: utf-8 -*-

from docopt import docopt
from matplotlib import pyplot as plt
from os.path import isfile

import sys
import re
import pysam
import numpy as np
import pandas as pd

multibedfile_h = pd.read_csv("feisthamelii_bedunion.txt", sep="\t")

male_id_list = ["IF_165", "IF_236", "IF_RVcoll11F366"]
female_id_list = ["IF_RVcoll12N508", "IF_RVcoll14B411"]
NA_id_list = ["IF_RVcoll17B439"]
reference_sex = "NA"


multibedfile_h["male_mean"] = multibedfile_h[male_id_list].mean(axis=1)
multibedfile_h["female_mean"] = multibedfile_h[female_id_list].mean(axis=1)

global_mean = multibedfile_h[['female_mean', 'male_mean']].mean().mean()

multibedfile_h['female_normalised'] = multibedfile_h["female_mean"]/global_mean
multibedfile_h['male_normalised'] = multibedfile_h["male_mean"]/global_mean
multibedfile_h["female_mean/male_mean"] = multibedfile_h["female_normalised"]/multibedfile_h["male_normalised"]
multibedfile_h['f_log2'] = np.log2(multibedfile_h['female_normalised'])
multibedfile_h['m_log2'] = np.log2(multibedfile_h['male_normalised'])
multibedfile_h['f/m_log2'] = np.log2(multibedfile_h['female_mean/male_mean'])
multibedfile_h['IF_165_mean_scaled_log2'] = np.log2(multibedfile_h['IF_165']/multibedfile_h['male_mean'])
multibedfile_h['IF_236_mean_scaled_log2'] = np.log2(multibedfile_h['IF_236']/multibedfile_h['male_mean'])
multibedfile_h['IF_RVcoll11F366_mean_scaled_log2'] = np.log2(multibedfile_h['IF_RVcoll11F366']/multibedfile_h['male_mean'])
multibedfile_h['IF_RVcoll12N508_mean_scaled_log2'] = np.log2(multibedfile_h['IF_RVcoll12N508']/multibedfile_h['male_mean'])
multibedfile_h['IF_RVcoll14B411_mean_scaled_log2'] = np.log2(multibedfile_h['IF_RVcoll14B411']/multibedfile_h['male_mean'])
multibedfile_h['IF_RVcoll17B439_mean_scaled_log2'] = np.log2(multibedfile_h['IF_165']/multibedfile_h['male_mean'])
multibedfile_h = multibedfile_h.replace([np.inf, -np.inf], np.nan)

kb_per_contig = multibedfile_h[multibedfile_h['f/m_log2']<-0.5]['chrom'].value_counts()*5

multibedfile_h['left_dist'] = 0
multibedfile_h.loc[multibedfile_h['f/m_log2']<-0.5, 'left_dist'] = 1

#plot male/female log 2 histogram
#multibedfile_h.hist(column = 'f/m_log2', bins = 1000)

#logy transformed male/female log 2 histogram
multibedfile_h.hist(column = 'f/m_log2', bins = 1000, log = True)

#log2 male vs female mean scatterplot
#multibedfile_h.plot.scatter(x='f_log2', y='m_log2')

#individual mean male scaled log2 plots
#multibedfile_h.hist(column = 'IF_165_mean_scaled_log2', bins = 1000)
#multibedfile_h.hist(column = 'IF_236_mean_scaled_log2', bins = 1000)
#multibedfile_h.hist(column = 'IF_RVcoll11F366_mean_scaled_log2', bins = 1000)
#multibedfile_h.hist(column = 'IF_RVcoll12N508_mean_scaled_log2', bins = 1000)
#multibedfile_h.hist(column = 'IF_RVcoll14B411_mean_scaled_log2', bins = 1000)
#multibedfile_h.hist(column = 'IF_RVcoll17B439_mean_scaled_log2', bins = 1000)

contigs_left_dist = kb_per_contig.index
contigs_all = multibedfile_h['chrom'].unique()

#print("contig\tn_windows\tcontig_size(kb)\tproportion_in_left_dist")
#for contig in contigs_left_dist:
#	right_dist_count = multibedfile_h[multibedfile_h['chrom'] == str(contig)]['left_dist'].value_counts()[0]
#	left_dist_count = multibedfile_h[multibedfile_h['chrom'] == str(contig)]['left_dist'].value_counts()[1]
#	n_windows = left_dist_count + right_dist_count
#	print(str(contig) + "\t" + str(n_windows) + "\t" + str(n_windows * 5) + "\t" + str(left_dist_count/n_windows))

#with open("left_dist_summary.tsv", "w") as f:
#        f.write("contig\tn_windows\tcontig_size(kb)\tproportion_windows_<_-0.5_log2_f/m_cov\n")
#
#for contig in contigs_all:
#	right_dist_count = multibedfile_h[multibedfile_h['chrom'] == str(contig)]['left_dist'].value_counts()[0]
#	try:
#		left_dist_count = multibedfile_h[multibedfile_h['chrom'] == str(contig)]['left_dist'].value_counts()[1]
#	except KeyError:
#		left_dist_count = 0
#	n_windows = left_dist_count + right_dist_count
#	f_line = str(contig) + "\t" + str(n_windows) + "\t" + str(n_windows * 5) + "\t" + str(left_dist_count/n_windows) + "\n"
#	with open("left_dist_summary.tsv", "a") as f:
#		f.write(f_line)

#left_dist_summary = pd.read_csv("left_dist_summary.tsv", sep="\t")
#left_dist_summary[['species','ref_id','version','contig']] = left_dist_summary['contig'].str.split('.',expand=True)
#scatter plot contig size vs proportion of half mapped
#left_dist_summary['log_contig_size(kb)'] = np.log(left_dist_summary['contig_size(kb)'])
#left_dist_summary.plot.scatter(x='proportion_in_left_dist', y='contig_size(kb)')
#left_dist_summary.plot.scatter(x='proportion_windows_<_-0.5_log2_f/m_cov', y='log_contig_size(kb)')


#print(left_dist_summary[left_dist_summary['proportion_windows_<_-0.5_log2_f/m_cov'] > 0.8])

#multibedfile_h.hist(column = 'male_mean', bins = 10000)
#multibedfile_h.hist(column = 'female_mean', bins = 10000)


#print(multibedfile_h[multibedfile_h['male_mean']<5]['chrom'].value_counts()*5)


#multibedfile_h[multibedfile_h['chrom']=='iphiclides_feisthamelii.IF_142.v1_1.ctg_116'].hist(column=['female_normalised','male_normalised'], bins = 10, sharex=True, sharey=True, layout=(2,1))
#multibedfile_h[multibedfile_h['chrom']=='iphiclides_feisthamelii.IF_142.v1_1.ctg_116'].hist(column=['f_log2','m_log2'], bins = 10, sharex=True, sharey=True, layout=(2,1))


plt.show()

