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

protein_name_df = pd.read_csv("protein_name_list_feisthamelii.txt")
protein_name_df[['isoform','protein']] = protein_name_df['isoform'].str.split('.',expand=True)
protein_name_df[['Trinity','DN','cluster','gene','isoform']] = protein_name_df['isoform'].str.split('_',expand=True)
protein_name_df['cluster_gene'] = protein_name_df['Trinity'] + '_' + protein_name_df['DN'] + '_' + protein_name_df['cluster'] + '_' + protein_name_df['gene']
protein_count_df = protein_name_df.groupby(['cluster_gene','isoform']).size().to_frame(name = 'n_proteins').reset_index()
isoform_dist_df = protein_count_df.groupby(['cluster_gene']).size().to_frame(name = 'n_isoforms').reset_index()
print(isoform_dist_df)
bin_size = isoform_dist_df['n_isoforms'].max()

#isoform_dist_df.hist(column = 'n_isoforms', bins = bin_size)
#isoform_dist_df.hist(column = 'n_isoforms', bins = bin_size, log=True)

#plt.show()

#add sex bias factor
#get isoform lengths

