#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

orthogroup_df = pd.read_csv("orthofinder/Orthogroups/Orthogroups.tsv", sep="\t")
sp1_ortho_df = orthogroup_df[['Orthogroup', 'iphiclides_feisthamelii.pep']]
sp2_ortho_df = orthogroup_df[['Orthogroup', 'iphiclides_podalirius.pep']]

sp1_ortho_df_split = sp1_ortho_df['iphiclides_feisthamelii.pep'].str.split(', ', expand = True)
sp1_ortho_df_split['Orthogroup'] = sp1_ortho_df['Orthogroup']
sp1_ortho_df_melt = sp1_ortho_df_split.melt(id_vars = 'Orthogroup').sort_values(by=['Orthogroup', 'variable'])
sp1_ortho_df_melt.columns = ['orthogroup','protein_number','protein_name']

sp2_ortho_df_split = sp2_ortho_df['iphiclides_podalirius.pep'].str.split(', ', expand = True)
sp2_ortho_df_split['Orthogroup'] = sp2_ortho_df['Orthogroup']
sp2_ortho_df_melt = sp2_ortho_df_split.melt(id_vars = 'Orthogroup').sort_values(by=['Orthogroup', 'variable'])
sp2_ortho_df_melt.columns = ['orthogroup','protein_number','protein_name']

with open("sp1_orthogroups.txt", "w") as f:
        sp1_ortho_df_melt.mask(sp1_ortho_df_melt.eq('None')).dropna().to_csv(f, sep="\t", index=False, na_rep="NA")

with open("sp2_orthogroups.txt", "w") as f:
        sp2_ortho_df_melt.mask(sp2_ortho_df_melt.eq('None')).dropna().to_csv(f, sep="\t", index=False, na_rep="NA")