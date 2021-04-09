#!/usr/bin/env python
# -*- coding: utf-8 -*-

from docopt import docopt
from matplotlib import pyplot
from os.path import isfile

import sys
import re
import pysam
import numpy as np
import pandas as pd

dup_df = pd.read_csv("orthofinder/Gene_Duplication_Events/Duplications.tsv", sep="\t")
#Orthogroup, species tree node, gene tree node, support, type, genes 1, genes 2

print(dup_df)