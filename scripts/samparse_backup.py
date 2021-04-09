#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Samparse

Usage: 
 samparse.py -i <FILE> -o <FILE>

Options:
 -i, --infile <FILE>          BAMFILE
 -o, --outfile <FILE>         Output tsv

"""
from docopt import docopt
from matplotlib import pyplot
from os.path import isfile

import sys
import re
import pysam
import numpy as np
import pandas as pd

def samparse(bamfile, outfile="samparse_outfile.tsv"):
    
    print("Parsing...")
    allowed_operations = set([0, 1, 2, 3, 4, 5])
    #0 = match, 1 = insertion, 2 =  deletion, 3 = ref skip, 4 = soft clip, 5 = hard clip, 6 = pad, 7 = equal, 8 = diff, 9 = back
    
    with open(str(outfile), "w") as f:
        #f.write("query_name\treference_name\tNM\tM\tI\tD\tN\tn_N\tS\tH\tidentity\tquery_coverage\taln_length\tref_start\tref_stop\tflag\tAS\tcigar\n")
        f.write("\t".join([x for x in ["query_name", "reference_name", "NM", "M", "I", "D", "N", "n_N", "S", "H", 
        "identity", "query_coverage", "aln_length", "ref_start", "ref_stop", "flag", "AS", "cigar"]]))
        
    for read in bamfile.fetch(until_eof=True):
        
        queryname = read.query_name
        
        if read.reference_name:
            referencename = read.reference_name
        else:
            referencename="*"
        
        cigar_string = read.cigarstring
        
        if cigar_string:
            N_introns = cigar_string.count('N')
        else:
            N_introns = 0

        cigarparse={'0':0,'1':0,'2':0,'3':0,'4':0, '5':0}
        if read.cigartuples:
            for operation,length in read.cigartuples:
                if operation in allowed_operations:
                    cigarparse[str(operation)] += length  

        try:
            nm = read.get_tag("NM")
        except KeyError:
            nm = 0
        
        aln_len = cigarparse['0'] + cigarparse['1']
            
        try:    
            qlen = aln_len + cigarparse['4']
        except KeyError:
            qlen = aln_len
        
        try:
            identity = (cigarparse['0'] - nm) / aln_len
        except ZeroDivisionError:
            identity = "NA"
    
        try:
            querycoverage = aln_len/qlen
        except ZeroDivisionError:
            querycoverage = "NA"
            
        try:
            as_tag = read.get_tag("AS")
        except KeyError:
            as_tag = "NA"
        
        flag = read.flag
        #secondary alignments flag 256+
        
        try:
            ref_start = read.get_reference_positions()[1]
        except IndexError:
            ref_start = "NA"
            
        try:
            ref_stop = read.get_reference_positions()[-1]
        except IndexError:
            ref_stop = "NA"
        
        with open(str(outfile), "a") as f:
            f.write("\t".join([str(x) for x in [queryname, referencename, nm, cigarparse['0'], cigarparse['1'], 
                    cigarparse['2'], cigarparse['3'], N_introns, cigarparse['4'], cigarparse['5'], identity, 
                    querycoverage, aln_len,  ref_start, ref_stop, flag, as_tag, cigar_string, "\n"]]))
        
    print("Done")
    

#need to think about multiple alignments aligning to the same place
#and the same alignment in multiple places
#ref stop-ref start != aln_len
#print blessed alignment list (id and cov > 0.9)
#can also filter query, reference, start, and stop for a ghetto bed file after filtering

#    def bed(self):
#        return "\t".join([str(x) for x in [self.isoform_id, self.start, self.end, self.gene_id, self.score, self.strand, self.cds_length, self.protein_id]])

def samplots(samparse_tsv):
    
    print("Plotting...")
    
    samparse_df = pd.read_csv(samparse_tsv, sep = "\t")

    #remove unaligned query names
    samparse_df = samparse_df[samparse_df['reference_name'] != '*']

    #log log length plot    
    samparse_df['log_aln_length'] = np.log10(samparse_df['aln_length'])
    pyplot.figure(figsize=(12,12))
    pyplot.hist(samparse_df['log_aln_length'], log = True, histtype= 'step')
    pyplot.xlabel("log aln length")
    pyplot.ylabel("log count")
    pyplot.savefig("log_aln_length.png")

    #NM plot
    pyplot.figure(figsize=(12,12))
    pyplot.hist(samparse_df[['reference_name', 'NM']].groupby('reference_name').mean()['NM'],
                         log = True, histtype = 'step')
    pyplot.xlabel("Mean mismatches per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("NM_per_contig_mean.png")

    #M plot
    pyplot.figure(figsize=(12,12))
    pyplot.hist(samparse_df[['reference_name', 'M']].groupby('reference_name').mean()['M'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Mean matches per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("M_per_contig_mean.png")
    
    #i plot
    pyplot.figure(figsize=(12,12))
    pyplot.hist(samparse_df[['reference_name', 'I']].groupby('reference_name').mean()['I'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Mean insertions per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("I_per_contig_mean.png")

    #d plot
    pyplot.figure(figsize=(12,12))
    pyplot.hist(samparse_df[['reference_name', 'D']].groupby('reference_name').mean()['D'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Mean deletions per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("D_per_contig_mean.png")
    
    #N mean plot
    pyplot.figure(figsize=(12,12))
    d_max_contig = pyplot.hist(samparse_df[['reference_name', 'N']].groupby('reference_name').mean()['N'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Mean intron length per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("N_per_contig_mean.png")
    
    #N max plot
    pyplot.figure(figsize=(12,12))
    d_max_contig = pyplot.hist(samparse_df[['reference_name', 'N']].groupby('reference_name').max()['N'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Max intron length per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("N_per_contig_max.png")
    
    #mean identity per contig plot
    pyplot.figure(figsize=(12,12))
    identity_contig = pyplot.hist(samparse_df[['reference_name', 'identity']].groupby('reference_name').mean()['identity'],
                    log = True, histtype = 'step')
    pyplot.xlabel("Mean identity per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("identity_per_contig_mean.png")

    #mean query coverage per contig plot
    pyplot.figure(figsize=(12,12))
    query_coverage_contig = pyplot.hist(samparse_df[['reference_name', 'query_coverage']].groupby('reference_name').mean()['query_coverage'],
                        log = True, histtype = 'step')
    pyplot.xlabel("Mean query coverage per contig")
    pyplot.ylabel("log count")
    pyplot.savefig("qcov_per_contig_mean.png")

    print("Done")

def natsort(l):

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split(_nsre, key)]

    return sorted(l, key=alphanum_key)

def readBed(infile):
    if not isfile(infile):
        sys.exit("[X] %s is not a file." % (infile))
    with open(infile) as fh:
        for l in fh:
            if not l.startswith("track"):
                field = l.split()
                isoform_id = field[0]
                protein_id = field[3].split("~")[2].split(";")[0]
                status = field[3].split(";")[2].split(":")[1].split("_")[0]
                protein_length = int(field[3].split(":")[2].split("_")[0].split("_")[0])
                cds_length = protein_length * 3
                score = float(field[3].split("=")[2].split(",")[0])
                strand = field[5]
                start = int(field[6])
                end = int(field[7])
                if status == 'complete':
                    if protein_length >= MIN_LENGTH:
                        if score >= MIN_SCORE:
                            bedObj = BedObj(isoform_id, protein_id, status, protein_length, cds_length, score, strand, start, end)
                            yield bedObj

class BedObj():
    def __init__(self, isoform_id, protein_id, status, protein_length, cds_length, score, strand, start, end):
        self.isoform_id = isoform_id
        self.gene_id = "_".join(isoform_id.split('_')[0:4])
        self.protein_id = protein_id
        self.status = status
        self.protein_length = protein_length
        self.cds_length = cds_length
        self.score = score
        self.strand = strand
        self.start = start
        self.end = end

    def bed(self):
        return "\t".join([str(x) for x in [self.isoform_id, self.start, self.end, self.gene_id, self.score, self.strand, self.cds_length, self.protein_id]])

class MainObj():
    def __init__(self):
        self.bedObjs_by_gene_id = {}
        self.parse_bed()
        self.write_bed()

    def parse_bed(self):
        for bedObj in readBed(BED):
            # print bedObj.__dict__
            try:
                self.bedObjs_by_gene_id[bedObj.gene_id].append(bedObj)
            except KeyError:
                self.bedObjs_by_gene_id[bedObj.gene_id] = [bedObj]

    def write_bed(self):
        for gene_id in natsort(self.bedObjs_by_gene_id):
            for idx, bedObj in enumerate(sorted(self.bedObjs_by_gene_id[gene_id], key=lambda x: (x.score), reverse=True)):
                if idx >= 1:
                    break
                with open('transdecoder_blessed.tsv', "a") as f:
                    f.write(bedObj.bed() + "\n")
                    
# maybe write function to pandas data frame and then combine with existing data before writing that as an option

def postprocess(samparse_tsv, transdecoder_tsv):

    print("postprocessing...")
    
    samparse_df = pd.read_csv(samparse_tsv, sep = "\t")
    
    with open("unaligned_query_names.txt", "w") as f:
        samparse_df[samparse_df['reference_name'] == '*']['query_name'].to_csv(f, index=False)
    
    #remove unaligned query names
    samparse_df = samparse_df[samparse_df['reference_name'] != '*']

    with open("samparse_output_nofails.tsv", "w") as f:
        samparse_df.to_csv(f, sep = "\t", index=False)
    
    if transdecoder_tsv:
        transdecoder_blessed_df = pd.read_csv("transdecoder_blessed.tsv", sep = "\t")
        merged_df = pd.merge(samparse_df, transdecoder_blessed_df, on = 'query_name', how = 'outer')
        with open("samparse_output_nofails_transdecoder.tsv", "w") as f:
            samparse_df.to_csv(f, sep = "\t", index=False)
    
    #save unique isoforms that successfully mapped
    np.savetxt("uniq_query_isoform_mapped.txt", samparse_df["query_name"].unique(), fmt='%s')
    
    samparse_df = pd.read_csv("samparse_output_nofails.tsv", sep = "\t", index_col=False)
    unaligned_query_names = samparse_df[samparse_df['reference_name'] == '*']['query_name']
    samparse_df = samparse_df[samparse_df['reference_name'] != '*']

    transdecoder_blessed_df = pd.read_csv("transdecoder_blessed.tsv", sep = "\t")
    merged_df_1 = pd.merge(samparse_df, transdecoder_blessed_df, on = 'query_name', how = 'outer')

    with open("samparse_output_nofails_transdecoder.tsv", "w") as f:
        merged_df_1.to_csv(f, sep = "\t", index=False, na_rep = "NA")
    
    #print files for 
    #queries with no best transdecoder   
    #queries with transdecoder output
    #transdecoder proteins with no queries

    #target tables obtained after processing sleuth output in R
    # EST counts is a count table
    # TPM count is TPM table
    # DE genes differentially expressed genes
    # sexantag file has appended info from other files

    sexantag_TPM = pd.read_csv("kallisto/Kallisto_TPM_table_gonepteryx.txt", sep = " ", index_col=False)
    sexantag_EST = pd.read_csv("kallisto/Kallisto_ESTcounts_table_gonepteryx.txt", sep = " ", index_col=False)
    sexantag_DE = pd.read_csv("kallisto/DEgenes_gonepteryx.txt", sep = " ")
    sexantag_DE = sexantag_DE.rename(columns={'target_id' : 'Transcript'})

    #replace x with tpm and y with est in header
    #test for correct merge
    #remove bits not needed
    #do calculations with what is needed
    #save output
    #append to previous by transcript
    #rename for merging

    temp = pd.merge(sexantag_TPM, sexantag_EST, on = 'Transcript', how = 'outer')
    sexantag_merge = pd.merge(temp, sexantag_DE, on = 'Transcript', how = 'outer')
    

        
    print("Done")
    
    #could add a tally to isoform count

_nsre = re.compile('([0-9]+)')
outfile = "samparse_output.tsv"
with pysam.AlignmentFile("minimap2/gonepteryx_rhamni.GR_164.transcriptome.bam","rb") as bamfile:
    samparse(bamfile, outfile)
    
BED = './transdecoder/gonepteryx_rhamni.trinity.Trinity.fasta.transdecoder.bed'
MIN_SCORE = float(0)
MIN_LENGTH = int(30)
with open('transdecoder_blessed.tsv', "w") as f:
    f.write("query_name\tstart_on_gene\tstop_on_gene\tgene_id\tscore\tstrand\tcds_length\tprotein_id" + "\n")
MainObj()
    
#samplots(outfile)
postprocess(outfile, 'transdecoder_blessed.tsv')

#def main():
#    args = docopt(__doc__)
#    infile = args['--infile']
#    outfile = args['--outfile']
#    
#    with pysam.AlignmentFile(str(infile),"rb") as bamfile:
#        samparse(bamfile, str(outfile))
#
#if __name__ == '__main__':
#    main()
#    
#if __name__ == "__main__":
#    _nsre = re.compile('([0-9]+)')
#    args = docopt(__doc__)
#
#    BED = args['--transdecoder_bed']
#    MIN_SCORE = float(args['--score'])
#    MIN_LENGTH = int(args['--length'])
#
#    MainObj()