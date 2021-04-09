#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Samparse

Usage: 
 samparse.py -i <FILE> -o <FILE>

Options:
 -i, --infile <FILE>          BAM file
 -o, --outfile <FILE>         parsed BAM tsv

"""
from docopt import docopt
from matplotlib import pyplot
from os.path import isfile

import sys
import re
import pysam
import numpy as np
import pandas as pd

def parse_bam(bam_fh, parsed_bam_f="bam_parsed.tsv"):

    #convert BAM file into tsv with selected fields and simple coverage metrics

    print("parsing...")
    allowed_operations = set([0, 1, 2, 3, 4, 5])
    # 0 = match, 1 = insertion, 2 =  deletion, 3 = ref skip, 4 = soft clip, 5 = hard clip, 6 = pad, 7 = equal, 8 = diff, 9 = back

    with open(str(parsed_bam_f), "w") as f:
        f.write(
            "\t".join(
                [
                    x
                    for x in [
                        "query_name",
                        "reference_name",
                        "NM",
                        "M",
                        "I",
                        "D",
                        "N",
                        "n_N",
                        "S",
                        "H",
                        "identity",
                        "query_coverage",
                        "aln_length",
                        "ref_start",
                        "ref_stop",
                        "flag",
                        "AS",
                        "cigar",
                        "\n"
                    ]
                ]
            )
        )

    for read in bam_fh.fetch(until_eof=True):

        queryname = read.query_name

        if read.reference_name:
            referencename = read.reference_name
        else:
            referencename = "*"

        cigar_string = read.cigarstring

        if cigar_string:
            N_introns = cigar_string.count("N")
        else:
            N_introns = 0

        cigarparse = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0, "5": 0}
        if read.cigartuples:
            for operation, length in read.cigartuples:
                if operation in allowed_operations:
                    cigarparse[str(operation)] += length

        try:
            nm = read.get_tag("NM")
        except KeyError:
            nm = 0

        aln_len = cigarparse["0"] + cigarparse["1"]

        try:
            qlen = aln_len + cigarparse["4"]
        except KeyError:
            qlen = aln_len

        try:
            identity = (cigarparse["0"] - nm) / aln_len
        except ZeroDivisionError:
            identity = "NA"

        try:
            querycoverage = aln_len / qlen
        except ZeroDivisionError:
            querycoverage = "NA"

        try:
            as_tag = read.get_tag("AS")
        except KeyError:
            as_tag = "NA"

        flag = read.flag
        # secondary alignments flag 256+

        try:
            ref_start = read.get_reference_positions()[1]
        except IndexError:
            ref_start = "NA"

        try:
            ref_stop = read.get_reference_positions()[-1]
        except IndexError:
            ref_stop = "NA"

        with open(str(parsed_bam_f), "a") as f:
            f.write(
                "\t".join(
                    [
                        str(x)
                        for x in [
                            queryname,
                            referencename,
                            nm,
                            cigarparse["0"],
                            cigarparse["1"],
                            cigarparse["2"],
                            cigarparse["3"],
                            N_introns,
                            cigarparse["4"],
                            cigarparse["5"],
                            identity,
                            querycoverage,
                            aln_len,
                            ref_start,
                            ref_stop,
                            flag,
                            as_tag,
                            cigar_string,
                            "\n",
                        ]
                    ]
                )
            )

    print("done")


# need to think about multiple alignments aligning to the same place
# and the same alignment in multiple places
# ref stop-ref start != aln_len
# print blessed alignment list (id and cov > 0.9)
# can also filter query, reference, start, and stop for a pseudobed file after filtering

def append_transdecoder_blessed(parsed_bam_f, transdecoder_blessed_tsv_f, bam_transdecoder_appended_f):

    #append proteins blessed by transdecoder

    print("appending transdecoder blessed")

    parsed_bam_df = pd.read_csv(parsed_bam_f, sep="\t", encoding = 'utf8')

    if transdecoder_blessed_tsv_f:
        transdecoder_blessed_df = pd.read_csv(transdecoder_blessed_tsv_f, sep="\t")
        merged_df = pd.merge(parsed_bam_df, transdecoder_blessed_df, on="query_name", how="outer")
        with open(bam_transdecoder_appended_f, "w") as f:
            merged_df.to_csv(f, sep="\t", index=False)

    #figure out how to drop unnamed columns        

    print("done")

def append_expression_sleuthed(bam_transdecoder_appended_f, TPM_f, EST_f, DE_f, bam_expression_appended_f):

    # append expression tables obtained after processing sleuth output in R
    # TPM count is TPM table
    # EST counts is a count table
    # DE genes differentially expressed genes
    # sexantag file has appended info from other files

    # replace x with tpm and y with est in header
    # test for correct merge
    # remove bits not needed
    # do calculations with what is needed
    # save output
    # append to previous by transcript
    # rename for merging

    print("appending expression sleuthed")

    bam_transdecoder_appended_df = pd.read_csv(bam_transdecoder_appended_f, sep = "\t")
    TPM_df = pd.read_csv(TPM_f, sep=" ", index_col=False)
    EST_df = pd.read_csv(EST_f, sep=" ", index_col=False)
    DE_df = pd.read_csv(DE_f, sep=" ")
    DE_df = DE_df.rename(columns={"target_id": "Transcript"})

    temp = pd.merge(TPM_df, EST_df, on="Transcript", how="outer")
    expression_sleuthed_df = pd.merge(temp, DE_df, on="Transcript", how="outer")

    expression_sleuthed_df = expression_sleuthed_df.rename(columns={'Transcript': 'query_name'})
    bam_expression_appended_df = pd.merge(bam_transdecoder_appended_df, expression_sleuthed_df, on="query_name", how="outer")

    with open(bam_expression_appended_f, "w") as f:
        bam_expression_appended_df.to_csv(f, sep="\t", index=False)

    print("done")

def append_orthogroups():
    #append orthogroups between species and its pair, and labelled single copy orthologs
    print("appending orthogroups")

def append_orthofinder_putative_duplications():
    #append orthofinder putative duplications by protein
    print("appending putatitve duplications")

def quality_control_filtering(df_fh):
    #filtering and quality control processing

    print("filtering and qc")

    df = pd.read_csv(df_fh, sep="\t")

    #remove unaligned query names
    df = df[["reference_name"] != "*"]

    #write file with no unaligned queries
    with open("samparse_output_nofails.tsv", "w") as f:
        df.to_csv(f, sep="\t", index=False)

    #write list of unaligned queries    
    with open("unaligned_query_names.txt", "w") as f:
        df[df["reference_name"] == "*"]["query_name"].to_csv(f, index=False)

    # save unique isoforms that successfully mapped
    np.savetxt("uniq_query_isoform_mapped.txt", samparse_df["query_name"].unique(), fmt="%s")

    print("done")

    # print files for
    # queries with no best transdecoder
    # queries with transdecoder output
    # transdecoder proteins with no queries
    # could add a tally to isoform count

bam_f = "minimap2/iphiclides_podalirius.IP_504.transcriptome.bam"
parsed_bam_f = "working_dir/bam_parsed.tsv"
transdecoder_blessed_tsv_f = "transdecoder/transdecoder_blessed.tsv"
bam_transdecoder_appended_f = "working_dir/bam_transdecoder_appended.tsv"
TPM_f = "kallisto/Kallisto_TPM_table_iphiclides.txt"
EST_f = "kallisto/Kallisto_ESTcounts_table_iphiclides.txt"
DE_f = "kallisto/DEgenes_iphiclides.txt"
bam_expression_appended_f = "working_dir/bam_expression_appended.tsv"

#with pysam.AlignmentFile(bam_f, "rb") as bam_fh:
#    parse_bam(bam_fh, parsed_bam_f)
#append_transdecoder_blessed(parsed_bam_f, transdecoder_blessed_tsv_f, bam_transdecoder_appended_f)
#append_expression_sleuthed(bam_transdecoder_appended_f, TPM_f, EST_f, DE_f, bam_expression_appended_f)
#append_orthogroups()
#quality_control_filtering(df_f)



#docopt code

# def main():
#    args = docopt(__doc__)
#    infile = args['--infile']
#    outfile = args['--outfile']
#
#    with pysam.AlignmentFile(str(infile),"rb") as bamfile:
#        samparse(bamfile, str(outfile))
#
# if __name__ == '__main__':
#    main()