from docopt import docopt
from matplotlib import pyplot
from os.path import isfile

import sys
import re
import pysam
import numpy as np
import pandas as pd

#natsort, readBed, BedObj, and MainObj create transdecoder blessed tsv
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
                if status == "complete":
                    if protein_length >= MIN_LENGTH:
                        if score >= MIN_SCORE:
                            bedObj = BedObj(
                                isoform_id,
                                protein_id,
                                status,
                                protein_length,
                                cds_length,
                                score,
                                strand,
                                start,
                                end,
                            )
                            yield bedObj


class BedObj:
    def __init__(
        self,
        isoform_id,
        protein_id,
        status,
        protein_length,
        cds_length,
        score,
        strand,
        start,
        end,
    ):
        self.isoform_id = isoform_id
        self.gene_id = "_".join(isoform_id.split("_")[0:4])
        self.protein_id = protein_id
        self.status = status
        self.protein_length = protein_length
        self.cds_length = cds_length
        self.score = score
        self.strand = strand
        self.start = start
        self.end = end

    def bed(self):
        return "\t".join(
            [
                str(x)
                for x in [
                    self.isoform_id,
                    self.start,
                    self.end,
                    self.gene_id,
                    self.score,
                    self.strand,
                    self.cds_length,
                    self.protein_id,
                ]
            ]
        )


class MainObj:
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
            for idx, bedObj in enumerate(
                sorted(
                    self.bedObjs_by_gene_id[gene_id],
                    key=lambda x: (x.score),
                    reverse=True,
                )
            ):
                if idx >= 1:
                    break
                with open("transdecoder_blessed.tsv", "a") as f:
                    f.write(bedObj.bed() + "\n")


# maybe convert function to pandas data frame and then combine with existing data before writing that as an option
#would be more contiguous with the rest of the script

_nsre = re.compile("([0-9]+)")
BED = "./transdecoder/iphiclides_podalirius.trinity.Trinity.fasta.transdecoder.bed"
MIN_SCORE = float(0)
MIN_LENGTH = int(30)
with open("transdecoder_blessed.tsv", "w") as f:
    f.write(
        "\t".join(
            [
                x
                for x in [
                    "query_name",
                    "start_on_gene",
                    "stop_on_gene",
                    "gene_id",
                    "score",
                    "strand",
                    "cds_length",
                    "protein_id",
                    "\n",
                ]
            ]
        )
    )

MainObj()

# if __name__ == "__main__":
#    _nsre = re.compile('([0-9]+)')
#    args = docopt(__doc__)
#
#    BED = args['--transdecoder_bed']
#    MIN_SCORE = float(args['--score'])
#    MIN_LENGTH = int(args['--length'])
#
#    MainObj()
