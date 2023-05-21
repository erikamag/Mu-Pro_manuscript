#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import sys
import logging
import numpy as np
import pandas as pd

from jcvi.apps.base import sh, mkdir
from jcvi.formats.base import must_open

def main(args):
    yid = args.project
    fi = f'{args.excel}/{yid}.xlsx'
    if op.isfile(fi):
        sh(f"excel.py tsv {fi} tmp.tsv")
        cvts = dict(SampleID=str,Tissue=str,Genotype=str,r0=str,r1=str,r2=str)
        sl = pd.read_csv('tmp.tsv', sep="\t", header=0, converters=cvts)
        for i in range(len(sl)):
            sid, paired = sl['SampleID'][i], sl['paired'][i]
            r0, r1, r2= sl['r0'][i], sl['r1'][i], sl['r2'][i]
            if args.download:
                r0n = f"{sid}_R0.fq.gz"
                r1n = f"{sid}_R1.fq.gz"
                r2n = f"{sid}_R2.fq.gz"
                isfile0 = op.isfile(f"{args.barn}/{yid}/{r0n}")
                isfile1 = op.isfile(f"{args.barn}/{yid}/{r1n}")
                isfile2 = op.isfile(f"{args.barn}/{yid}/{r2n}")
                isfile12 = isfile1 and isfile2
                if (paired == 'SE' and not isfile0) or (paired == 'PE' and not isfile12):
                    cmd = f"nf.fastq.py --source {args.source} --paired {paired} --cpu {args.cpu} --mem {args.mem} --tmp {args.tmp} {sid} --r0 '{r0}' --r1 '{r1}' --r2 '{r2}'"
                    sh(cmd)
                    mkdir(f"{args.barn}/{yid}")
                    sh(f"mv {sid}*.fq.gz {args.barn}/{yid}")
                else:
                    logging.info(f"{sid} already downloaded")
                r0, r1, r2 = r0n, r1n, r2n
            if paired == 'SE':
                f0 = f"{args.barn}/{yid}/{r0}"
                if not op.isfile(f0):
                    f0 = f"{args.barn}/{yid}/{sid}_R0.fq.gz"
                    if not op.isfile(f0):
                        logging.warning(f"{f0} not found")
                sl.at[i,'r0'] = f0
            else:
                f1 = f"{args.barn}/{yid}/{r1}"
                f2 = f"{args.barn}/{yid}/{r2}"
                if not op.isfile(f1):
                    f1 = f"{args.barn}/{yid}/{sid}_R1.fq.gz"
                    if not op.isfile(f1):
                        logging.warning(f"{f1} not found")
                if not op.isfile(f2):
                    f2 = f"{args.barn}/{yid}/{sid}_R2.fq.gz"
                    if not op.isfile(f2):
                        logging.warning(f"{f2} not found")
                sl.at[i,'r1'] = f1
                sl.at[i,'r2'] = f2
        sl.to_csv(args.output, sep="\t", header=True, index=False)
        sh("rm tmp.tsv")
    else:
        sys.exit(f"cannot read from excel [{fi}]")

if __name__ == "__main__":
    import argparse
    ps = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "prepare sequence (fastq) meta table to run nextflow"
    )

    ps.add_argument('project', help = 'project ID')
    ps.add_argument('--output', '-o', default='design.tsv', help = 'output file')
    ps.add_argument('--excel', default=f'{os.environ["proj"]}/barn/data/05_input', help='input excel')
    ps.add_argument('--barn', default=f'{os.environ["s3"]}/zhoup-barn', help='s3 barn directory')
    ps.add_argument('--download', action="store_true", help='download from sra/s3?')
    ps.add_argument('--source', default='sra', choices=['sra','s3'], help='sequence source')
    ps.add_argument('--cpu', type=int, default=1, help='number processors/threads to use')
    ps.add_argument('--mem', type=int, default=20, help='size of memeroy (in GBs) to use')
    ps.add_argument('--tmp', default=os.environ['NXF_TEMP'], help='temporary directory to use')

    args = ps.parse_args()
    main(args)

