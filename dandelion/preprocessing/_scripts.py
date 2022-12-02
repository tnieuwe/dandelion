#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-12-02 18:50:42

import os
import pandas as pd
import shutil

from pathlib import Path
from tqdm import tqdm
from typing import List, Optional

from dandelion.utilities._io import fasta_iterator, Write_output


def prepare_non10x_fasta(
    fasta: str,
    outdir: Optional[str] = None,
):
    """
    Prepare a non-10x fasta so that it can be ingested like for downstream analysis.

    Parameters
    ----------
    fasta : str
        path to fasta file.
    outdir : Optional[str], optional
        path to output location. `None` defaults to 'dandelion'.
    """
    fh = open(fasta, "r")
    seqs = {}
    for header, sequence in fasta_iterator(fh):
        seqs[header] = sequence
    fh.close()
    basedir = os.path.dirname(fasta)
    if outdir is None:
        outdir = basedir.rstrip("/") + "/" + Path(os.path.basename(fasta)).stem
    if not outdir.endswith("/"):
        outdir = outdir + "/"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out_fasta = outdir + "all_contig.fasta"
    out_anno_path = outdir + "all_contig_annotations.csv"
    fh1 = open(out_fasta, "w")
    fh1.close()
    out = ""
    anno = []
    for l in seqs:
        out = ">" + l + "-1_contig-1" + "\n" + seqs[l] + "\n"
        Write_output(out, out_fasta)
        # also create a dummy contig_annotations.csv
        defaultrow = {
            "barcode": l,
            "is_cell": "TRUE",
            "contig_id": l + "-1_contig-1",
            "high_confidence": "TRUE",
            "length": str(len(seqs[l])),
            "chain": "None",
            "v_gene": "None",
            "d_gene": "None",
            "j_gene": "None",
            "c_gene": "None",
            "full_length": "TRUE",
            "productive": "TRUE",
            "cdr3": "None",
            "cdr3_nt": "None",
            "reads": 1,
            "umis": 1,
            "raw_clonotype_id": "None",
            "raw_consensus_id": "None",
        }
        anno.append(defaultrow)
    anno = pd.DataFrame(anno)
    anno.to_csv(out_anno_path, index=False)


def prepare_non10x_fastas(
    fastas: List[str],
    outdir: Optional[str] = None,
):
    """
    Prepare a non-10x fastas so that it can be ingested like for downstream analysis.

    Parameters
    ----------
    fastas : List[str]
        list of paths to fasta files.
    outdir : Optional[str], optional
        path to out put location.
    """
    if type(fastas) is not list:
        fastas = [fastas]

    for i in tqdm(
        range(0, len(fastas)),
        desc="Formating fasta(s) ",
        bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}",
    ):
        prepare_non10x_fasta(
            fastas[i],
            outdir=outdir,
        )
