#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:23:35 2022

@author: lendler
"""

import sys
import os
import warnings
import re
import argparse
import glob
import gzip
import bz2
import lzma
from pysam import VariantFile

def open_file(filename, mode='rb'):
    """

    Parameters
    ----------
    filename : string
        Name of the file to open.
    mode : string, optional
        opening mode. The default is 'rb'.

    Returns
    -------
    filehandle, depending on mode and ending of filename (gz, xz of bz2).

    """
    if re.search("\.gz$", filename):
        return(gzip.open(filename, mode=mode))
    elif re.search("\.xz$", filename):
        return(lzma.open(filename, mode=mode))
    elif re.search("\.bz2$", filename):
        return(bz2.open(filename, mode=mode))
    else:
        return(open(filename, mode = mode))

parser = argparse.ArgumentParser(description="""
   This script takes one or multiple vcf files with format fields a list of format types (default AF and DP) and creates a long tsv with the following entries:
    SAMPLEID, CHROM, POS, REF, ALT, GENE, AA, AF, DP 
    multiple ALTs will be split to only one value per line
    Type will be AF, DP and others and taken from the FORMAT fields
    only alleles in a sample fullfiling a criteria (default AF . 0.01) are output
    for samples, you can either given an index, a name, or "" which will use all samples in each filezQ?SA<<>X< A 
    If there are annotations in an INFO ANN field, the script extracts the first one for each allele, 
    translates 3 letter AA chagnes to 1 letter, gets rid of the trailing p. and adds them in the gene and AA columns
    """)
parser.add_argument("-i", dest="vcf_file", help="vcf input file (extension .vcf or vcf.gz), file with list of files (file not ending with vcf/vcf.gz) or directory containing vcfs", required=True)
parser.add_argument("-t", "--types", dest="types", 
                    help="comma delimited list of names of value to gather (default: AF, DP)", default="AF,DP")
parser.add_argument("--cond_type", dest="cond_type", 
                    help="type string for condition (minimal allele frequency), if empty string, no check (default: AF)", default="AF")
parser.add_argument("--field", dest="field", help="get values from FORMAT of INFO field, only format implemented!", default="FORMAT")
parser.add_argument("-m","--min_af", dest="min_af", type=float, help="minimal allele frequency, default: 0.01", default=0.01)
parser.add_argument("--append", dest="append",  action='store_true', 
                    help="append values to an existing file, default: FALSE", default=False)
parser.add_argument("--sampleid", dest="sampleid", 
                    help="which sample name to take (integer to indicate position in sample list or string or empty for all samples in file), default: \"\"", default="")
parser.add_argument("-o", dest="outname", help="output file name, default: stdout", default="stdout")
parser.add_argument("--na_string", dest="na_string", help="string for None or missing values, default: NA", default="NA")
args = parser.parse_args()
#args = parser.parse_args("-i /Users/lendler/Cemm_sync/Test_VCFs/CoV_24998_S109501_samp.lofreq.vcf.gz -t AF,DP -o stdout".split())
#args = parser.parse_args("-i /Users/lendler/Cemm_sync/Test_VCFs/CoV_24995_S109562_call_dp.vcf.gz -t AF,DP -o stdout".split())

#for translating AA three letter to two
onelet = ["stop", "A", "C", "D", "E", "F",  "G", "H", "I", "K", "L", "M", "N",
          "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]
threelet = ['\*', "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
            "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
            "Val", "Trp", "Tyr", "\?\?\?"]
aa_map = dict(zip(threelet, onelet))
aa_pat = re.compile("|". join(aa_map.keys()))

val_types = args.types.split(',')
vcf_fns = list()
# check if directory and then get vcf files
if os.path.isdir(args.vcf_file):    
    vcf_fns.extend(glob.glob(args.vcf_file + '/*.vcf'))
    vcf_fns.extend(glob.glob(args.vcf_file + '/*.vcf.gz'))
    vcf_fns = set(vcf_fns)
#just one vcf file: 
elif os.path.isfile(args.vcf_file) and  re.search('\.vcf(?:.gz)*$',args.vcf_file):
    vcf_fns.append(args.vcf_file)
# file with list of files:
elif os.path.isfile(args.vcf_file):
    with open(args.vcf_file,"r") as f:
        for vcf_fn in f.readlines():
            vcf_fn = vcf_fn.rstrip()
            if os.path.isfile(vcf_fn) and re.search('\.vcf(?:.gz)*$',vcf_fn):
                vcf_fns.append(vcf_fn)
    vcf_fns = set(vcf_fns)
else:
    sys.exit()
    
if args.cond_type == "":
    args.cond_type = False
    
if len(vcf_fns) == 0:
    sys.exit(f"no vcf files found in {args.vcf_file}")
#set field to format!
args.field = "FORMAT"
out_file = False
for vcf_fn in vcf_fns:
    try:
        vcf_in = VariantFile(vcf_fn)
    except:
        warnings.warn(f"vcf reader could not open file {vcf_fn}, passing file")
        continue
    #check if fields in header
    if args.field == "FORMAT":
        vcf_keys = { key : value.number for key, value in vcf_in.header.formats.iteritems() }
    elif args.field == "INFO":
        vcf_keys = { key : value.number for key, value in vcf_in.header.info.iteritems() }
    else:
        sys.exit(f"unknown field type {args.field}, should be INFO or FORMAT")
    if len(set(val_types).intersection(vcf_keys)) == 0:
        warnings.warn(f"{args.types} not found in header {args.field} of file {vcf_fn}, passing file")
        continue
    samples = list(vcf_in.header.samples)
    if args.sampleid == "":
        sampleids = samples
    elif args.sampleid.isnumeric() and int(args.sampleid) < len(samples):
        sampleids = [ samples[int(args.sampleid)] ]
    elif args.sampleid in samples:
        sampleids = [ args.sampleid ]
    else:
        warnings.warn(f"sample {args.sampleid} not found in file {vcf_fn}, skipping!")
        continue
    ann_found = False
    if "ANN" in vcf_in.header.info.keys():
        ann_found = True
    # check if vcf_out has been opened, if not do it now
    if not out_file:
        header_str = "SAMPLEID\tCHROM\tPOS\tREF\tALT\tGENE\tAA\t"+"\t".join(val_types) 
        if args.outname == "stdout":
            out_file = sys.stdout
            print(header_str)
        elif args.append is True and os.path.exists(args.outname):
            out_file = open_file(args.outname, "at")
        else:
            args.append = False
            out_file = open_file(args.outname, "wt")
            print(header_str, file=out_file)
    for rec in vcf_in:
        #rec = next(vcf_in)
        rec_ann = dict()
        try:
            alleles = rec.alts
            alts = len(rec.alts)
        except:
            warnings.warn(f"no alts in file {vcf_fn} with entry:\n{rec}")
            continue
        if alts is None or alts == 0:       
            continue
        # get annotationsa 
        if ann_found and "ANN" in rec.info.keys():
            for anno in  rec.info['ANN']:
                anno = anno.split("|")
                if not anno[0] in rec_ann:
                    anno[10] = re.sub("^p\.","", anno[10])
                    rec_ann[anno[0]] = [anno[3], aa_pat.sub(lambda m: aa_map[re.escape(m.group(0))], anno[10])]
        for i in range(0,alts):
            for sampleid in sampleids:
                if rec_ann and alleles[i] in rec_ann.keys():
                    ann_vals = "\t".join(rec_ann[alleles[i]])
                else:
                    ann_vals = "\t"
                if args.cond_type:
                    # if threshold type defined and not there, None, or too low, skip record
                    if not args.cond_type in rec.samples[sampleid].keys():
                        continue
                    elif (( vcf_keys[ args.cond_type ] == 'A' and 
                           ( rec.samples[sampleid][args.cond_type][i] is None or
                           rec.samples[sampleid][args.cond_type][i] < args.min_af )) or 
                          ( vcf_keys[ args.cond_type ] == 1 and 
                           ( rec.samples[sampleid][args.cond_type] is None or 
                           rec.samples[sampleid][args.cond_type] < args.min_af ) )):
                        continue
                vals = []
                for val_type in val_types:
                    if not val_type in rec.samples[sampleid].keys():
                        vals.append("") 
                        # empty value if value not there
                        continue
                    elif vcf_keys[ val_type ] == 'A':
                        # if entries for each allele 
                        vals.append(str(rec.samples[sampleid][ val_type][i]))
                    elif vcf_keys[ val_type ] == 'R':
                        # if entries for each allele 
                        vals.append(str(rec.samples[sampleid][ val_type][i+1]))
                    else:
                        vals.append(str(rec.samples[sampleid][ val_type]))
                vals = "\t".join([ str( x or args.na_string) for x in vals])
                print(f'{sampleid}\t{rec.chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts[i]}\t{ann_vals}\t{vals}', file = out_file)
        
if out_file:
    out_file.close()
