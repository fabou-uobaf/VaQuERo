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
#import bz2
#import lzma
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
    if re.search(r"\.gz$", filename):
        return(gzip.open(filename, mode=mode))
#    elif re.search("\.xz$", filename):
#        return(lzma.open(filename, mode=mode))
#    elif re.search("\.bz2$", filename):
#        return(bz2.open(filename, mode=mode))
    else:
        return(open(filename, mode=mode))


parser = argparse.ArgumentParser(description="""
   This script takes one or multiple vcf files with format or info fields a list of format types (default AF and DP) and creates a long tsv with the following entries:
    SAMPLEID, CHROM, POS, REF, ALT, GENE, AA, AF, DP, PQ
    The field used can be FORMAT (for vcfs with sample info), INFO or AUTO.
    For FORMAT the sample name is derived from the sample name in the vcf file, for INFO frpm the file name without vcf.gz.
    For AUTO (default), if a vcf contains one or more FORMAT sample fields, FORMAT is used, otherwise INFO.
    multiple ALT alleles will be split to only one allele and value per line
    Type will be AF, DP and others and taken from the FORMAT fields
    only alleles in a sample fullfiling a criteria (default AF . 0.01) are output
    The output is written to stdtout or a filename (gzipped if ending with .gz)
    For samples, you can either give an index, a name, or "" which will use all samples in each file. 
    If there are annotations in an INFO ANN field, the script extracts the first one for each allele, 
    translates 3 letter AA changes to 1 letter, gets rid of the trailing p. and adds them in the gene and AA columns
    If you take the raw lowfreq output (only one sample per file), set field to INFO and the script will only use the filename (without vcf.gz) as the sample name.
    """)
parser.add_argument("-i", dest="vcf_file", help="vcf input file (extension .vcf or vcf.gz), file with list of files (file not ending with vcf/vcf.gz) or directory containing vcfs", required=True)
parser.add_argument("-t", "--types", dest="types",
                    help="comma delimited list of names of value to gather (default: AF, DP)", default="AF,DP,PQ")
parser.add_argument("--cond_type", dest="cond_type",
                    help="type string for condition (minimal allele frequency), if empty string, no check (default: AF)", default="AF")
parser.add_argument("--field", dest="field", help="get values from FORMAT of INFO field, AUTO assumes format is vaild, if there are samples, otherwise uses INFO fields", default="AUTO")
parser.add_argument("-m", "--min_af", dest="min_af", type=float, help="minimal allele frequency, default: 0.01", default=0.01)
parser.add_argument("--append", dest="append", action='store_true',
                    help="append values to an existing file, default: FALSE", default=False)
parser.add_argument("--sampleid", dest="sampleid",
                    help="which sample name to take from vcf with format fields, not info! (integer to indicate position in sample list or string or empty for all samples in file), default: \"\"", default="")
parser.add_argument("-o", dest="outname", help="output file name, default: stdout", default="stdout")
parser.add_argument("--na_string", dest="na_string", help="string for None or missing values, default: NA", default="NA")
args = parser.parse_args()
#args = parser.parse_args("-i /Users/lendler/Cemm_sync/Test_VCFs/CoV_24998_S109501_samp.lofreq.vcf.gz -t AF,DP -o stdout".split())
#args = parser.parse_args("-i /Users/lendler/Cemm_sync/Test_VCFs/CoV_24995_S109562_call_dp.vcf.gz -t AF,DP -o stdout".split())
#args = parser.parse_args("-i /Users/lendler/Cemm_sync/Test_VCFs/sample40_variants.vcf.gz --field INFO -t AF,DP -o stdout".split())

#for translating AA three letter to two
onelet = ["stop", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
          "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]
threelet = [r'\*', "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
            "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
            "Val", "Trp", "Tyr", r"\?\?\?"]
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
elif os.path.isfile(args.vcf_file) and re.search(r'\.vcf(?:.gz)*$', args.vcf_file):
    vcf_fns.append(args.vcf_file)
# file with list of files:
elif os.path.isfile(args.vcf_file):
    with open(args.vcf_file, "r") as f:
        for vcf_fn in f.readlines():
            vcf_fn = vcf_fn.rstrip()
            if os.path.isfile(vcf_fn) and re.search(r'\.vcf(?:\.gz)*$', vcf_fn):
                vcf_fns.append(vcf_fn)
    vcf_fns = set(vcf_fns)
else:
    sys.exit()
#uppercase field
args.field = args.field.upper()
# set condition to false if none given
if args.cond_type == "":
    args.cond_type = False
if len(vcf_fns) == 0:
    sys.exit(f"no vcf files found in {args.vcf_file}")
#set field to format!
#args.field = "FORMAT"
out_file = False
for vcf_fn in vcf_fns:
    try:
        # pysam croaks on gz files, so just open gz/bgzip like this:
        vcf_in = VariantFile(open_file(vcf_fn))
    except:
        warnings.warn(f"vcf reader could not open file {vcf_fn}, passing file")
        continue
    # check if samples in vcf
    samp_vcf = [x for x in vcf_in.header.samples]
    f_field = args.field
    if f_field == "AUTO" and len(samp_vcf) > 0:
        f_field = "FORMAT"
    elif f_field == "AUTO" and len(samp_vcf) == 0:
        f_field = "INFO"
    #check if fields in header
    if f_field == "FORMAT":
        vcf_keys = {key: value.number for key, value in vcf_in.header.formats.iteritems()}
    elif f_field == "INFO":
        vcf_keys = {key: value.number for key, value in vcf_in.header.info.iteritems()}
    else:
        sys.exit(f"unknown field type {args.field}, should be INFO or FORMAT")
    if len(set(val_types).intersection(vcf_keys)) == 0:
        warnings.warn(f"{args.types} not found in header {args.field} of file {vcf_fn}, passing file")
        continue
    # set sample names
    if f_field == "FORMAT":
        samples = list(vcf_in.header.samples)
        if args.sampleid == "":
            sampleids = samples
        elif args.sampleid.isnumeric() and int(args.sampleid) < len(samples):
            sampleids = [samples[int(args.sampleid)]]
        elif args.sampleid in samples:
            sampleids = [args.sampleid]
        else:
            warnings.warn(f"sample {args.sampleid} not found in file {vcf_fn}, skipping!")
            continue
    elif f_field == "INFO":
        # for info filed, derive sample id from filename!
        sampleids = [re.sub(r"\.vcf(?:\.gz)*$", "", os.path.basename(vcf_fn))]
    ann_found = False
    if "ANN" in vcf_in.header.info.keys():
        ann_found = True
    # check if vcf_out has been opened, if not do it now
    if not out_file:
        header_str = "SAMPLEID\tCHROM\tPOS\tREF\tALT\tGENE\tAA\t" + "\t".join(val_types)
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
        # get annotations
        if ann_found and "ANN" in rec.info.keys():
            for anno in rec.info['ANN']:
                anno = anno.split("|")
                if not anno[0] in rec_ann:
                    anno[10] = re.sub(r"^p\.", "", anno[10])
                    rec_ann[anno[0]] = [anno[3], aa_pat.sub(lambda m: aa_map[re.escape(m.group(0))], anno[10])]
        for i in range(0, alts):
            for sampleid in sampleids:
                if rec_ann and alleles[i] in rec_ann.keys():
                    ann_vals = "\t".join(rec_ann[alleles[i]])
                else:
                    ann_vals = "\t"
                if f_field == "FORMAT":
                    rec_entry = rec.samples[sampleid]
                elif f_field == "INFO":
                    rec_entry = rec.info
                if args.cond_type:
                    # if threshold type defined and not there, None, or too low, skip record
                    if args.cond_type not in rec_entry.keys():
                        continue
                    elif ((vcf_keys[args.cond_type] == 'A'
                          and (len(rec_entry[args.cond_type]) <= i
                               or rec_entry[args.cond_type][i] is None
                               or rec_entry[args.cond_type][i] < args.min_af))
                          or (vcf_keys[args.cond_type] == 1
                              and (rec_entry[args.cond_type] is None
                                   or rec_entry[args.cond_type] < args.min_af))):
                        continue
                vals = []
                for val_type in val_types:
                    if val_type not in rec_entry.keys():
                        vals.append("")
                        # empty value if value not there
                        continue
                    elif vcf_keys[val_type] == 'A':
                        # if entries for each alternative allele
                        vals.append(str(rec_entry[val_type][i]))
                    elif vcf_keys[val_type] == 'R':
                        # if entries for ref and each allele
                        vals.append(str(rec_entry[val_type][i + 1]))
                    else:
                        vals.append(str(rec_entry[val_type]))
                vals = "\t".join([str(x or args.na_string) for x in vals])
                print(f'{sampleid}\t{rec.chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts[i]}\t{ann_vals}\t{vals}', file=out_file)

if out_file:
    out_file.close()
