'''
This a variant of a script given by the COMBINE lab to convert the binary bootstrap file to a TSV.
-initially the data was given in long format(tx as columns), So I'm, using pandas to make it a little cleaner
-was writing all bootstrap values, but I'm considering just writing the variance, really depends on whether or not
 I want to  use sleuth


'''
import gzip
import struct
import argparse
import os
import logging
import logging.handlers
import sys
import errno
import json
import pandas as pd

# from: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def main(args):
    logging.basicConfig(level=logging.INFO)
    quantDir = args.quantDir
    auxDir = "aux"

    # Check for a custom auxDir
    with open(os.path.sep.join([quantDir, "cmd_info.json"])) as cmdFile:
        dat = json.load(cmdFile)
        if 'auxDir' in dat:
            auxDir = dat['auxDir']

    bootstrapFile = os.path.sep.join([quantDir, auxDir, "bootstrap", "bootstraps.gz"])
    nameFile = os.path.sep.join([quantDir, auxDir, "bootstrap", "names.tsv.gz"])
    if not os.path.isfile(bootstrapFile):
       logging.error("The required bootstrap file {} doesn't appear to exist".format(bootstrapFile))
       sys.exit(1)
    if not os.path.isfile(nameFile):
       logging.error("The required transcript name file {} doesn't appear to exist".format(nameFile))
       sys.exit(1)

    txpNames = None
    with gzip.open(nameFile) as nf:
        txpNames = nf.read().decode().strip().split('\t')

    ntxp = len(txpNames)
    logging.info("Expecting bootstrap info for {} transcripts".format(ntxp))

    with open(os.path.sep.join([quantDir, auxDir, "meta_info.json"])) as fh:
        meta_info = json.load(fh)

    if meta_info['samp_type'] == 'gibbs':
        #s = struct.Struct('<' + 'i' * ntxp)
        s = struct.Struct('@' + 'd' * ntxp)
    elif meta_info['samp_type'] == 'bootstrap':
        s = struct.Struct('@' + 'd' * ntxp)
    else:
        logging.error("Unknown sampling method: {}".format(meta_info['samp_type']))
        sys.exit(1)

    numBoot = 0
    outDir = args.outDir
    if os.path.exists(outDir):
        if os.path.isfile(outDir):
            logging.error("The requested output directory {} already exists, but is a file".format(outDir))
            sys.exit(1)
        else:
            logging.warn("The requested output directory {} already exists --- any existing bootstraps may be overwritten".format(outDir))
    else:
        mkdir_p(outDir)

    outFile = os.path.sep.join([outDir, 'quant_bootstraps.tsv.gz'])
#changes start here
    tmpfile='/tmp/garbage.smoob'
    with open(tmpfile,'w') as ofile:
        # write the header
        out_list=[]
        out_list.append(txpNames)
        # Now, iterate over the bootstrap samples and write each
        with gzip.open(bootstrapFile) as bf:
            while True:
                try:
                    x = s.unpack_from(bf.read(s.size))
                    xs =[str(i) for i in x]
                    out_list.append(xs)
                    numBoot += 1
                except:
                    logging.info("read all bootstrap values")
                    break
        full_tab=pd.DataFrame(out_list)
        print(full_tab.shape)
        full_tab.transpose().to_csv(outFile, sep='\t', header=False, index=False, compression='gzip')
    logging.info("wrote {} bootstrap samples".format(numBoot))
    logging.info("converted bootstraps successfully.")

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Convert bootstrap results to text format")
   parser.add_argument('quantDir', type=str, help="path to (sailfish / salmon) quantification directory")
   parser.add_argument('outDir', type=str, help="path to directory where results should be written")
   args = parser.parse_args()
   main(args)
