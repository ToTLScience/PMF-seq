'''
Deconvolutes MiSeq sequencing files
Inputs: 1. Path to folder containing read and barcode fastq files
        2. Path to reference file (.csv with spacer sequence, spacer name)
        3. Path to conditions file (.csv with barcode1_barcode2, condition)
        4. Path to outputfile
'''
from pandas import read_table, DataFrame, read_csv
import csv, argparse, os, itertools, gzip
from itertools import izip

def getParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-folder',
        type=str,
        help='Folder containing read and barcode fastq files')
    parser.add_argument('--ref-file',
        type=str,
        help='Reference file')
    parser.add_argument('--cond-file',
        type=str,
        help='CSV conditions file, dual barcodes should be separated by a _')
    parser.add_argument('--outputfile-path',
        type=str, 
        help='Outputfile')
    return parser

def getCountHash(cond, ref):
    count_hash = {}
    for i,r in enumerate(ref):
        count_hash[r] = count_hash.get(r,{})
        for c in cond:
            count_hash[r][c] = 0
    return count_hash

def getCondHash(bcs_df):
    cond_hash = {}
    for i,b in enumerate(bcs_df['Barcode']):
        cond_hash[b]=bcs_df['Conditions'][i]+'_'+b
    return cond_hash

def getRefHash(ref_df):
    ref_hash = {}
    sp = list(ref_df['Spacer Sequence'])
    for i,r in enumerate(sp):
        ref_hash[r] = ref_df['Spacer'][i]
    return ref_hash

def getOutput(output_hash, cond_hash):
    final_df=DataFrame(output_hash).T.fillna(0)
    colnames = final_df.columns
    new_colnames = [cond_hash[x] for x in colnames]
    final_df.columns = new_colnames
    final_df = final_df[sorted(final_df.columns)]
    return final_df

if __name__ == '__main__':
#    args = getParser().parse_args()
#    folder = args.input_folder
    folder = os.getcwd()
    files = os.listdir(folder)
    ref_file = 'Mitoplus_reference.csv'
    ref_df = read_csv(ref_file)
    ref_df.columns = ['Spacer Sequence', 'Spacer']
    cond_file = 'conditions.csv'
    cond_df = read_csv(cond_file)
    cond_df.columns = ['Barcode', 'Conditions']
    cond_hash = getCondHash(cond_df)
    ref_hash = getRefHash(ref_df)
    count_hash = getCountHash(cond_hash.keys(), ref_hash.keys())
    count=-1
    bc_count=0
    sg_count=0
    accg_count=0
    for f in files:
        if ('L001_I' not in f) and ('fastq' in f):
            print f
            seq_file = folder+'/'+f
            if f[-2:] == 'gz':
                f_open = gzip.open
            else:
                f_open = open
            bc = f.split('_')
            bc[3] = 'I1'
            bc1_name = '_'.join(bc)
            bc1_file = folder+'/'+bc1_name
            bc[3] = 'I2'
            bc2_name = '_'.join(bc)
            bc2_file = folder+'/'+bc2_name
            if bc1_name[-2:] == 'gz':
                bc_open = gzip.open
            else:
                bc_open = open
            with f_open(seq_file, 'r') as s, bc_open(bc1_file, 'r') as b1, bc_open(bc2_file, 'r') as b2:
                for read, r_bar1, r_bar2 in izip(s,b1,b2):
                    count=count+1
                    if count%4 == 1:
                        #print count
                        bc = r_bar1[:-1]+'_'+r_bar2[:-1]
                        if cond_hash.has_key(bc):
                            bc_count+=1
                            pos = read.find('ACCG')
                            if pos != -1:
                                accg_count+=1
                                sgrna = read[pos+4:pos+24]
                                if ref_hash.has_key(sgrna):
                                    sg_count+=1
                                    count_hash[sgrna][bc] += 1
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue

    print '    Total number of reads: '+str(count/4)
    print 'Number of barcode matches: '+str(bc_count)
    print ' Number of vector matches: '+str(accg_count)
    print '  Number of sgrna matches: '+str(sg_count)

    count_df = getOutput(count_hash, cond_hash)
    outputfile = 'MitoPlus_decon_summary.txt'
    count_df.to_csv(outputfile, sep='\t')
