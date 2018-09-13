"""merge the 42 clinical files downloaded from gdc, adding a column maf_type.

Usage: PROG tsvfiles > outcsvfile
"""
import sys, pandas as pd

def main(infiles, outfile):
    mafs = []
    for f in infiles:
        maf_group = f.split('.')[1]
        print (maf_group, file=sys.stderr)
        curr = pd.read_csv(f, delimiter='\t')
        curr['maf_group'] = maf_group
        mafs.append(curr)
    res = pd.concat(mafs)
    res.to_csv(outfile, index=False)

if __name__ == '__main__':
    res = main(sys.argv[1:], sys.stdout)
