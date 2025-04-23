import sys

infile = sys.argv[1] ###assume to be data.site_proba.csv

isoforms, genes = set(), set()
for line in open(infile):
    isoid = line.split(',')[0]
    gene = isoid.split('_')[-1]
    isoforms.add(isoid)
    genes.add(gene)

print('number of isoforms with mod information: ', len(isoforms))
print('number of genes with mod information: ', len(genes))
print('average number of isoforms per gene with mod information: ', round(len(isoforms)/len(genes), 2))

