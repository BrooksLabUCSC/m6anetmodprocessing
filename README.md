# m6anetmodprocessing
Scripts to process m6anet output

To get the number of isoforms/genes in m6anet output:

```bash
python3 count_transcripts_and_genes.py data.site_proba.csv
```

Example output:
```bash
number of isoforms with mod information:  18103
number of genes with mod information:  7189
average number of isoforms per gene with mod information:  2.52
```

How to run Selam's mod ratio making table:
```bash
python3 mod_ratio_table.py -i manifest.tsv -o output.csv
```

selam script sample manifest format:
```
sample_name	file_path
mt1dmso		/path/to/data.site_proba.csv
mt1csc		/path/to/data.site_proba.csv
...
```
