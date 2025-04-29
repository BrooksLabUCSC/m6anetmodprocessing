import pandas as pd 
import numpy as np 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--output','-o', type=str, action='store', help = 'output file')
parser.add_argument('--input','-i', type=str, action='store', help = 'input manifest')
args = parser.parse_args()
output = args.output
manifest = args.input
print(output)
sample_manifest = pd.read_csv(manifest, sep='\t')


all_samples = []
sample_names = []

for i, row in sample_manifest.iterrows():
    sample_name = row['sample_name']
    file_path = row['file_path']
    df = pd.read_csv(file_path)
    # print(df)
    # df[['transcript_id', 'gene_id']] = df['transcript_id'].str.split('_', expand=True)
    split_cols = df['transcript_id'].str.split('_', expand=True)
    if split_cols.shape[1] >= 2:
        df['transcript_id'] = split_cols[0]
        df['gene_id'] = split_cols[1]

    df = df[['transcript_id', 'gene_id', 'transcript_position', 'n_reads', 'probability_modified', 'kmer', 'mod_ratio']]
    all_samples.append(df)
    sample_names.append(sample_name)

shared_transcripts = set(zip(all_samples[0]['transcript_id'], all_samples[0]['transcript_position']))
for df in all_samples[1:]:
    # Get the set of (transcript_id, transcript_position) pairs from the current dataframe
    current_transcripts = set(zip(df['transcript_id'], df['transcript_position']))

    # Keep only the common pairs (intersection) between the current set and the shared set
    shared_transcripts &= current_transcripts

# Now, shared_transcripts contains the common (transcript_id, transcript_position) pairs
# print(len(shared_transcripts))
# total = 4951
print(len(shared_transcripts))

filtered_dfs = []
for df in all_samples:
    df_filtered = df[df.set_index(['transcript_id', 'transcript_position']).index.isin(shared_transcripts)]
    filtered_dfs.append(df_filtered)


filtered_dfs_3plus = []
for df in filtered_dfs:
    # Count the number of unique transcript_positions for each transcript_id
    transcripts_with_3plus = df.groupby('transcript_id')['transcript_position'].nunique()
    # Keep only those transcripts with 3 or more unique mod positions
    df_filtered_3plus = df[df['transcript_id'].isin(transcripts_with_3plus[transcripts_with_3plus >= 3].index)]
    filtered_dfs_3plus.append(df_filtered_3plus)


filtered_dfs_final = []
for df in filtered_dfs_3plus:
    transcripts_per_gene = df.groupby('gene_id')['transcript_id'].nunique()
    # Keep only those genes with 2 or more transcripts
    df_filtered_gene = df[df['gene_id'].isin(transcripts_per_gene[transcripts_per_gene >= 2].index)]
    filtered_dfs_final.append(df_filtered_gene)
filtered_dfs_final[0]

combined_df = []
for df, sample_name in zip(filtered_dfs_final, ['mt1dmso', 'mt1csc', 'mt2csc', 'mt2dmso', 'wt1csc', 'wt1dmso', 'wt2csc', 'wt2dmso']):
    df['transcript_pos'] = df['gene_id'] + "_" + df['transcript_id'] + "_" + df['transcript_position'].astype(str)
    df_sample = df[['transcript_pos', 'mod_ratio']]
    df_sample.rename(columns={'mod_ratio': sample_name}, inplace=True)
    combined_df.append(df_sample)

final_combined_df = combined_df[0]
for df in combined_df[1:]:
    final_combined_df = pd.merge(final_combined_df, df, on='transcript_pos', how='outer')

final_combined_df.reset_index(drop=True, inplace=True)
# now i want to save to my computer as a csv file
final_combined_df.to_csv(output, index=False)
