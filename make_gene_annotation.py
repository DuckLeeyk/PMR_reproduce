import pandas as pd

pd.set_option('display.max_columns', None)

col_names = [
    'chrom',      # Chromosome
    'source',     # Annotation source, e.g., HAVANA
    'feature',    # Feature type, e.g., gene, exon
    'start',      # Start position
    'end',        # End position
    'score',      # Score, usually '.'
    'strand',     # Strand (+/-)
    'frame',      # Reading frame, usually '.'
    'attributes'  # Last column, like key "value"; key2 "value2";
]

df = pd.read_csv(
    'data/Gex_summary/gencode.v12.annotation.gtf.gz',
    sep='\t',
    comment='#',        # Ignore comment lines
    header=None,
    names=col_names,
    compression='gzip',
    dtype={'chrom': str}
)

# 2. Parse the attributes column
def parse_attributes(attr_str):
    """
    Parse 'key1 "val1"; key2 "val2"; …' into a dictionary
    """
    d = {}
    # Split by semicolon, then extract key/value pairs
    for field in attr_str.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        key, val = field.split(' ', 1)
        d[key] = val.strip().strip('"')
    return d

# Convert the attributes string in each row to a dict, then expand into multiple columns
attr_df = df['attributes']\
    .apply(parse_attributes)\
    .apply(pd.Series)

# 3. Merge the parsed attribute columns back into the original DataFrame
df = pd.concat([df.drop(columns=['attributes']), attr_df], axis=1)
df = df[['chrom', 'feature', 'start', 'end', 'gene_id', 'gene_type', 'gene_name']]

df = df[~df['chrom'].isin(['chrM', 'chrX', 'chrY'])].reset_index(drop=True)
df = df[df['feature'] == 'gene'].reset_index(drop=True)
df = df[df['gene_type'].isin(['protein_coding', 'lincRNA'])].reset_index(drop=True)

df = df.rename(columns={'chrom': 'chr'})
df['chr'] = df['chr'].str.replace(r'^chr', '', regex=True)
df.drop(columns=['feature'], inplace=True)

df.to_csv("data/Gex_summary/gencode_v12_gene_annotation.csv", index=False)
