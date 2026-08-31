# takes all the previous count matrices and performs PCA on the data before and after batch correction
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import yaml

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader = yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']

def pcafigure(final_df, tissue_name_lower):
    fig = plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', hue='sequencingBatch', data=final_df, palette='gist_rainbow')
    plt.title(f'PCA of {tissue_name_lower} gene counts')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend(title='Sequencing Batch', loc='best')
    return fig

def run_pca_for_tissue(tissue_abbreviation, tissue_name, projidmapping_main, counts_suffix):
    rnaseq_df = pd.read_csv(
        os.path.join(scratch_path, 'processed_data', f'{tissue_abbreviation}_RNAseq_counts{counts_suffix}.csv'),
        index_col=0
    )
    metadata_df = pd.read_csv(
        os.path.join(scratch_path, 'processed_data', f'{tissue_abbreviation}_RNAseq_metadata.csv')
    )
    projidmapping = projidmapping_main[
        (projidmapping_main['tissue'] == tissue_name) & (projidmapping_main['assay'] == 'rnaSeq')
    ]

    rnaseq_transposed = rnaseq_df.T
    rnaseq_transposed['projid'] = rnaseq_transposed.index
    rnaseq_transposed = rnaseq_transposed[rnaseq_transposed['projid'] != 'NA']
    rnaseq_transposed['projid'] = rnaseq_transposed['projid'].map(projidmapping.set_index('specimenID')['projid'])

    rnaseq_transposed['projid'] = rnaseq_transposed['projid'].astype(str).str.split('.').str[0]
    metadata_df['projid'] = metadata_df['projid'].astype(str).str.split('.').str[0]

    metadata_df = metadata_df[metadata_df['projid'].notna()]
    metadata_df['projid'] = metadata_df['projid'].apply(lambda x: x.zfill(8))

    rnaseq_transposed = rnaseq_transposed.drop_duplicates(subset='projid', keep='first')
    metadata_df = metadata_df.drop_duplicates(subset='projid', keep='first')

    numeric_columns = rnaseq_transposed.drop(columns=['projid'])

    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(numeric_columns)

    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(10)])
    pca_df['projid'] = rnaseq_transposed['projid'].values

    final_df = pd.merge(pca_df, metadata_df, on='projid', how='left').sort_values(by='sequencingBatch')

    tissue_name_lower = tissue_name.lower()
    return pcafigure(final_df, tissue_name_lower)

projidmapping_main = pd.read_csv(
    os.path.join(data_path, 'ROSMAP_RNA/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv')
)

tissue_list = {
    'dlpfc': 'dorsolateral prefrontal cortex',
    'tc': 'temporal cortex',
    'fc': 'frontal cortex',
    'pcc': 'posterior cingulate cortex',
    'hcn': 'Head of caudate nucleus'
}

output_dir = os.path.join(scratch_path, 'figures/supplemental')
os.makedirs(output_dir, exist_ok=True)

# Before batch correction — one PDF per tissue
for tissue_abbreviation, tissue_name in tissue_list.items():
    fig = run_pca_for_tissue(tissue_abbreviation, tissue_name, projidmapping_main, counts_suffix='')
    out_path = os.path.join(output_dir, f'03_PCA_before_{tissue_abbreviation}.pdf')
    with PdfPages(out_path) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# After batch correction — one PDF per tissue
for tissue_abbreviation, tissue_name in tissue_list.items():
    fig = run_pca_for_tissue(tissue_abbreviation, tissue_name, projidmapping_main, counts_suffix='_filtered_corrected')
    out_path = os.path.join(output_dir, f'03_PCA_after_{tissue_abbreviation}.pdf')
    with PdfPages(out_path) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)