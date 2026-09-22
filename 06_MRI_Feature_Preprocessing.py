import os
import pandas as pd
import re
import yaml
import pip
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
# adding this because neurocombat isn't on anaconda
try:
  import neuroCombat
except:
  pip.main(['install neuroCombat'])
  import neuroCombat

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader = yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']

# pulls out the name and values of global features from the stats output files (same in LH and RH)
def parse_lines(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("# Measure"):
                parts = line.split(", ")
                measure_name = parts[1] 
                value = parts[3]
                data[measure_name] = float(value)
    return data

# processes all the text files 
def process_folder(folder_path, output_csv):
    all_data = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            filepath = os.path.join(folder_path, filename)
            # Use regex to extract Subject, Session, and Hemi information from filename
            subject_match = re.search(r"sub-(\d+)", filename)
            session_match = re.search(r"ses-(\d+)", filename)
            hemi_match = re.search(r"_(lh|rh)\b", filename) 
            
            subject_id = f"sub-{subject_match.group(1)}" if subject_match else None
            session_id = f"ses-{session_match.group(1)}" if session_match else None
            hemi = hemi_match.group(1) if hemi_match else None
            
            # using above function to get stats
            data = parse_lines(filepath)
            data['Subject'] = subject_id
            data['Session'] = session_id
            data['Hemi'] = hemi
            all_data.append(data)
                
    combined_df = pd.DataFrame(all_data)
    combined_df.to_csv(output_csv, index=False)

# uses the two functions above to make the metadata and global feature files
folder_path = os.path.join(scratch_path, "MRI_Processing/DKTstats")
# this doesn't include a mean thickness estimate - global stats need to be derived from aparc (below)
output_csv = os.path.join(scratch_path, "processed_data/DKT_stats_global_fromDKT.csv")
process_folder(folder_path, output_csv)


# 
all_data = []
# Loop through each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)

        # Extract subject, session, and hemisphere information from filename
        subject = re.search(r"sub-(\d+)", filename)
        session = re.search(r"ses-(\d+)", filename)
        hemi = re.search(r"_(lh|rh)\b", filename)
        subject_id = subject.group(1) if subject else "Unknown"
        session_id = session.group(1) if session else "Unknown"
        hemisphere = hemi.group(1) if hemi else "Unknown"

        # Read the file content
        with open(file_path, "r") as file:
            lines = file.readlines()

            in_col_headers_section = False
            
            # Loop through lines to find and parse the table section
            for line in lines:
                # Identify start of table section
                if "ColHeaders" in line:
                    in_col_headers_section = True
                    continue

                if in_col_headers_section:
                    # Split line by whitespace
                    values = re.split(r"\s+", line.strip())
                    
                    # columns check
                    if len(values) == 10:
                        struct_name = values[0]
                        num_vert, surf_area, gray_vol, thick_avg, thick_std, mean_curv, gaus_curv, fold_ind, curv_ind = values[1:]

                        # Append data to list
                        all_data.append({
                            "Subject": subject_id,
                            "Session": session_id,
                            "Hemi": hemisphere,
                            "StructName": struct_name,
                            "NumVert": int(num_vert),
                            "SurfArea": int(surf_area),
                            "GrayVol": int(gray_vol),
                            "ThickAvg": float(thick_avg),
                            "ThickStd": float(thick_std),
                            "MeanCurv": float(mean_curv),
                            "GausCurv": float(gaus_curv),
                            "FoldInd": float(fold_ind),
                            "CurvInd": float(curv_ind)
                        })

DKT_stats = pd.DataFrame(all_data)
output_csv_path = os.path.join(scratch_path, "processed_data/DKT_stats.csv")
DKT_stats.to_csv(output_csv_path, index=False)

# Get global stats from aparc.stats
lh_data = []
rh_data = []
fs_path = os.path.join(scratch_path, "MRI_Processing/fsoutput")
for dirpath, _, filenames in os.walk(fs_path):
    for hemi in ['lh', 'rh']:
        stats_filename = f'{hemi}.aparc.stats'
        if stats_filename in filenames:
            print(f'processing global stats from {dirpath}')
            stats_path = os.path.join(dirpath, stats_filename)

            match = re.search(r'/sub-(\d+)/ses-(\d+)/', stats_path)
            if match:
                Subject = match.group(1)
                Session = match.group(2)

                data = parse_lines(stats_path)
                data.update({'Subject': Subject, 'Session': Session})

                if hemi == 'lh':
                    lh_data.append(data)
                else:
                    rh_data.append(data)

lh_df = pd.DataFrame(lh_data)
rh_df = pd.DataFrame(rh_data)
merged_df = pd.merge(lh_df, rh_df, on=['Subject', 'Session'], suffixes=('_lh', '_rh'))
avg_data = {
    'Subject': merged_df['Subject'],
    'Session': merged_df['Session']
}
values_to_add = ['NumVert', 'WhiteSurfArea']
values_to_average = ['MeanThickness']
same_in_both_files = ['BrainSegVol', 'BrainSegVolNotVent', 'BrainSegVolNotVentSurf', 'CortexVol', 'SupraTentorialVol', 'SupraTentorialVolNotVent', 'eTIV']
for i in values_to_add:
    avg_data[i] = (merged_df[f'{i}_lh'] + merged_df[f'{i}_rh'])
for i in values_to_average:
    avg_data[i] = (merged_df[f'{i}_lh'] + merged_df[f'{i}_rh']) / 2
for i in same_in_both_files:
    avg_data[i] = merged_df[f'{i}_lh']

global_stats_df = pd.DataFrame(avg_data)
global_output_path = os.path.join(scratch_path, "processed_data/DKT_stats_global.csv")
global_stats_df.to_csv(global_output_path, index=False)

# should now have DKT_stats.csv and DKT_stats_global.csv in the processed_data folder

# pull surface holes info
base_dir = os.path.join(scratch_path, "MRI_Processing/fsoutput")
results = []

for sub in os.listdir(base_dir):
    if not sub.startswith('sub-'):
        continue
    projid = sub.split('-')[1]
    sub_path = os.path.join(base_dir, sub)
    ses_dirs = [d for d in os.listdir(sub_path) if d.startswith('ses-')]
    for ses in ses_dirs:
        stats_path = os.path.join(sub_path, ses, 'output', sub, 'stats', 'aseg.stats')
        if not os.path.exists(stats_path):
            continue
        surfaceholes = None
        with open(stats_path, 'r') as f:
            for line in f:
                if line.startswith('# Measure SurfaceHoles'):
                    parts = line.strip().split(',')
                    if len(parts) >= 5:
                        try:
                            surfaceholes = int(parts[3].strip())
                        except ValueError:
                            surfaceholes = None
                    break
        results.append({'projid': projid, 'surfaceholes': surfaceholes, 'session': re.sub('ses-', '', ses)})

surfaceholes_df = pd.DataFrame(results)
holes_csv_path = os.path.join(scratch_path, "processed_data/surfaceholes.csv")
surfaceholes_df.to_csv(holes_csv_path, index=False)

# Some workarounds happen below due to weirdness with the metadata. We need to reconcile the confound dataframe not having 
# session numbers (only dates) with the actual data, which I'm doing by defining which sessions have MRI data using DKT_stats
# then sequentially mapping the metadata onto them

DKT_stats.columns = DKT_stats.columns.str.strip()
DKT_stats = DKT_stats.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
DKT_stats['StructName'] = DKT_stats['Hemi'] + '.' + DKT_stats['StructName']
DKT_stats = DKT_stats.drop(columns=['Hemi'])

transformed_data = DKT_stats[['Subject', 'Session']].drop_duplicates().reset_index(drop=True)
for struct_name in DKT_stats['StructName'].unique():
    struct_pivot = DKT_stats[DKT_stats['StructName'] == struct_name].pivot_table(
        index=['Subject', 'Session'], values=['SurfArea', 'GrayVol', 'ThickAvg', 'CurvInd'], aggfunc='first'
    )
    struct_pivot.columns = [f"{struct_name}.{m}" for m in struct_pivot.columns]
    transformed_data = transformed_data.merge(struct_pivot.reset_index(), on=['Subject', 'Session'], how='left')

transformed_data = transformed_data.loc[:, ~transformed_data.columns.str.contains('temporalpole')] #this is from an earlier dkt version and needs to get dropped
transformed_data['Subject'] = transformed_data['Subject'].astype(int)
transformed_data['Session'] = transformed_data['Session'].astype(int)

global_DKT_stats = pd.read_csv(os.path.join(scratch_path, "processed_data/DKT_stats_global.csv"))
global_DKT_stats['Subject'] = global_DKT_stats['Subject'].astype(int)
global_DKT_stats['Session'] = global_DKT_stats['Session'].astype(int)

combined_data = pd.merge(transformed_data, global_DKT_stats, on=['Subject', 'Session'], how='inner')

for metric in ['NumVert', 'SurfArea', 'GrayVol', 'CortexVol']:
    for col in combined_data.columns:
        if metric in col:
            combined_data[f"{col}_norm"] = combined_data[col] / combined_data['eTIV']
            combined_data = combined_data.drop(columns=[col])
combined_data = combined_data.drop(columns=['eTIV'])
combined_data.to_csv(os.path.join(scratch_path, 'processed_data/DKT_stats_normalized.csv'), index=False)

# this one is where the actual scanner info comes from

metadata_path = os.path.join(data_path, 'ROSMAP_MRI/mri_metadata.csv')
metadata_df = pd.read_csv(metadata_path)
metadata_df['Subject'] = metadata_df['Subject'].str.replace('sub-', '', regex=False).astype(int)
metadata_df['Session'] = metadata_df['Session'].str.replace('ses-', '', regex=False).astype(int)

batcheffect = 'DeviceSerialNumber'
batch_info = metadata_df[['Subject', 'Session', batcheffect]].drop_duplicates(subset=['Subject', 'Session'])

combat_input = pd.merge(combined_data, batch_info, on=['Subject', 'Session'], how='inner')
combat_input = combat_input.dropna(subset=[batcheffect])

dkt_pairs = combined_data[['Subject', 'Session']].drop_duplicates()
matched_pairs = combat_input[['Subject', 'Session']].drop_duplicates()
unmatched_batch = pd.merge(dkt_pairs, matched_pairs, on=['Subject', 'Session'], how='left', indicator=True)
unmatched_batch = unmatched_batch.query('_merge == "left_only"').drop(columns='_merge')

print("Batch effect (metadata_df) matching")
print(f"DKT (Subject, Session) pairs: {len(dkt_pairs)}")
print(f"Matched with {batcheffect}: {len(matched_pairs)}")
print(f"Unmatched: {len(unmatched_batch)}")
if not unmatched_batch.empty:
    print(unmatched_batch.to_string(index=False))
print()

# matching with the age_scan data - need this later so we know age at scan (rosmaster stops at 90 for some reason)

imaging_metadata_path = os.path.join(data_path, 'ROSMAP_MRI/rosmapAD_data_imaging_confounds_ageScan.csv')
mri_info_df = pd.read_csv(imaging_metadata_path)
mri_info_df.rename(columns={mri_info_df.columns[0]: 'Subject', mri_info_df.columns[3]: 'Scanner'}, inplace=True)
mri_info_df['Subject'] = mri_info_df['Subject'].astype(str)
mri_info_df['mri_date'] = mri_info_df['mri_date'].astype(str)
mri_info_df['mri_date'] = mri_info_df['mri_date'].apply(lambda x: '0' + x if len(x) == 5 else x)
bad_dates = mri_info_df.loc[mri_info_df['mri_date'].str.len() != 6, 'mri_date']
assert bad_dates.empty, f"Found mri_date values that aren't 5 or 6 digits after padding:\n{bad_dates.unique()}"
mri_info_df['mri_date'] = pd.to_datetime(mri_info_df['mri_date'], format='%y%m%d')
mri_info_df = mri_info_df.sort_values(['Subject', 'mri_date']).reset_index(drop=True)

sessions_with_mri = DKT_stats[['Subject', 'Session']].drop_duplicates().sort_values(['Subject', 'Session'])
sessions_with_mri['Subject'] = sessions_with_mri['Subject'].astype(str)

assigned_groups = []
count_mismatches = []
missing_from_info = []

dkt_subjects = set(sessions_with_mri['Subject'])
info_subjects = set(mri_info_df['Subject'])

for subject in dkt_subjects - info_subjects:
    n_sessions = len(sessions_with_mri.loc[sessions_with_mri['Subject'] == subject])
    missing_from_info.append({'Subject': subject, 'n_dkt_sessions': n_sessions})

for subject, group in mri_info_df.groupby('Subject', group_keys=False):
    known_sessions = sessions_with_mri.loc[sessions_with_mri['Subject'] == subject, 'Session'].tolist()
    if not known_sessions:
        continue
    if len(group) != len(known_sessions):
        count_mismatches.append({
            'Subject': subject, 'n_scans_in_info': len(group),
            'n_dkt_sessions': len(known_sessions), 'dkt_sessions': known_sessions
        })
        n = min(len(group), len(known_sessions))
        group = group.copy().sort_values('mri_date').iloc[:n]
        group['Session'] = known_sessions[:n]
        assigned_groups.append(group)
        continue
    group = group.copy()
    group['Session'] = known_sessions
    assigned_groups.append(group)

mri_info_assigned = pd.concat(assigned_groups, ignore_index=True) if assigned_groups else mri_info_df.iloc[0:0]
mri_info_assigned['Subject'] = mri_info_assigned['Subject'].astype(int)
mri_info_assigned['Session'] = mri_info_assigned['Session'].astype(int)

age_info = mri_info_assigned[['Subject', 'Session', 'age_scan', 'mri_date']].drop_duplicates(subset=['Subject', 'Session'])
matched_age_pairs = age_info[['Subject', 'Session']].drop_duplicates()
unmatched_age = pd.merge(dkt_pairs, matched_age_pairs, on=['Subject', 'Session'], how='left', indicator=True)
unmatched_age = unmatched_age.query('_merge == "left_only"').drop(columns='_merge')

age_info_path = os.path.join(scratch_path, "processed_data/MRI_age_info.csv")
age_info.to_csv(age_info_path, index=False)

print("Age info (mri_info_df) matching")
print(f"DKT (Subject, Session) pairs: {len(dkt_pairs)}")
print(f"Matched with age: {len(matched_age_pairs)}")
print(f"Unmatched: {len(unmatched_age)}\n")

# Drop ones with a lot of holes (assuming that's due to lots of movement)

SURFACEHOLES_THRESHOLD = 250 

surfaceholes_df['Subject'] = surfaceholes_df['projid'].astype(int)
surfaceholes_df['Session'] = surfaceholes_df['session'].astype(int)

combat_input = pd.merge(
    combat_input, surfaceholes_df[['Subject', 'Session', 'surfaceholes']],
    on=['Subject', 'Session'], how='left'
)
n_before_qc = len(combat_input)
combat_input = (combat_input[combat_input['surfaceholes'] <= SURFACEHOLES_THRESHOLD].drop(columns=['surfaceholes']).reset_index(drop=True))
print(f"Surface holes QC: dropped {n_before_qc - len(combat_input)} / {n_before_qc} scans (threshold={SURFACEHOLES_THRESHOLD})\n")
assert not combat_input.duplicated(subset=['Subject','Session']).any()

# make PCA plot before combat
combat_input = combat_input.drop(['BrainSegVol','BrainSegVolNotVent','BrainSegVolNotVentSurf','SupraTentorialVol','SupraTentorialVolNotVent'], axis ='columns')
precombat_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_before_combat.csv")
combat_input.to_csv(precombat_data_path, index=False)
meta_cols = ['Subject', 'Session', batcheffect]
feature_cols = [c for c in combat_input.columns if c not in meta_cols]
numeric_features = combat_input[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

scaled_before = StandardScaler().fit_transform(numeric_features)
pca_before = PCA(n_components=10).fit_transform(scaled_before)

plot_before = combat_input[[batcheffect]].copy()
plot_before['PC1'], plot_before['PC2'] = pca_before[:, 0], pca_before[:, 1]

fig = plt.figure(figsize=(10, 8))
sns.scatterplot(x='PC1', y='PC2', hue=batcheffect, data=plot_before)
plt.title('PCA of MRI Features Before Batch Correction')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title=batcheffect, loc='best')
with PdfPages(os.path.join(scratch_path, 'figures/supplemental/06_MRI_PCA_before.pdf')) as pdf:
    pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)



# Combat
covars = combat_input[[batcheffect]]
data_combat = neuroCombat.neuroCombat(dat=scaled_before.T, covars=covars, batch_col=batcheffect)["data"].T
combat_df = pd.DataFrame(data_combat, columns=feature_cols).apply(pd.to_numeric, errors='coerce')
n_nan = combat_df.isna().sum().sum()
n_inf = np.isinf(data_combat).sum()
print("NA and infinite value count:")
print(n_nan, n_inf)
combat_df = combat_df.fillna(combat_df.mean())

# PCA after

pca_after = PCA(n_components=10).fit_transform(combat_df)

plot_after = combat_input[[batcheffect]].copy()
plot_after['PC1'], plot_after['PC2'] = pca_after[:, 0], pca_after[:, 1]

fig = plt.figure(figsize=(10, 8))
sns.scatterplot(x='PC1', y='PC2', hue=batcheffect, data=plot_after)
plt.title(f'PCA of MRI Features After {batcheffect} Removal')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title=batcheffect, loc='best')
with PdfPages(os.path.join(scratch_path, 'figures/supplemental/06_MRI_PCA_after.pdf')) as pdf:
    pdf.savefig(fig, bbox_inches='tight')
plt.close(fig)

for i in range(10):
    combat_input[f'PC{i+1}'] = pca_after[:, i]


combat_df = pd.concat([combat_input[['Subject', 'Session']], combat_df], axis=1)
combat_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_combat.csv")
combat_df.to_csv(combat_data_path, index=False)

# Identify the columns to preserve (Subject and Session)
preserve_columns = ['Subject', 'Session']

# Filter lh and rh columns
lh_columns = [col for col in combat_df.columns if col.startswith('lh.')]
rh_columns = [col for col in combat_df.columns if col.startswith('rh.')]

# Create separate DataFrames for lh and rh
lh_df = combat_df[preserve_columns + lh_columns]
rh_df = combat_df[preserve_columns + rh_columns]
non_lr_columns = [c for c in combat_df.columns if not c.startswith('lh.') and not c.startswith('rh.') and c not in preserve_columns]
global_df = combat_df[preserve_columns + non_lr_columns]

# Save to separate CSV files
lh_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_combat_lh.csv")
lh_df.to_csv(lh_data_path, index=False)
rh_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_combat_rh.csv")
rh_df.to_csv(rh_data_path, index=False)
global_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_combat_global.csv")
global_df.to_csv(global_data_path, index=False)

# Create a new DataFrame to store the averaged values
averaged_df = combat_df[preserve_columns].copy()

# Iterate over each unique region/feature pair
for lh_col, rh_col in zip(lh_columns, rh_columns):
    # Extract the region/feature name by removing the 'lh.' or 'rh.' prefix
    region_feature = lh_col.split('.', 1)[1]
    
    # Compute the average between lh and rh
    averaged_df[region_feature] = (lh_df[lh_col] + rh_df[rh_col]) / 2

# Save the averaged values to a new CSV file
averaged_data_path = os.path.join(scratch_path, "processed_data/DKT_stats_combat_averaged.csv")
averaged_df.to_csv(averaged_data_path, index=False)

total_scans = len(averaged_df)
total_individuals = averaged_df['Subject'].nunique()

print(f'processed {total_scans} scans from {total_individuals} individuals')