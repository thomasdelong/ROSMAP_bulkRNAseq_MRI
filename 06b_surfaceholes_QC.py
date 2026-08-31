import os
import glob
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import imageio

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader=yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']

surfaceholes_path = os.path.join(scratch_path, "processed_data/surfaceholes.csv")
surfaceholesdf = pd.read_csv(surfaceholes_path)

if 'surfaceholes' not in surfaceholesdf.columns:
    raise ValueError("Column 'surfaceholes' not found in the CSV file.")

# histogram of surfaceholes in all scans
plt.figure(figsize=(8, 6))
plt.hist(surfaceholesdf['surfaceholes'].dropna(), bins=20, color='skyblue', edgecolor='black')
plt.title('Histogram of Surface Holes (All Scans)')
plt.xlabel('Number of Surface Holes')
plt.ylabel('Frequency')
plt.grid(True)
plt.savefig(os.path.join(scratch_path, 'figures/supplemental/06b_surfaceholes_all.pdf'), bbox_inches='tight')
plt.show()

# get just the final scans (ones used in analyses)
highest_session_df = (
    surfaceholesdf.sort_values(['projid', 'session'])
    .groupby('projid', as_index=False)
    .last()
)

# find intersection with RNA data
processed_data_dir = os.path.join(scratch_path, "processed_data")
rnaseq_files = glob.glob(os.path.join(processed_data_dir, "*RNAseq_metadata.csv"))

rnaseq_projids = set()
for f in rnaseq_files:
    df = pd.read_csv(f)
    rnaseq_projids.update(df['projid'].astype(str))

highest_session_df['projid'] = highest_session_df['projid'].astype(str)
rnaseq_subset = highest_session_df[highest_session_df['projid'].isin(rnaseq_projids)]

print(f"Total MRI subjects: {len(highest_session_df)}")
print(f"Intersected with RNAseq projids: {len(rnaseq_subset)}")

# histo for just the scans used in the paper
plt.figure(figsize=(8, 6))
plt.hist(rnaseq_subset['surfaceholes'].dropna(), bins=20, color='skyblue', edgecolor='black')
plt.title('Histogram of Surface Holes (Used Scans Only)')
plt.xlabel('Number of Surface Holes')
plt.ylabel('Frequency')
plt.grid(True)
plt.savefig(os.path.join(scratch_path, 'figures/supplemental/06b_surfaceholes_rnaseq_highest_session.pdf'), bbox_inches='tight')
plt.show()

# gifs of the gray matter for the 3 scans with the most and least holes

def create_gray_matter_gif(folder_path, output_gif='gray_matter.gif', save_path=None):
    brain_path = os.path.join(folder_path, 'brain.mgz')
    aseg_path = os.path.join(folder_path, 'aseg.mgz')
    brain_img = nib.load(brain_path)
    brain_data = brain_img.get_fdata()

    aseg_img = nib.load(aseg_path)
    aseg_data = aseg_img.get_fdata()

    cerebral_gm_labels = [3, 42]
    gray_matter_mask = np.isin(aseg_data, cerebral_gm_labels)
    frames = []

    for slice_index in range(brain_data.shape[2]):
        brain_slice = brain_data[:, :, slice_index]
        mask_slice = gray_matter_mask[:, :, slice_index]
        if np.max(brain_slice) == 0:
            continue
        brain_rot = np.rot90(brain_slice, k=3)
        mask_rot = np.rot90(mask_slice, k=3)

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(brain_rot, cmap='gray', alpha=0.6)
        ax.imshow(mask_rot, cmap='Blues', alpha=0.5)
        ax.set_title(f'Cerebral Gray Matter (Slice {slice_index})')
        ax.axis('off')

        fig.canvas.draw()
        image = np.asarray(fig.canvas.buffer_rgba())
        image = image[:, :, :3]
        frames.append(image)
        plt.close(fig)

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        gif_path = save_path
    else:
        gif_path = output_gif
    imageio.mimsave(gif_path, frames, duration=0.1, loop=0)
    print(f"Saved GIF to: {gif_path}")

highest_3 = surfaceholesdf.nlargest(3, 'surfaceholes')
lowest_3 = surfaceholesdf.nsmallest(3, 'surfaceholes')

gif_output_dir = os.path.join(scratch_path, 'figures/supplemental/06b_graymatter_gifs')
os.makedirs(gif_output_dir, exist_ok=True)

fsoutput_dir = os.path.join(scratch_path, "MRI_Processing/fsoutput")

def make_gif_for_row(row, label):
    subject_id = f"sub-{row['projid']}"
    session_id = f"ses-{row['session']}"
    mri_dir = os.path.join(fsoutput_dir, subject_id, session_id, "output", subject_id, "mri")
    save_path = os.path.join(gif_output_dir, f"{label}_{subject_id}_{session_id}.gif")
    if not os.path.exists(os.path.join(mri_dir, 'brain.mgz')):
        print(f"WARNING: brain.mgz not found for {subject_id} {session_id}, skipping.")
        return None
    create_gray_matter_gif(mri_dir, save_path=save_path)
    return save_path
    
highest_paths = []
for i, (_, row) in enumerate(highest_3.iterrows()):
    highest_paths.append(make_gif_for_row(row, f"highest{i+1}"))

lowest_paths = []
for i, (_, row) in enumerate(lowest_3.iterrows()):
    lowest_paths.append(make_gif_for_row(row, f"lowest{i+1}"))

# middle frames of the gifs to a figure for the supplemental (do any journals take gifs?)

def get_middle_frame(gif_path):
    if gif_path is None or not os.path.exists(gif_path):
        return None
    frames = imageio.mimread(gif_path)
    if not frames:
        return None
    return frames[len(frames) // 2]

fig, axes = plt.subplots(3, 2, figsize=(8, 12))

for row_idx in range(3):
    low_frame = get_middle_frame(lowest_paths[row_idx])
    ax = axes[row_idx, 0]
    if low_frame is not None:
        ax.imshow(low_frame)
        ax.set_title(f"Lowest {row_idx+1}: {lowest_3.iloc[row_idx]['projid']} "
                     f"({int(lowest_3.iloc[row_idx]['surfaceholes'])} holes)", fontsize=9)
    else:
        ax.text(0.5, 0.5, "N/A", ha='center', va='center')
    ax.axis('off')

    high_frame = get_middle_frame(highest_paths[row_idx])
    ax = axes[row_idx, 1]
    if high_frame is not None:
        ax.imshow(high_frame)
        ax.set_title(f"Highest {row_idx+1}: {highest_3.iloc[row_idx]['projid']} "
                     f"({int(highest_3.iloc[row_idx]['surfaceholes'])} holes)", fontsize=9)
    else:
        ax.text(0.5, 0.5, "N/A", ha='center', va='center')
    ax.axis('off')

fig.suptitle('Gray Matter QC: Lowest vs Highest Surface Hole Counts (Middle Slice)', fontsize=12)
plt.tight_layout()

comparison_path = os.path.join(scratch_path, 'figures/supplemental/06b_graymatter_comparison.pdf')
os.makedirs(os.path.dirname(comparison_path), exist_ok=True)
plt.savefig(comparison_path, bbox_inches='tight')
plt.show()
print(f"Saved comparison figure to: {comparison_path}")