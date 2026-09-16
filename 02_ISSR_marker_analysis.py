#!/usr/bin/env python3
"""
=============================================================================
SCRIPT 02: ISSR Molecular Marker Analysis
Author:    gh (Postgraduate Researcher, Bioinformatics & Molecular Biology)
Data:      Binary gel scoring data from ISSR PCR bands
Method:    Jaccard & Dice similarity, Genetic distance, UPGMA clustering
Output:    Similarity matrices, distance matrix, Newick tree, summary stats
=============================================================================

DESCRIPTION:
    This script processes binary band-scoring data from ISSR (Inter-Simple
    Sequence Repeat) gel electrophoresis. It:
    1. Reads a binary presence/absence matrix (1 = band, 0 = no band)
    2. Combines replicate rows per sample using logical OR
    3. Calculates Jaccard and Dice similarity coefficients
    4. Calculates genetic distance (1 - Jaccard)
    5. Performs UPGMA hierarchical clustering
    6. Outputs a Newick-format tree for use in MEGA or FigTree
    7. Calculates band frequency and PIC (Polymorphism Information Content)
    8. Exports all results to CSV files

USAGE:
    python3 02_ISSR_marker_analysis.py

    To use your own data, replace the RAW_DATA dictionary below with your
    binary scoring table. Format: {sample_name: [[row1_scores], [row2_scores]]}

DEPENDENCIES:
    numpy, pandas, scipy (install with: pip install numpy pandas scipy)
=============================================================================
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving figures
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# 1. INPUT DATA
# =============================================================================
# Format: {sample_name: [[row1_band_scores], [row2_band_scores]]}
# Each row = one gel lane (replicate). 1 = band present, 0 = band absent.
# Replace this with your own data.

# Loci 1-4 dataset
RAW_DATA_L1_4 = {
    'A': [[1, 1, 1, 1], [0, 0, 0, 0]],
    'B': [[1, 1, 1, 1], [1, 1, 0, 0]],
    'C': [[1, 1, 1, 1], [0, 1, 0, 0]],
    'D': [[1, 1, 1, 1], [1, 1, 0, 0]],
    'E': [[1, 1, 1, 1], [0, 0, 0, 0]],
    'F': [[1, 1, 1, 1], [0, 1, 0, 0]],
    'G': [[1, 0, 1, 1], [1, 0, 0, 0]],
    'H': [[1, 1, 1, 1], [1, 1, 0, 0]],
    'I': [[1, 0, 1, 1], [1, 0, 0, 0]],
}
LOCUS_NAMES_L1_4 = ['Locus1', 'Locus2', 'Locus3', 'Locus4']

# Loci 5-8 dataset
RAW_DATA_L5_8 = {
    'A': [[0, 0, 0, 1], [1, 1, 1, 1]],
    'B': [[0, 0, 1, 1], [1, 0, 1, 1]],
    'C': [[0, 1, 1, 1], [1, 1, 1, 1]],
    'D': [[0, 0, 1, 1], [0, 0, 0, 1]],
    'E': [[1, 1, 1, 1], [1, 1, 0, 1]],
    'F': [[1, 1, 1, 0], [1, 1, 0, 0]],
    'G': [[1, 0, 0, 1], [1, 0, 0, 1]],
    'H': [[1, 0, 1, 1], [1, 0, 0, 0]],
    'I': [[0, 0, 1, 1], [1, 0, 0, 0]],
}
LOCUS_NAMES_L5_8 = ['Locus5', 'Locus6', 'Locus7', 'Locus8']


# =============================================================================
# 2. CORE FUNCTIONS
# =============================================================================

def combine_replicates(raw_data):
    """
    Combine replicate rows per sample using logical OR.
    A band is scored as present (1) if it appears in either replicate row.

    Parameters:
        raw_data (dict): {sample: [[row1], [row2]]}

    Returns:
        dict: {sample: [combined_scores]}
    """
    combined = {}
    for sample, rows in raw_data.items():
        combined[sample] = [max(rows[0][i], rows[1][i]) for i in range(len(rows[0]))]
    return combined


def jaccard_similarity(a, b):
    """
    Calculate Jaccard similarity between two binary vectors.
    Jaccard = a / (a + b + c)
    where:
        a = bands present in BOTH samples
        b = bands present in sample 1 ONLY
        c = bands present in sample 2 ONLY

    Parameters:
        a, b (list or array): binary vectors

    Returns:
        float: Jaccard similarity coefficient (0 to 1)
    """
    a, b = np.array(a), np.array(b)
    shared  = np.sum((a == 1) & (b == 1))
    either  = np.sum((a == 1) | (b == 1))
    return round(float(shared / either), 4) if either > 0 else 0.0


def dice_similarity(a, b):
    """
    Calculate Dice similarity coefficient between two binary vectors.
    Dice = 2a / (2a + b + c)
    Dice weights shared bands more heavily than Jaccard.
    Preferred for dominant markers (ISSR, RAPD, AFLP).

    Parameters:
        a, b (list or array): binary vectors

    Returns:
        float: Dice similarity coefficient (0 to 1)
    """
    a, b = np.array(a), np.array(b)
    shared = np.sum((a == 1) & (b == 1))
    total  = np.sum(a == 1) + np.sum(b == 1)
    return round(float(2 * shared / total), 4) if total > 0 else 0.0


def build_similarity_matrix(combined, sim_func):
    """
    Build a pairwise similarity matrix for all samples.

    Parameters:
        combined (dict): {sample: [combined_scores]}
        sim_func (function): similarity function (jaccard or dice)

    Returns:
        tuple: (matrix as 2D list, sample_names as list)
    """
    samples = list(combined.keys())
    n       = len(samples)
    matrix  = [[0.0] * n for _ in range(n)]

    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            matrix[i][j] = sim_func(combined[s1], combined[s2])

    return matrix, samples


def band_frequency_and_pic(combined, locus_names):
    """
    Calculate band frequency and PIC for each locus.

    Band frequency = proportion of samples showing a band at this locus.
    PIC (Polymorphism Information Content) = 1 - sum(p_i^2)
    For a dominant binary marker: PIC = 1 - (p^2 + q^2)
    where p = band frequency and q = 1 - p.
    PIC = 0 means the locus is monomorphic (no information).
    PIC > 0.5 indicates high informativeness.

    Parameters:
        combined (dict): {sample: [combined_scores]}
        locus_names (list): names of loci

    Returns:
        pd.DataFrame: summary statistics per locus
    """
    samples  = list(combined.keys())
    n        = len(samples)
    freq_lst = []
    pic_lst  = []

    for li in range(len(locus_names)):
        p = sum(combined[s][li] for s in samples) / n
        q = 1 - p
        pic = round(1 - (p**2 + q**2), 4)
        freq_lst.append(round(p, 4))
        pic_lst.append(pic)

    df = pd.DataFrame({
        'Locus'          : locus_names,
        'Band_Frequency' : freq_lst,
        'PIC'            : pic_lst,
        'Polymorphic'    : ['Yes' if p > 0 else 'No' for p in pic_lst]
    })
    return df


def upgma_newick(sample_names, dist_matrix):
    """
    Perform UPGMA hierarchical clustering and return a Newick format tree.
    Uses scipy's average linkage (UPGMA) on a condensed distance matrix.

    Parameters:
        sample_names (list): list of sample names
        dist_matrix (list of lists): pairwise distance matrix

    Returns:
        str: Newick format tree string
    """
    dist_array = np.array(dist_matrix)
    # Convert square matrix to condensed form for scipy
    condensed  = squareform(dist_array, checks=False)
    # UPGMA = average linkage in scipy
    linkage_matrix = linkage(condensed, method='average')

    # Convert linkage matrix to Newick format
    def get_newick(node, parent_dist, leaf_names, newick=''):
        if node.is_leaf():
            return f"{leaf_names[node.id]}:{round(parent_dist - node.dist, 4)}"
        else:
            left  = get_newick(node.get_left(),  node.dist, leaf_names)
            right = get_newick(node.get_right(), node.dist, leaf_names)
            return f"({left},{right}):{round(parent_dist - node.dist, 4)}"

    tree = to_tree(linkage_matrix, rd=False)
    newick_str = get_newick(tree, tree.dist, sample_names) + ";"
    return newick_str


def plot_dendrogram(sample_names, dist_matrix, title, filename):
    """
    Plot and save a UPGMA dendrogram.

    Parameters:
        sample_names (list): list of sample names
        dist_matrix (list of lists): pairwise distance matrix
        title (str): plot title
        filename (str): output filename (.png)
    """
    dist_array = np.array(dist_matrix)
    condensed  = squareform(dist_array, checks=False)
    linkage_matrix = linkage(condensed, method='average')

    fig, ax = plt.subplots(figsize=(10, 5))
    dendrogram(
        linkage_matrix,
        labels      = sample_names,
        orientation = 'top',
        leaf_rotation = 0,
        ax          = ax,
        color_threshold = 0.3,
        above_threshold_color = '#2C6B2F',
    )
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel('Sample', fontsize=11)
    ax.set_ylabel('Genetic Distance (1 - Jaccard)', fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Dendrogram saved: {filename}")


def export_matrix_csv(matrix, sample_names, filename, label="Similarity"):
    """
    Export a pairwise similarity or distance matrix to CSV.

    Parameters:
        matrix (list of lists): pairwise matrix
        sample_names (list): row/column labels
        filename (str): output filename
        label (str): type of matrix for logging
    """
    df = pd.DataFrame(matrix, index=sample_names, columns=sample_names)
    df.to_csv(filename, float_format='%.4f')
    print(f"{label} matrix saved: {filename}")


def run_analysis(raw_data, locus_names, dataset_label):
    """
    Run the full ISSR analysis pipeline for one dataset.

    Parameters:
        raw_data (dict): raw binary scoring data
        locus_names (list): names of loci
        dataset_label (str): label for output files (e.g., 'L1_4')
    """
    print(f"\n{'='*60}")
    print(f" ISSR ANALYSIS: {dataset_label.replace('_', ' ').upper()}")
    print(f"{'='*60}")

    # Step 1: Combine replicates
    combined = combine_replicates(raw_data)
    samples  = list(combined.keys())
    print(f"\nSamples: {', '.join(samples)}")
    print(f"Loci:    {', '.join(locus_names)}")
    print("\nCombined binary matrix (OR of replicate rows):")
    for s in samples:
        print(f"  {s}: {combined[s]}   total bands = {sum(combined[s])}")

    # Step 2: Similarity matrices
    jac_matrix, _ = build_similarity_matrix(combined, jaccard_similarity)
    dice_matrix, _ = build_similarity_matrix(combined, dice_similarity)

    # Step 3: Genetic distance (1 - Jaccard)
    dist_matrix = [[round(1 - jac_matrix[i][j], 4)
                    for j in range(len(samples))]
                   for i in range(len(samples))]

    print("\nJaccard Similarity Matrix:")
    header = f"{'':6}" + "".join(f"{s:8}" for s in samples)
    print(header)
    for i, si in enumerate(samples):
        row = f"{si:6}" + "".join(f"{jac_matrix[i][j]:8.4f}" for j in range(len(samples)))
        print(row)

    # Step 4: Band frequency and PIC
    freq_df = band_frequency_and_pic(combined, locus_names)
    print(f"\nBand Frequency and PIC:")
    print(freq_df.to_string(index=False))
    print(f"\nMean PIC: {freq_df['PIC'].mean():.4f}")
    print(f"Polymorphic loci: {(freq_df['PIC'] > 0).sum()} of {len(locus_names)}")

    # Step 5: UPGMA Newick tree
    newick = upgma_newick(samples, dist_matrix)
    print(f"\nNewick tree:\n{newick}")
    newick_file = f"ISSR_{dataset_label}_newick.nwk"
    with open(newick_file, 'w') as f:
        f.write(newick)
    print(f"Newick saved: {newick_file}")

    # Step 6: Dendrogram plot
    plot_dendrogram(
        sample_names = samples,
        dist_matrix  = dist_matrix,
        title        = f"UPGMA Dendrogram — ISSR {dataset_label.replace('_',' ')} (Jaccard Distance)",
        filename     = f"ISSR_{dataset_label}_dendrogram.png"
    )

    # Step 7: Export matrices and stats
    export_matrix_csv(jac_matrix,  samples, f"ISSR_{dataset_label}_jaccard_similarity.csv",  "Jaccard similarity")
    export_matrix_csv(dice_matrix, samples, f"ISSR_{dataset_label}_dice_similarity.csv",     "Dice similarity")
    export_matrix_csv(dist_matrix, samples, f"ISSR_{dataset_label}_genetic_distance.csv",    "Genetic distance")

    freq_df.to_csv(f"ISSR_{dataset_label}_band_frequency_PIC.csv", index=False)
    print(f"Band frequency/PIC saved: ISSR_{dataset_label}_band_frequency_PIC.csv")

    # Step 8: Total bands per sample summary
    summary = pd.DataFrame({
        'Sample'      : samples,
        'Total_Bands' : [sum(combined[s]) for s in samples],
        'Pct_Presence': [round(sum(combined[s]) / len(locus_names) * 100, 1)
                         for s in samples]
    })
    summary.to_csv(f"ISSR_{dataset_label}_sample_summary.csv", index=False)
    print(f"Sample summary saved: ISSR_{dataset_label}_sample_summary.csv")

    print(f"\nAnalysis for {dataset_label} complete.")
    return combined, jac_matrix, dist_matrix, freq_df


# =============================================================================
# 3. COMBINED ANALYSIS (Loci 1-8 together)
# =============================================================================

def combine_datasets(raw1, raw2):
    """
    Merge two datasets (e.g., Loci 1-4 and Loci 5-8) into one.
    Samples must be the same in both datasets.

    Parameters:
        raw1, raw2 (dict): raw data dictionaries

    Returns:
        dict: combined raw data with concatenated locus scores
    """
    combined_raw = {}
    for sample in raw1.keys():
        if sample in raw2:
            combined_raw[sample] = [
                raw1[sample][0] + raw2[sample][0],  # row 1 concatenated
                raw1[sample][1] + raw2[sample][1],  # row 2 concatenated
            ]
    return combined_raw


# =============================================================================
# 4. MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":

    print("ISSR Molecular Marker Analysis Pipeline")
    print("========================================")
    print("Author: gh | Postgraduate Researcher, Bioinformatics")
    print()

    # Run analysis for Loci 1-4
    combined_L1_4, jac_L1_4, dist_L1_4, freq_L1_4 = run_analysis(
        raw_data      = RAW_DATA_L1_4,
        locus_names   = LOCUS_NAMES_L1_4,
        dataset_label = "Loci1_4"
    )

    # Run analysis for Loci 5-8
    combined_L5_8, jac_L5_8, dist_L5_8, freq_L5_8 = run_analysis(
        raw_data      = RAW_DATA_L5_8,
        locus_names   = LOCUS_NAMES_L5_8,
        dataset_label = "Loci5_8"
    )

    # Run combined analysis (all 8 loci)
    combined_all_raw = combine_datasets(RAW_DATA_L1_4, RAW_DATA_L5_8)
    all_locus_names  = LOCUS_NAMES_L1_4 + LOCUS_NAMES_L5_8

    combined_ALL, jac_ALL, dist_ALL, freq_ALL = run_analysis(
        raw_data      = combined_all_raw,
        locus_names   = all_locus_names,
        dataset_label = "All_Loci"
    )

    print("\n" + "="*60)
    print(" ALL ANALYSES COMPLETE")
    print("="*60)
    print("\nOutput files generated:")
    print("  - *_jaccard_similarity.csv   (Jaccard pairwise matrix)")
    print("  - *_dice_similarity.csv      (Dice pairwise matrix)")
    print("  - *_genetic_distance.csv     (1 - Jaccard distance matrix)")
    print("  - *_band_frequency_PIC.csv   (Band freq + PIC per locus)")
    print("  - *_sample_summary.csv       (Total bands per sample)")
    print("  - *_newick.nwk               (Newick tree for MEGA/FigTree)")
    print("  - *_dendrogram.png           (UPGMA dendrogram figure)")
    print("\nTo visualise the Newick tree in MEGA 11:")
    print("  File → Open a File/Session → select .nwk file")
    print("  Tree Explorer opens automatically → View → Rectangular")
