'''
**Goal:** Compare error robustness across different genetic codes.
**Codes analyzed:** Standard Genetic Code (SGC), vertebrate mitochondrial code, CTG→Ser (*Candida*) code.
**Metric:** Error cost—average change in amino acid properties caused by single-nucleotide mutations.
**Approach:** For each codon, evaluate all one-base mutants, calculate error costs, and compare across codes using a distance matrix.

'''
import pandas as pd
from variant import sgc_code, mito_code, candida_code
import random
import copy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def compare_codes(sgc_code, mito_code, candida_code):
    codons = sorted(sgc_code.keys())
    df = pd.DataFrame({
        'SGC': [sgc_code[c] for c in codons],
        'Mito': [mito_code[c] for c in codons],
        'Candida': [candida_code[c] for c in codons]
    }, index=codons)

    #find diffrencec
    mask = ~((df['SGC'] == df['Mito']) & (df['SGC'] == df['Candida']))
    df_diff = df[mask]
    print("Codons with differences between codes:")
    print(df_diff)
    return df_diff

def grantham_distance(aa1, aa2):
    if aa1 == '*' or aa2 == '*':
        return 0  # או תן float('inf') אם תרצה להחמיר עבור stop codons
    return grantham_matrix.get((aa1, aa2), 0)

def error_cost(code_dict, grantham_matrix):
    bases = ['A', 'C', 'G', 'T']
    total = 0.0
    count = 0

    for codon, aa in code_dict.items():
        if aa == '*':
            continue  # דלג על קודוני עצירה

        codon = codon.upper()
        for i in range(3):  # עבור כל פוזיציה בקודון
            for base in bases:
                if base == codon[i]:
                    continue  # דלג אם אותו בסיס
                neighbor = codon[:i] + base + codon[i+1:]
                aa2 = code_dict.get(neighbor, '*')
                if aa2 == '*' or aa2 == aa:
                    continue  # דלג על neighbor שהוא stop או סינונימי
                dist = grantham_matrix.get((aa, aa2), grantham_matrix.get((aa2, aa), 0))
                total += dist
                count += 1
    if count == 0:
        return 0
    return total / count

def generate_one_swap_neighbors(code_dict, n_neighbors=1000, seed=None):
    if seed is not None:
        random.seed(seed)
    neighbors = []
    codons = [c for c in code_dict.keys() if code_dict[c] != '*']
    for _ in range(n_neighbors):
        code_copy = copy.deepcopy(code_dict)
        # בחר שני קודונים אקראיים שונים (לא עצירה)
        c1, c2 = random.sample(codons, 2)
        # החלף ביניהם
        code_copy[c1], code_copy[c2] = code_copy[c2], code_copy[c1]
        neighbors.append(code_copy)
    return neighbors

def error_costs_for_neighbors(neighbor_list, grantham_matrix):
    costs = []
    for code_variant in neighbor_list:
        cost = error_cost(code_variant, grantham_matrix)
        costs.append(cost)
    return costs

def plot_neighbor_costs(sgc_costs, mito_costs, candida_costs, orig_sgc, orig_mito, orig_candida):
    plt.figure(figsize=(8, 5))
    data = [sgc_costs, mito_costs, candida_costs]
    labels = ['SGC', 'Mitochondrial', 'Candida']
    box = plt.boxplot(data, labels=labels, patch_artist=True, showmeans=True)


    plt.scatter([1], [orig_sgc], color='red', marker='D', s=60, label='SGC original')
    plt.scatter([2], [orig_mito], color='red', marker='D', s=60, label='Mito original')
    plt.scatter([3], [orig_candida], color='red', marker='D', s=60, label='Candida original')

    plt.ylabel('Error Cost (Grantham)')
    plt.title('Distribution of Error Costs in One-Swap Neighbors')
    plt.legend()
    plt.tight_layout()
    plt.show()

def build_substitution_heatmap(code_dict, neighbors, grantham_matrix):
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    aa_idx = {aa: i for i, aa in enumerate(amino_acids)}
    matrix = np.zeros((20, 20))
    count_matrix = np.zeros((20, 20))

    for neighbor in neighbors:
        for codon in code_dict:
            orig_aa = code_dict[codon]
            new_aa = neighbor[codon]
            if orig_aa != '*' and new_aa != '*' and orig_aa != new_aa:
                i, j = aa_idx[orig_aa], aa_idx[new_aa]
                matrix[i, j] += grantham_matrix[(orig_aa, new_aa)]
                count_matrix[i, j] += 1

    # אפשר לחשב ממוצע או סכום; הנה ממוצע חילוף
    with np.errstate(divide='ignore', invalid='ignore'):
        avg_matrix = np.where(count_matrix > 0, matrix / count_matrix, 0)

    plt.figure(figsize=(10, 8))
    sns.heatmap(avg_matrix, xticklabels=amino_acids, yticklabels=amino_acids, cmap='YlOrRd', annot=False)
    plt.xlabel('To Amino Acid')
    plt.ylabel('From Amino Acid')
    plt.title('Average Grantham Distance of Amino Acid Substitutions\n(Across One-Swap Neighbors)')
    plt.tight_layout()
    plt.show()
    return avg_matrix, count_matrix


amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]

grantham_table = [
    #       A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V
    [   0, 112, 111, 126, 195,  91, 107,  60,  86,  94,  96, 106,  84, 113,  27,  99,  58, 148, 112,  64],  # A
    [ 112,   0,  86,  96, 180,  43,  54, 125,  29,  97, 102,  26,  91,  97, 103, 110,  71, 101,  77,  96],  # R
    [ 111,  86,   0,  23, 139,  46,  42,  80,  68, 149, 153,  94, 142, 158,  91,  46,  65, 174, 143, 133],  # N
    [ 126,  96,  23,   0, 154,  61,  45,  94,  81, 168, 172, 101, 160, 177, 108,  65,  85, 181, 160, 152],  # D
    [ 195, 180, 139, 154,   0, 154, 170, 159, 174, 198, 198, 202, 196, 205, 169, 112, 149, 215, 194, 192],  # C
    [  91,  43,  46,  61, 154,   0,  29,  87,  24, 109, 113,  53, 101, 116,  76,  68,  42, 130,  99,  96],  # Q
    [ 107,  54,  42,  45, 170,  29,   0,  98,  40, 134, 138,  56, 126, 140,  93,  80,  65, 152, 122, 121],  # E
    [  60, 125,  80,  94, 159,  87,  98,   0,  98, 135, 138, 127, 127, 153,  42,  56,  59, 184, 147, 109],  # G
    [  86,  29,  68,  81, 174,  24,  40,  98,   0,  94,  99,  32,  87, 100,  77,  89,  47, 115,  83,  84],  # H
    [  94,  97, 149, 168, 198, 109, 134, 135,  94,   0,   5, 102,  10,  21,  95, 142,  89,  61,  33,  29],  # I
    [  96, 102, 153, 172, 198, 113, 138, 138,  99,   5,   0, 107,  15,  22,  98, 145,  92,  61,  36,  32],  # L
    [ 106,  26,  94, 101, 202,  53,  56, 127,  32, 102, 107,   0,  95, 102, 103, 121,  78, 110,  85,  97],  # K
    [  84,  91, 142, 160, 196, 101, 126, 127,  87,  10,  15,  95,   0,  28,  87, 135,  81,  67,  36,  21],  # M
    [ 113,  97, 158, 177, 205, 116, 140, 153, 100,  21,  22, 102,  28,   0, 114, 155, 103,  40,  22,  50],  # F
    [  27, 103,  91, 108, 169,  76,  93,  42,  77,  95,  98, 103,  87, 114,   0,  74,  38, 147, 110,  68],  # P
    [  99, 110,  46,  65, 112,  68,  80,  56,  89, 142, 145, 121, 135, 155,  74,   0,  58, 177, 144, 124],  # S
    [  58,  71,  65,  85, 149,  42,  65,  59,  47,  89,  92,  78,  81, 103,  38,  58,   0, 128,  92,  69],  # T
    [ 148, 101, 174, 181, 215, 130, 152, 184, 115,  61,  61, 110,  67,  40, 147, 177, 128,   0,  37,  88],  # W
    [ 112,  77, 143, 160, 194,  99, 122, 147,  83,  33,  36,  85,  36,  22, 110, 144,  92,  37,   0,  55],  # Y
    [  64,  96, 133, 152, 192,  96, 121, 109,  84,  29,  32,  97,  21,  50,  68, 124,  69,  88,  55,   0]   # V
]

grantham_matrix = {}
for i, aa1 in enumerate(amino_acids):
    for j, aa2 in enumerate(amino_acids):
        grantham_matrix[(aa1, aa2)] = grantham_table[i][j]



# דוגמה לשימוש:

compare_codes(sgc_code, mito_code, candida_code)

print(grantham_distance('R','A'))

compare_codes(sgc_code, mito_code, candida_code)

sgc_ec = error_cost(sgc_code, grantham_matrix)
mito_ec = error_cost(mito_code, grantham_matrix)
candida_ec = error_cost(candida_code, grantham_matrix)

print("SGC error cost:", sgc_ec)
print("Mito error cost:", mito_ec)
print("Candida error cost:", candida_ec)

#make 1000 neighbors
sgc_neighbors = generate_one_swap_neighbors(sgc_code)
mito_neighbors = generate_one_swap_neighbors(mito_code)
candida_neighbors = generate_one_swap_neighbors(candida_code)

# score for neighbors
sgc_neighbor_costs = error_costs_for_neighbors(sgc_neighbors, grantham_matrix)
mito_neighbor_costs = error_costs_for_neighbors(mito_neighbors, grantham_matrix)
candida_neighbor_costs = error_costs_for_neighbors(candida_neighbors, grantham_matrix)

plot_neighbor_costs(sgc_neighbor_costs, mito_neighbor_costs, candida_neighbor_costs,
                   sgc_ec, mito_ec, candida_ec)


build_substitution_heatmap(sgc_code, sgc_neighbors, grantham_matrix)
build_substitution_heatmap(mito_code, mito_neighbors, grantham_matrix)
build_substitution_heatmap(candida_code, candida_neighbors, grantham_matrix)

