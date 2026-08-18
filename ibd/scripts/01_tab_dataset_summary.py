import os
import csv
import sys

# Results directory may be given on the command line; the revision run is the default.
root = sys.argv[1] if len(sys.argv) > 1 else 'results/revision_v2'
directories = [f'{root}/{name}' for name in ['uc_vs_hc', 'uninf_vs_hc', 'cd_vs_hc', 'uc_vs_cd']]

# Dictionary to store the aggregated counts for each dataset
dataset_counts = {}

for directory in directories:
    file_path = os.path.join(directory, "dataset_summary.csv")
    
    # Infer conditions from the directory name
    # e.g., "uc_vs_hc" -> condition1 = "uc", condition2 = "hc"
    dir_name = os.path.basename(directory)
    if "_vs_" in dir_name:
        condition1, condition2 = dir_name.split("_vs_")
    else:
        continue  # Skip if directory name doesn't match expected format
    
    if os.path.exists(file_path):
        with open(file_path, mode='r', encoding='utf-8') as f:
            reader = csv.reader(f)
            
            # Skip the header row
            next(reader, None)
            
            for row in reader:
                # Ensure the row has at least 4 columns (dataset, n_total, n_cond1, n_cond2)
                if len(row) < 4:
                    continue
                
                dataset = row[0]
                count_cond1 = row[2] # 3rd column
                count_cond2 = row[3] # 4th column
                
                if dataset not in dataset_counts:
                    dataset_counts[dataset] = {}
                
                # Assign counts based on the inferred conditions
                dataset_counts[dataset][condition1] = count_cond1
                dataset_counts[dataset][condition2] = count_cond2
    else:
        print(f"% Warning: File not found: {file_path}")

# Platform per dataset, from the acquisition tables
platforms = {}
for table in ['data/DAQ.csv', 'data/DAQ_rnaseq.csv']:
    if os.path.exists(table):
        with open(table, encoding='utf-8') as f:
            for row in csv.DictReader(f, delimiter=';'):
                platforms[row['Accession number']] = row['Platform']

# Generate the LaTeX table
print("% Remember to include \\usepackage{booktabs} in your LaTeX preamble")
print("\\begin{table}[htbp]")
print("\\centering")
print("\\begin{tabular}{llccccc}") 
print("\\toprule")
print("Dataset & Platform & $n_{UC}$ & $n_{uninfl}$ & $n_{CD}$ & $n_{HC}$ \\\\")
print("\\midrule")

# Variables to keep track of totals
total_uc = 0
total_uninfl = 0
total_cd = 0
total_hc = 0

for dataset, counts in dataset_counts.items():
    # Fetch counts, defaulting to '-' if the dataset doesn't have that specific group
    n_uc = counts.get('uc', '-')
    n_uninfl = counts.get('uninf', '-')
    n_cd = counts.get('cd', '-')
    n_hc = counts.get('hc', '-')
    
    # Add to totals (converting string to int, ignoring '-')
    if n_uc != '-': total_uc += int(n_uc)
    if n_uninfl != '-': total_uninfl += int(n_uninfl)
    if n_cd != '-': total_cd += int(n_cd)
    if n_hc != '-': total_hc += int(n_hc)
    
    print(f"{dataset} & {platforms.get(dataset, '')} & {n_uc} & {n_uninfl} & {n_cd} & {n_hc} \\\\")

# Print the Total row
print("\\midrule")
print(f"\\textbf{{Total}} & & \\textbf{{{total_uc}}} & \\textbf{{{total_uninfl}}} & \\textbf{{{total_cd}}} & \\textbf{{{total_hc}}} \\\\")

print("\\bottomrule")
print("\\end{tabular}")
print("\\caption{Summary of datasets and sample counts per group.}")
print("\\label{tab:dataset_summary}")
print("\\end{table}")