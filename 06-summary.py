import pandas as pd
import sys
import glob

# ---------------------------
# Define method order based on context
# If '--pacbio' is passed, use a reduced set of phasing tools
# ---------------------------
if (len(sys.argv) == 2) and (sys.argv[1] == '--pacbio'):
    methods = ["TinkerHap", "WhatsHap", "HapCUT2"]
else:
    methods = ["TinkerHap+SI", "TinkerHap", "WhatsHap", "HapCUT2", "ShapeIT"]

# ---------------------------
# Load detailed evaluation data
# ---------------------------
li = []
df = pd.read_csv(glob.glob("./*-detailed.txt")[0], delimiter='\t')  # Read the first matching file
li.append(df)
df = pd.concat(li, axis=0, ignore_index=True)  # Combine all dataframes into one

# Ensure 'method' column is treated as categorical with a fixed order
df['method'] = pd.Categorical(df['method'], methods)

# ---------------------------
# Aggregate summary statistics per method
# ---------------------------
grouped = df.groupby('method', observed=False).agg({
    'method': 'count',                  # Count samples
    'sites': 'sum',                     # Total number of sites across all samples
    'errors': 'sum',                    # Total phasing errors
    'missing': 'sum',                   # Total missing sites
    'haplotypes_size_med': 'median'     # Median haplotype block size
})
grouped = grouped.rename(columns={'method': 'samples'})  # Rename 'method' to 'samples'

# ---------------------------
# Compute success metrics
# ---------------------------
grouped['Phased'] = 100 * ((grouped['sites'] - grouped['missing']) / grouped['sites'])
grouped['Phased-correctly'] = 100 * ((grouped['sites'] - grouped['missing'] - grouped['errors']) / grouped['sites'])

# Extract relevant columns for the output table
output = grouped[['samples', 'sites', 'Phased', 'Phased-correctly', 'errors', 'missing', 'haplotypes_size_med']]
data = output.transpose()  # Transpose so rows are metrics, columns are methods

# ---------------------------
# Add median common-site error rate per method
# ---------------------------
df = pd.read_csv(glob.glob("./*-common-sites.txt")[0], sep='\t')  # Load common-site comparison
common_sites = df.iloc[:, 2:].median()  # Median error rate per method
data.loc['common_errors'] = pd.Series(common_sites) * 100  # Convert to percent

# ---------------------------
# Add runtime and memory usage stats
# ---------------------------
df = pd.read_csv(glob.glob("./*-time-mem.txt")[0], sep='\t')  # Load timing/memory usage logs
time_mem_df_median = df.iloc[:, 1:].median()
time_mem_df_sum = df.iloc[:, 1:].sum()
time_mem_df_max = df.iloc[:, 1:].max()

runtime = {}
memory_usage = {}

# Parse and convert runtime/memory depending on context
for method in methods:
    if (len(sys.argv) == 2) and (sys.argv[1] == '--pacbio'):
        # Sum total runtime and peak memory (divided by 2 for paired runs)
        runtime[method] = (time_mem_df_sum[method + '-time'] / 2).round(1)
        memory_usage[method] = (time_mem_df_max[method + '-mem'] / 1024).round()  # Convert KB to MB
    else:
        # Use median runtime and memory
        runtime[method] = time_mem_df_median[method + '-time'].round(1)
        memory_usage[method] = (time_mem_df_median[method + '-mem'] / 1024).round()  # Convert KB to MB

data.loc['runtime'] = pd.Series(runtime, name='runtime')
data.loc['memory_usage'] = pd.Series(memory_usage)

# ---------------------------
# Clean up and format output
# ---------------------------
pd.options.mode.copy_on_write = True  # Ensures pandas doesn't warn on inplace operations
output.rename(columns={'method': 'samples'}, inplace=True)

# Set display options for cleaner console output
pd.set_option('display.max_colwidth', None)
pd.set_option('display.float_format', lambda x: '%.2f' % x)
pd.set_option('display.max_columns', None)

# ---------------------------
# Print final summary to stdout
# ---------------------------
print(data.to_string())
