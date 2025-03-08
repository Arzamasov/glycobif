import os
import pandas as pd

# Paths
annotation_folder = 'tmp/annotation'
statistics_path = 'database/fr_length_stats.tsv'
output_folder = 'tmp/annotation_filt'
non_passing_output_path = 'tmp/filtered_annotations.tsv'
default_coverage_threshold = 0.3

# Custom coverage thresholds for specific roles observed from manual curation of reference genomes
custom_coverage_thresholds = {
    "Putative sucrose MFS permease (TC 2.A.1.5)": 0.5,
    "Sucrose phosphorylase (EC 2.4.1.7)": 0.5,
    "Putative mannitol MFS permease (TC 2.A.1.18)": 0.5,
}

# Ensure the output directory exists
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Load functional role length data
stats_df = pd.read_csv(statistics_path, sep='\t', index_col='role')
print("Functional role length data loaded successfully.")

# DataFrame to hold all entries not passing the threshold
non_passing_entries = pd.DataFrame()

# Roles to always include
always_include_roles = [
    "Putative 2-dehydro-3-deoxygluconate kinase (EC 2.7.1.45)",
    "Putative N-acetylglucosamine-specific PTS system IIA component (TC 4.A.1)",
    "Putative N-acetylglucosamine-specific PTS system IIB component (TC 4.A.1)",
    "Putative N-acetylglucosamine-specific PTS system-2 IIA component (TC 4.A.1)",
    "Putative N-acetylglucosamine-specific PTS system-2 IIB component (TC 4.A.1)",
    "Putative transcriptional regulator of fucosylated HMO utilization (LacI family)",
    "Putative galactose MFS permease (TC 2.A.2)",
    "Lactose MFS permease (TC 2.A.2)",
    "Alpha-1,6-glucosidase (GH13_31)",
]

# Function to filter rows based on criteria
def filter_rows(df, stats_df):
    # Sum Region_len by Winner
    sum_len = df.groupby('Winner')['Region_len'].sum()
    
    # Calculate thresholds with custom values where applicable
    thresholds = {
        winner: custom_coverage_thresholds.get(winner, default_coverage_threshold) * stats_df.at[winner, 'median']
        for winner in sum_len.index if winner in stats_df.index
    }

    # Determine valid winners, including always include roles regardless of their sums
    valid_winners = [
        winner for winner in sum_len.index 
        if winner in stats_df.index and (sum_len[winner] >= thresholds.get(winner, 0) or winner in always_include_roles)
    ]

    passing_df = df[df['Winner'].isin(valid_winners)]
    non_passing_df = df[~df['Winner'].isin(valid_winners)]
    return passing_df, non_passing_df, sum_len, valid_winners

# Process each file in the annotation folder
for filename in os.listdir(annotation_folder):
    file_path = os.path.join(annotation_folder, filename)
    output_path = os.path.join(output_folder, filename)  # Save files in the output folder with the same name
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        passing_df, non_passing_df, sum_lengths, passed_winners = filter_rows(df, stats_df)
        
        # Append non-passing entries to the combined DataFrame
        non_passing_entries = pd.concat([non_passing_entries, non_passing_df], ignore_index=True)
        
        if not passing_df.empty:
            passing_df.to_csv(output_path, sep='\t', index=False)
            print(f"File {filename} processed: data passed the filtering criteria and saved.")
        else:
            print(f"File {filename} processed: no data passed the filtering criteria.")
    except UnicodeDecodeError:
        print(f"Skipping {filename} due to a Unicode decoding error.")
    except Exception as e:
        print(f"Failed to process {filename} due to an error: {e}")

# Save all non-passing entries to a separate combined file
if not non_passing_entries.empty:
    non_passing_entries.to_csv(non_passing_output_path, sep='\t', index=False)
    print(f"Combined file for non-passing entries saved at {non_passing_output_path}")

print("All files have been processed")
