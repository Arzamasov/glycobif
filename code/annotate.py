import diamondann
from collections import defaultdict
from glob import glob
from os import path
import argparse
import csv
import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from multiprocessing import Pool
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import multiprocessing

multiprocessing.set_start_method('fork')

parser = argparse.ArgumentParser(description='Annotate DIAMOND outfmt 6 output with domain fusion handling')
parser.add_argument('-d', '--dmnd', help='Pathway to the DIAMOND outfmt 6 result')
parser.add_argument('-a', '--ann', help='file with annotations')
parser.add_argument('-f', '--fmt', nargs='*',
                    default=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'],
                    type=str,
                    help='Specify outfmt 6 columns for DIAMOND output. Default is the DIAMOND outfmt 6 default'
                    )
parser.add_argument('-p', '--processes', default=1, help='Number of processes. Default: 1.')
parser.add_argument('-i', '--ident', default=0.8, nargs='?', help='Identity precentile that should be set as a threshold. Defauls: 0.8 (80%%).')
parser.add_argument('-m', '--min_samples', default=5, help='Minimum number of samples for DBSCAN clustering')
parser.add_argument('-o', '--out', help='output file')
parser.add_argument('--dir', action='store_true', help='If flag is set then inpit pathway from --dmnd parameter is a'
                                                       ' directory with DIAMOND search results. Output files would be '
                                                       'named as input files in the folder specified in --out '
                                                       'parameter.')
parser.add_argument('--skip_present', action='store_false', help='Skip file if the corresponding output exists. Default: False.')
parser.add_argument('--kde', action='store_true', help='Do estimation of the threshold based on a highest value of minimum of gaussian Kernel Density Estimator of identities.')
args = parser.parse_args()

def annotate(query_hit_table):
    # new instance of the object for a hit
    dmnd_ann = diamondann.DiamondAnn(query_hit_table)

    # check if we should continue (if any hit above identity threshold exists)
    if dmnd_ann.check():
        eps = 30
        labels = dmnd_ann.DBSCAN(eps=eps)

        # check if DBSCAN clustered borders in at least 2 groups
        while len(np.unique(labels)) == 1 and eps - 5 > 0:
            eps = eps - 5
            labels = dmnd_ann.DBSCAN(eps=eps)

        # check that DBSCAN was able to separate borders, excluding -1 as outliers
        if len(np.unique(labels[labels > -1])) == 1:
            return None

        # test if any good hits are kept
        if dmnd_ann.hit_table.shape[0] == 0:
            return None

        # AssignClusters return 1 if No hits are above identity threshold
        if dmnd_ann.AssignClusters(left_cluster=labels[0, :], right_cluster=labels[1, :]):
            return None

        annotations = pd.DataFrame()
        id_a = np.array([])
        winners_a = np.array([])
        score_a = np.array([])
        region_len_a = np.array([])
        win_votes_a = np.array([])
        all_votes_a = np.array([])
        region_max_ident_a = np.array([])
        rleft_a = np.array([])
        rright_a = np.array([])
        ident_thr_a = np.array([])

        for i in range(dmnd_ann.regions.shape[0]):
            hits = dmnd_ann.CollectHitsByRegion(left=dmnd_ann.regions[i][0], right=dmnd_ann.regions[i][1])
            # if region has no hits above identity threshold we should skip it
            if hits.shape[0] == 0:
                break
            winners, n, n_winners, score, max_ident, ident_thr = dmnd_ann.Voting(hit_ids=hits)

            if winners.size:
                # Query ID
                id_a = np.append(id_a, np.full(winners.size, query_hit_table['qseqid'].iloc[0]))
                # Winner functions
                winners_a = np.append(winners_a, winners)
                # Score
                score_a = np.append(score_a, np.full(winners.size, score))
                # Length of the region
                region_len_a = np.append(region_len_a, np.full(winners.size, dmnd_ann.regions[i][1] - dmnd_ann.regions[i][0]))
                # Number of win functions
                win_votes_a = np.append(win_votes_a, n_winners)
                # Number of all votes
                all_votes_a = np.append(all_votes_a, np.full(winners.size, n))
                # Max ident of winners
                region_max_ident_a = np.append(region_max_ident_a, max_ident)
                # Left boundary of region
                rleft_a = np.append(rleft_a, np.full(winners.size, dmnd_ann.regions[i, 0]))
                # Right boundary of region
                rright_a = np.append(rright_a, np.full(winners.size, dmnd_ann.regions[i, 1]))
                # Identity threshold
                ident_thr_a = np.append(ident_thr_a, np.full(winners.size, ident_thr))

        # Remove unassigned winners
        if winners_a.size:
            drop_empty = np.where(winners_a == '')
            id_a = np.delete(id_a, drop_empty)
            winners_a = np.delete(winners_a, drop_empty)
            score_a = np.delete(score_a, drop_empty)
            region_len_a = np.delete(region_len_a, drop_empty)
            win_votes_a = np.delete(win_votes_a, drop_empty)
            all_votes_a = np.delete(all_votes_a, drop_empty)
            region_max_ident_a = np.delete(region_max_ident_a, drop_empty)
            rleft_a = np.delete(rleft_a, drop_empty)
            rright_a = np.delete(rright_a, drop_empty)
            ident_thr_a = np.delete(ident_thr_a, drop_empty)

            if winners_a.size:
                annotations['ID'] = id_a
                annotations['Winner'] = winners_a
                annotations['Score'] = score_a
                annotations['Region_len'] = region_len_a
                annotations['Win_votes'] = win_votes_a
                annotations['All_votes'] = all_votes_a
                annotations['Region_max_ident'] = region_max_ident_a
                annotations['R_left'] = rleft_a
                annotations['R_right'] = rright_a
                annotations['Identity_thr'] = ident_thr_a

                return annotations


# Estimate the highest minimum of a KDE function for identities
def KDE_h_minimum(query_hit_table):
    ident = query_hit_table['pident'].to_numpy()
    kde = KernelDensity(kernel='gaussian', bandwidth=3).fit(ident.reshape(-1, 1))
    s = np.linspace(0, 100, 100)
    e = np.exp(kde.score_samples(s.reshape(-1, 1)))
    mi = argrelextrema(e, np.less)[0]
    mi = np.append(np.min(ident), mi)
    return pd.DataFrame({'qseqid': [query_hit_table['qseqid'].iloc[0]], 'threshold': [mi[-1]]})


# load mcSEED functions table as default dictionary
# File required to have structure:
# - No header
# Columns: Organism, PegID, Functional role, Name, Classificator1, C2, C3, C4
# Default value of dictionary is the empty string
# Dictionary [PegID] => [Functional role]
print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Reading reference functonal roles")

mcseed = pd.read_csv(args.ann,
                     header=None,
                     names=['sseqid', 'function', 'name', 'subsystem', 'm1', 'm2', 'm3', 'organism'],
                     usecols=['sseqid', 'function'],
                     sep='\t')

mcseed['function'] = mcseed['function'].apply(lambda x: x.split(' / '))

if mcseed.shape[0] == 0:
    print("WARNING: Dictionary of functional roles was not provided")

print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Initializing search")
diamondann.DiamondAnn.min_samples = int(args.min_samples)

# Collect files with DIAMOND output
if args.dir:
    dmnd_files = glob(path.join(args.dmnd, '*'))
else:
    dmnd_files = [args.dmnd]

for dmnd_file in dmnd_files:
    # set output file
    if args.dir:
        out_file = path.join(args.out, path.basename(dmnd_file))
    else:
        out_file = args.out

    # Skip file
    if args.skip_present:
        if path.exists(out_file):
            print(f'{out_file} exists. Skip annotation.')
            continue

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Reading DIAMOND output")
    # Load DIAMOND table with hits
    dmnd_table = pd.read_csv(dmnd_file,
                             header=None,
                             names=args.fmt,
                             sep='\t',
                             usecols=['qseqid','sseqid','pident','qstart','qend'],
                             dtype={'qseqid': str, 'sseqid': str, 'pident': np.float32, 'qstart': np.uint16, 'qend': np.uint16}
                             )
	
	
    if len(dmnd_table) == 0:
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": No hits found")
        sys.exit()

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Assign function to DIAMOND hits")

    # Assign functions to hits
    dmnd_table = dmnd_table.merge(mcseed, how='left', on='sseqid')
    dmnd_table = dmnd_table.fillna('')

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Calculate thresholds for each protein")

    # Split dataframe by query names
    query_split = dmnd_table.groupby('qseqid')

    # Assign thresholds for ident
    if args.kde:
        print("Use KDE algorithm to find threshold")
        with Pool(processes=int(args.processes)) as p:
            max_ = query_split.ngroups
            thr_table = list(tqdm(p.imap(KDE_h_minimum, [group for name, group in query_split]), total=max_))
            p.close()
            p.join()
        thr_table = pd.concat(thr_table, ignore_index=True)
    else:
        print(f"Use {args.ident} precentile for threshold")
        thr_table = dmnd_table[['qseqid', 'pident']].groupby('qseqid', as_index=False). \
            quantile(float(args.ident)). \
            rename(columns={'pident': 'threshold'})
    dmnd_table = dmnd_table.merge(thr_table, how='left', on='qseqid')

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Start annotation")
    annotations_out = list()

    # Repeat grouping to include thresholds
    query_split = dmnd_table.groupby('qseqid')

    with Pool(processes=int(args.processes)) as p:
        max_ = query_split.ngroups
        annotations_out = list(tqdm(p.imap(annotate, [group for name, group in query_split]), total=max_))
        p.close()
        p.join()

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Filter empty annotations")
    annotations_out = [ann for ann in annotations_out if ann is not None]
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Write output to file")
    if len(annotations_out):
        annotations_out = pd.concat(annotations_out, ignore_index=True)
        annotations_out = annotations_out.astype({'Region_len': np.uint16, 'Win_votes': np.uint16,
                                                  'All_votes': np.uint16, 'R_left': np.uint16, 'R_right': np.uint16})
        annotations_out.to_csv(out_file, sep='\t', index=False, quoting=csv.QUOTE_NONE, float_format='%.2f')
    else:
        pd.DataFrame(
            {
                'ID': [],
                'Winner': [],
                'Score': [],
                'Region_len': [],
                'Win_votes': [],
                'All_votes': [],
                'Region_max_ident': [],
                'R_left': [],
                'R_right': [],
                'Identity_thr': []
            }
        ).to_csv(out_file + '.ann', sep='\t', index=False, quoting=csv.QUOTE_NONE, float_format='%.2f')

    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ": Finished")