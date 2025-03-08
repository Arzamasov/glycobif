import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from itertools import chain

class DiamondAnn:
    min_samples = 5             # Set minimum number of samples for DBSCAN clustering

    def __init__(self, hit_table):
        self.hit_table = hit_table.copy()       # Hits
        self.regions = np.empty((0, 2), int)     # Regions

    # Check if the potentially true hits exist and it is worth to continue
    def check(self):
        if self.hit_table.shape[0] < self.min_samples:         # will be no good support for DBSCAN
            return 0            # fail the check

        check_portion = self.hit_table.query('pident >= threshold')
        mcseed_hits_nullcheck = np.invert(check_portion['function'].isnull().values)
        # Mark all NaN values (no functional annotation) as False and annotated hits as true.
        # If at least one True - there are annotations for vote
        if np.any(mcseed_hits_nullcheck):
            return 1            # pass the check
        else:
            return 0            # fail the check

    # perform clustering of borders
    def DBSCAN(self, eps):
        # Make a list of boundaries
        aln_boundaries = np.concatenate([self.hit_table['qstart'].values, self.hit_table['qend'].values], axis=0)
        aln_boundaries = aln_boundaries.reshape(-1, 1)
        ## neighbors = NearestNeighbors(n_neighbors=5, algorithm='ball_tree').fit(aln_boundaries)
        ## max_dist = neighbors.kneighbors(aln_boundaries)[0].max()
        # eps could be selected based on kNN, but empirically it seems that 30 is good (approximately smallest domain)
        dbscan = DBSCAN(eps=eps, min_samples=self.min_samples).fit(aln_boundaries)
        # Reshape labels so labels for qstart will be in the first row and labels for qend in the second
        return dbscan.labels_.reshape(2, -1)

    # assign clusters to DBSCAN
    def AssignClusters(self, left_cluster, right_cluster):
        self.hit_table['clabel_left'] = left_cluster
        self.hit_table['clabel_right'] = right_cluster

        # Drop hits that DBSCAN predicted as outliers
        self.hit_table = self.hit_table.query('(clabel_left != -1) & (clabel_right != -1)')

        # Drop hits with ident below threshold
        self.hit_table = self.hit_table.query('pident >= threshold')

        # break if no hits left
        if self.hit_table.shape[0] == 0:
            return 1

        # Calculate mean boundary coordinate for each cluster
        cluster_list = np.concatenate((
            self.hit_table['clabel_left'].to_numpy(),
            self.hit_table['clabel_right'].to_numpy())
        )
        clusters = np.unique(cluster_list)
        coords = np.concatenate((
            self.hit_table['qstart'].to_numpy(),
            self.hit_table['qend'].to_numpy())
        )
        means = {cluster: np.mean(coords[cluster_list == cluster]).astype(int) for cluster in clusters}

        # Assign clusters to hits
        self.hit_table['left_mean'] = self.hit_table['clabel_left'].apply(lambda x: means[x])
        self.hit_table['right_mean'] = self.hit_table['clabel_right'].apply(lambda x: means[x])

        # DBSCAN was not able to separate borders
        if len(means) == 1:
            return 1

        # Assign clusters to regions
        boundaries = np.sort(np.array(list(means.values())))
        for i in range(len(boundaries) - 1):
            self.regions = np.append(self.regions, boundaries[i:i+2])

        # Remove regions less than 35aa
        self.regions = self.regions.reshape(-1, 2)
        save_regions = np.apply_along_axis(lambda x: (x[1] - x[0]) >= 35, 1, self.regions)
        self.regions = self.regions[save_regions, :]
        return 0

    # collect hits the overlap with the region
    def CollectHitsByRegion(self, left, right):
        return self.hit_table[['sseqid', 'pident', 'function']][
            np.logical_and(
                self.hit_table['left_mean'].values <= left,
                self.hit_table['right_mean'].values >= right
            )]

    # Voting for the most represented. Use Ident as weight
    def Voting(self, hit_ids):
        # Make annotations list
        hit_ids = hit_ids.explode('function')
        votes = hit_ids[['pident', 'function']]. \
            groupby('function', as_index=False). \
            agg({'pident': ['count', 'sum', 'max']})
        votes.columns = ['function', 'count', 'sum', 'max']

        # Select winners
        winner_votes = votes[votes['sum'] == votes['sum'].max()]

        # Output
        # 1. Winner functions
        winners = winner_votes['function'].values
        # 2. Number of functions in the voting
        n = votes['count'].sum()
        # 3. number of winners
        n_winners = winner_votes['count'].values
        # 4. Score
        score = votes['sum'].max() / votes['sum'].sum()
        # 5. Max ident of winner functions
        max_ident = votes[votes['function'].isin(winners)]['max'].values
        # 6. Idenetity threshold
        ident_thr = self.hit_table['threshold'].iloc[0]
        return winners, n, n_winners, score, max_ident, ident_thr
