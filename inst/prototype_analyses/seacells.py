#!/opt/homebrew/Caskroom/mambaforge/base/envs/seacells/bin/python


import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib as plt
import re
import os

ad = sc.read_h5ad('/Users/mcgaugheyd/git/plaeApp/inst/example_analyses/scEiaD_all_anndata.fix01.h5ad')

#cut down to one sample
pd.set_option('display.max_rows', None)
ad.obs[['Tissue','sample_accession','organism']].value_counts()
sample_counts = ad[ad.obs['Organ'].isin(['Eye'])].obs[['sample_accession']].value_counts(ascending = True)
existing = os.listdir('/Users/mcgaugheyd/data/eiad_seacells/')
existing_samples = [re.sub(r'.seacell_aggr.csv.gz|.obs.csv.gz', "", string = x) for x in existing]


for sample in sample_counts[sample_counts > 700].index.to_list():

    sample = str(sample)[2:len(sample)-4]
    
    if sample in existing_samples:
        print(sample + " skip!")
    else:
        print(sample)

        ad_sample = ad[ad.obs['sample_accession'].isin([sample])]

        # Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
        # This step should be performed after filtering 
        raw_ad = sc.AnnData(ad_sample.X)
        raw_ad.obs_names, raw_ad.var_names = ad_sample.obs_names, ad_sample.var_names
        ad_sample.raw = raw_ad

        # Normalize cells, log transform and compute highly variable genes
        sc.pp.normalize_per_cell(ad_sample)
        sc.pp.log1p(ad_sample)
        sc.pp.highly_variable_genes(ad_sample, n_top_genes=1500)

        # Compute principal components - 
        sc.tl.pca(ad_sample, n_comps=20, use_highly_variable=True)

        ## Core parameters 
        ## want ~ 75 cells per meta
        n_SEACells = int(ad_sample.shape[0] / 75)
        build_kernel_on = 'X_pca' 
        ## Additional parameters
        n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

        model = SEACells.core.SEACells(ad_sample, 
                    build_kernel_on=build_kernel_on, 
                    n_SEACells=n_SEACells, 
                    n_waypoint_eigs=n_waypoint_eigs,
                    convergence_epsilon = 1e-5)


        model.construct_kernel_matrix()
        M = model.kernel_matrix

        # Initialize archetypes
        model.initialize_archetypes()

        # fit model

        model.fit(min_iter=5, max_iter=200)
        #model.plot_convergence()

        labels,weights = model.get_soft_assignments()

        meta_aggr = SEACells.core.summarize_by_soft_SEACell(ad_sample, model.A_,summarize_layer='raw', minimum_weight=0.05)

        meta_aggr_df = pd.DataFrame(meta_aggr.X.toarray())
        meta_aggr_df.columns = meta_aggr.var_names
        meta_aggr_df.index = sample + '__' + meta_aggr_df.index.astype('str')

        pd.DataFrame(ad_sample.obs).to_csv('~/data/eiad_seacells/' + sample + '.obs.csv.gz')
        meta_aggr_df.to_csv('~/data/eiad_seacells/' + sample + '.seacell_aggr.csv.gz')

