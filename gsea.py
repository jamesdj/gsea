import itertools

import numpy as np
import pandas as pd


def consecutive_pairs(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def max_abs(x):
    return x[np.argmax(np.abs(x))]


def _base_gsea(ranked_genes, sets_to_genes, collect_func, alpha=0.75):
    """
    Basic idea:
    -make weighted ecdf for hits
    -make ecdf for misses
    -take elementwise difference
    -collect_func() to get result (max_abs() for gsea, np.sum() for ssgsea)

    Might still be able to speed up using array funcs rather than iterating?
    """
    n_genes = len(ranked_genes)
    ranks = list(range(n_genes))
    gene_to_rank = dict(zip(ranked_genes, ranks))
    enrichment_scores = {}
    for set_name, set_genes in sets_to_genes.items():
        ranked_set_genes = [gene for gene in set_genes if gene in gene_to_rank]
        n_non_set_genes = float(n_genes - len(ranked_set_genes))
        hit_ranks = [gene_to_rank[gene] for gene in ranked_set_genes]
        misses = np.ones_like(ranks)
        misses[hit_ranks] = 0
        cum_misses = np.cumsum(misses)
        miss_ecdf = cum_misses / n_non_set_genes
        cum_hits = np.zeros_like(ranks)
        if len(hit_ranks) > 0:
            cum_hit_sum = 0
            sorted_hit_ranks = sorted(hit_ranks)
            # add one so ranks to weight start from 1, not zero
            # however, it's convenient to start at zero otherwise so I can index using the ranks
            weighted_ranks = (np.array(sorted_hit_ranks) + 1) ** alpha
            hit_rank_pairs = list(consecutive_pairs(sorted_hit_ranks))  # given [a, b, c, d] yields (a, b), (b, c), (c, d)
            for i, (idx1, idx2) in enumerate(hit_rank_pairs):
                cum_hit_sum += weighted_ranks[i]
                cum_hits[idx1:idx2] = cum_hit_sum
            cum_hit_sum += weighted_ranks[-1]
            cum_hits[sorted_hit_ranks[-1]:] = cum_hit_sum
            weighted_hit_ecdf = cum_hits / cum_hit_sum
        else:
            weighted_hit_ecdf = cum_hits  # still np.zeros_like(ranks)
        ecdf_dif = np.subtract(weighted_hit_ecdf, miss_ecdf)
        enrichment_score = collect_func(ecdf_dif)
        enrichment_scores[set_name] = enrichment_score
    return pd.Series(enrichment_scores)


def ssgsea(exp_data, sets_to_genes, alpha=0.75):
    """
    Single-sample GSEA as described in Barbie et al. (2009)
    :param exp_df: Pandas DataFrame or Series of expression values, (n_samples, n_genes) or (n_genes,)
    :param sets_to_genes: dictionary with set names as keys and sets of genes as values, e.g. {'set1': {'g1', 'g2'}}
    :param alpha: float, weighting factor between zero and one.
        Smaller values give more weight to top/bottom of list.
    :return: Pandas DataFrame or Series of expression projected onto gene sets
    """
    if isinstance(exp_data, pd.Series):
        return ssgsea_per_sample(exp_data, sets_to_genes, alpha=alpha)
    elif isinstance(exp_data, pd.DataFrame):
        return exp_data.apply(ssgsea_per_sample, axis=1, args=(sets_to_genes, alpha))
    else:
        raise ValueError("exp_data must be Pandas DataFrame or Series")


def ssgsea_per_sample(exp_series, sets_to_genes, alpha=0.75):
    sorted_exp_series = exp_series.sort_values(ascending=False)
    enrichment_scores = _base_gsea(sorted_exp_series.index, sets_to_genes, collect_func=np.sum, alpha=alpha)
    enrichment_scores /= 0.5 * len(exp_series)  # maximum possible score
    return enrichment_scores


def gsea(ranked_genes, sets_to_genes, alpha=0.75):
    """
    GSEA as described in Subramanian et al. (2005)
    :param ranked_genes: list of gene names
    :param sets_to_genes: dictionary with set names as keys and sets of genes as values, e.g. {'set1': {'g1', 'g2'}}
    :param alpha: float, weighting factor between zero and one.
        Smaller values give more weight to top/bottom of list.
    :return: Pandas Series of enrichment scores for each gene set.
    """
    enrichment_scores = _base_gsea(ranked_genes, sets_to_genes, collect_func=max_abs, alpha=alpha)
    return enrichment_scores


def gsea_nes_fdr(es, shuffled_phenotype_gseas):
    # according to Pablo, the most correct permutation is of the phenotype upstream of producing the ranked list
    res = pd.concat(shuffled_phenotype_gseas, axis=1).astype(float)
    nneg_res_means = res[res >= 0].mean(axis=1)
    neg_res_means = res[res <= 0].mean(axis=1)
    nres = res.copy().T
    nres[res.T >= 0] /= nneg_res_means
    nres[res.T <= 0] /= neg_res_means.abs()
    nres = nres.T
    nes = es.loc[nres.index]
    nes[nes >= 0] /= nneg_res_means
    nes[nes < 0] /= neg_res_means.abs()
    nes.name = 'NES'
    ge = nres.ge(es, axis=0).loc[es.index]
    more_extreme = ge.copy()
    le = nres.le(es, axis=0).loc[es.index]
    more_extreme[es < 0] = le[es < 0]
    pvals = more_extreme.sum(axis=1) / more_extreme.shape[1]
    pvals.name = 'p'
    wpval = pd.concat([nes, pvals], axis=1)
    wpval_pos = wpval[wpval.NES >= 0.0].sort_values(by='p')
    wpval_neg = wpval[wpval.NES < 0.0].sort_values(by='p', ascending=False)
    return pd.concat([wpval_pos, wpval_neg], axis=0)


