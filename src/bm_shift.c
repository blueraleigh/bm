#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <phy.h>


/* struct to track sufficient statistics for computing the
** likelihood under Brownian motion of a set of phylogenetic
** independent contrasts */
struct bm {
    // number of contrasts
    int n;
    // sum of squared contrasts scaled by their variances (or, rather, by a
    // scalar that is proportional to their variance)
    double su;
    // sum of the logarithms of contrast variances (or, rather, a
    // scalar that is proportional to their variance)
    double sv;
};


static void bm_push(
    double u,       /* a contrast */
    double vu,      /* sum of adjusted branch lengths for contrast */
    struct bm *bm)
{
    bm->n += 1;
    if (bm->n == 1)
    {
        bm->su = (u*u)/vu;
        bm->sv = log(vu);
    }
    else
    {
        bm->su += (u*u)/vu;
        bm->sv += log(vu);
    }
}


static void bm_merge(struct bm *a, struct bm *b, struct bm *c)
{
    c->n = a->n + b->n;
    c->su = a->su + b->su;
    c->sv = a->sv + b->sv;
}


static double bm_loglikelihood(struct bm *bm)
{
    // MLE of the Brownian motion variance parameter
    double var = bm->su / bm->n;
    if (var > 0)
    {
        double lnL = bm->n * M_LN_2PI + bm->n * log(var) + bm->sv + bm->su / var;
        return -0.5 * lnL;
    }
    return 0;
}


struct dp {
    /* number of terminals */
    int ntl;

    /* only filled entries are valid */
    int filled;

    /* interpolated node state under Brownian motion */
    double x;

    /* variance of x under Brownian motion */
    double vx;

    /* log likelihood of shifted processes */
    double logp;

    /* log likelihood root process and shifted processes */
    double score;

    /* sufficient statistics log for the root process */
    struct bm bm;

    /* index of item in the dp array of left descendant that was
    ** used to form this entry */
    int j;

    /* index of item in the dp array of right descendant that was
    ** used to form this entry */
    int k;
};


static void downpass(struct phy *phy, int *n_edge)
{
    int i;
    int j;
    int k;

    double u;
    double vu;
    double x;
    double vx;
    double lvx;
    double rvx;
    double score;

    struct phy_node *lf;
    struct phy_node *rt;
    struct phy_node *node;
    struct phy_cursor *cursor;

    struct dp *vec;
    struct dp *l_vec;
    struct dp *r_vec;

    struct bm bm;
    struct bm empty = {0, 0, 0};

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        lf = phy_node_lfdesc(node);
        rt = phy_node_rtdesc(node);

        vec = (struct dp *)phy_node_data(node);
        l_vec = (struct dp *)phy_node_data(lf);
        r_vec = (struct dp *)phy_node_data(rt);

        u = l_vec[0].x - r_vec[0].x;
        vu = l_vec[0].vx + r_vec[0].vx;
        x = (r_vec[0].vx * l_vec[0].x + l_vec[0].vx * r_vec[0].x) / vu;
        vx = phy_node_brlen(node) + (l_vec[0].vx * r_vec[0].vx) / vu;
        bm_merge(&(l_vec[0].bm), &(r_vec[0].bm), &(vec[0].bm));
        bm_push(u, vu, &(vec[0].bm));
        vec[0].ntl = 1;
        vec[0].x = x;
        vec[0].vx = vx;
        vec[0].logp = bm_loglikelihood(&(vec[0].bm));
        vec[0].score = vec[0].logp;
        vec[0].filled = 1;
        vec[0].j = -1;
        vec[0].k = -1;

        for (j = 0; j <= n_edge[phy_node_index(lf)]; ++j)
        {
            for (k = 0; k <= n_edge[phy_node_index(rt)]; ++k)
            {
                if (!l_vec[j].filled || !r_vec[k].filled)
                    continue;

                i = j + k + 2;

                // if j or k is 0 we want to use the original
                // rather than adjusted branch length. that's because
                // a value of 0 means the descendant clade does not
                // belong to the root process and we are treating the
                // the MLE of the root state of the descendant clade
                // as fixed.
                lvx = j == 0 ? phy_node_brlen(lf) : l_vec[j].vx;
                rvx = k == 0 ? phy_node_brlen(rt) : r_vec[k].vx;

                u = l_vec[j].x - r_vec[k].x;
                vu = lvx + rvx;

                x = (rvx * l_vec[j].x + lvx * r_vec[k].x) / vu;
                vx = phy_node_brlen(node) + (lvx * rvx) / vu;

                if (j == 0 && k == 0)
                    bm_merge(&empty, &empty, &bm);
                else if (j == 0 && k > 0)
                    bm_merge(&empty, &(r_vec[k].bm), &bm);
                else if (j > 0 && k == 0)
                    bm_merge(&(l_vec[j].bm), &empty, &bm);
                else
                    bm_merge(&(l_vec[j].bm), &(r_vec[k].bm), &bm);

                bm_push(u, vu, &bm);
                score = l_vec[j].logp + r_vec[k].logp + bm_loglikelihood(&bm);

                if (!vec[i].filled || score > vec[i].score)
                {
                    vec[i].ntl = l_vec[j].ntl + r_vec[k].ntl;
                    vec[i].filled = 1;
                    vec[i].score = score;
                    vec[i].j = j;
                    vec[i].k = k;
                    vec[i].logp = l_vec[j].logp + r_vec[k].logp;
                    vec[i].x = x;
                    vec[i].vx = vx;
                    vec[i].bm = bm;
                }
            }
        }
    }
}


static void backtrack(int index, struct phy_node *node, struct phy *phy,
    double aic_w, double bg_rate, double *rate, int *i, int *shift)
{
    struct dp *vec = (struct dp *)phy_node_data(node);

    if (vec[index].j != -1 && vec[index].k != -1)
    {
        rate[phy_node_index(node)] += aic_w * bg_rate;
        backtrack(vec[index].j, phy_node_lfdesc(node), phy, aic_w, bg_rate, 
            rate, i, shift);
        backtrack(vec[index].k, phy_node_rtdesc(node), phy, aic_w, bg_rate, 
            rate, i, shift);
    }
    else
    {
        if (phy_node_istip(node))
        {
            rate[phy_node_index(node)] += aic_w * bg_rate;
        }
        else
        {
            struct phy_cursor *cursor;

            if (i && shift && node != phy_root(phy))
                shift[(*i)++] = phy_node_index(node)+1;

            rate[phy_node_index(node)] += aic_w * bg_rate;
            cursor = phy_cursor_prepare(phy, node, ALL_NODES, PREORDER);
            phy_cursor_step(cursor);
            while ((node = phy_cursor_step(cursor)) != 0)
                rate[phy_node_index(node)] += aic_w * (vec[index].bm.su / vec[index].bm.n);
        }
    }
}


static void dp_init(double *x, int *n_edge, struct phy *phy)
{
    int i;
    struct phy_node *node;
    struct phy_cursor *cursor;
    struct dp *dp;

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        if (phy_node_istip(node))
        {
            n_edge[phy_node_index(node)] = 0;
        }
        else
        {
            n_edge[phy_node_index(node)] =
                n_edge[phy_node_index(phy_node_lfdesc(node))] +
                n_edge[phy_node_index(phy_node_rtdesc(node))] + 2;
        }
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        dp = calloc(n_edge[phy_node_index(node)] + 1, sizeof(*dp));
        phy_node_set_data(node, (void *)dp, &free);
    }

    for (i = 0; i < phy_ntip(phy); ++i)
    {
        node = phy_node_get(phy, i);
        dp = (struct dp *)phy_node_data(node);
        dp[0].x = x[i];
        dp[0].vx = phy_node_brlen(node);
        dp[0].ntl = 1;
        dp[0].filled = 1;
        dp[0].j = -1;
        dp[0].k = -1;
    }
}


static double aic(double lnL, int k)
{
    return -2*lnL + 2*k;
}


static double bic(double lnL, int k, int n)
{
    return -2*lnL + log(n)*k;
}


SEXP C_bm_shift(SEXP x, SEXP rtree)
{
    int i;
    int n;
    int best = 0;
    double w = 0;
    double score = R_NegInf;
    struct dp *dp;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int n_edge[phy_nnode(phy)];
    double rate[phy_nnode(phy)];

    dp_init(REAL(x), n_edge, phy);
    downpass(phy, n_edge);
    dp = (struct dp *)phy_node_data(phy_root(phy));

    double aic_score[n_edge[phy_node_index(phy_root(phy))]+1];
    double aic_w[n_edge[phy_node_index(phy_root(phy))]+1];
    double min_aic_score = aic_score[0] = aic(dp[0].score, 2*dp[0].ntl);

    for (i = 1; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        if (dp[i].filled)
        {
            aic_score[i] = aic(dp[i].score, 2*dp[i].ntl);
            if (aic_score[i] < min_aic_score)
            {
                min_aic_score = aic_score[i];
                best = i;
            }
        }
    }

    for (i = 0; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        aic_w[i] = 0;
        if (dp[i].filled)
            w += aic_w[i] = exp(-0.5 * (aic_score[i] - min_aic_score));
    }

    memset(rate, 0, phy_nnode(phy) * sizeof(double));
    for (i = 0; i <= n_edge[phy_node_index(phy_root(phy))]; ++i)
    {
        if (dp[i].filled)
        {
            aic_w[i] /= w;
            backtrack(i, phy_root(phy), phy, aic_w[i], dp[i].bm.su / dp[i].bm.n, 
                rate, 0, 0);
        }
    }

    SEXP ans = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, phy_nnode(phy)));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, phy_nnode(phy)));

    memcpy(REAL(VECTOR_ELT(ans, 0)), rate, phy_nnode(phy) * sizeof(double));
    memcpy(REAL(VECTOR_ELT(ans, 1)), aic_w, phy_nnode(phy) * sizeof(double));

    UNPROTECT(1);
    return ans;
}


SEXP C_bm_shift_backtrack(SEXP index, SEXP rtree)
{
    int i = *INTEGER(index);
    int j = 0;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct dp *dp = (struct dp *)phy_node_data(phy_root(phy));
    int shift[phy_nnode(phy)];
    SEXP rate = PROTECT(allocVector(REALSXP, phy_nnode(phy)));
    memset(REAL(rate), 0, phy_nnode(phy) * sizeof(double));
    if (dp && dp[i].filled) {
        backtrack(i, phy_root(phy), phy, 1, dp[i].bm.su / dp[i].bm.n, 
            REAL(rate), &j, shift);
    }

    SEXP shifts = PROTECT(allocVector(INTSXP, j));

    memcpy(INTEGER(shifts), shift, j*sizeof(int));

    setAttrib(rate, install("shifts"), shifts);

    UNPROTECT(2);
    return rate;
}
