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


static void downpass(double *x, double *vx, double *u, struct phy *phy, 
    struct bm *bm)
{
    int i;
    int j;
    int k;
    struct phy_node *node;
    struct phy_cursor *cursor;

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        i = phy_node_index(node);
        j = phy_node_index(phy_node_lfdesc(node));
        k = phy_node_index(phy_node_rtdesc(node));

        u[i] = x[j] - x[k];
        x[i] = (vx[k]*x[j] + vx[j]*x[k]) / (vx[j] + vx[k]);
        vx[i] += (vx[j] * vx[k]) / (vx[j] + vx[k]);
    
        bm_push(u[i], vx[j] + vx[k], bm);
    }
}
    

SEXP C_bm_pic(SEXP y, SEXP rtree)
{
    struct bm bm = {0, 0, 0};
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    
    int i;
    int j;
    int k;
    int ntip = phy_ntip(phy);
    int nnode = phy_nnode(phy);

    double x[nnode];
    double u[nnode];
    double vx[nnode];

    struct phy_node *node;
    struct phy_cursor *cursor;

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
        vx[phy_node_index(node)] = phy_node_brlen(node);

    memcpy(x, REAL(y), ntip * sizeof(double));
    
    downpass(x, vx, u, phy, &bm);

    SEXP ret = PROTECT(allocVector(VECSXP, 2));

    SET_VECTOR_ELT(ret, 0, allocMatrix(REALSXP, nnode-ntip, 2));
    SET_VECTOR_ELT(ret, 1, allocVector(REALSXP, nnode-ntip));

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        i = phy_node_index(node);
        j = phy_node_index(phy_node_lfdesc(node));
        k = phy_node_index(phy_node_rtdesc(node));
        REAL(VECTOR_ELT(ret, 0))[i-ntip] = u[i];
        REAL(VECTOR_ELT(ret, 0))[(i-ntip) + (nnode-ntip)] = vx[j] + vx[k];
        REAL(VECTOR_ELT(ret, 1))[i-ntip] = x[i];
    }

    setAttrib(VECTOR_ELT(ret, 0), install("rate"), ScalarReal(bm.su / bm.n));

    UNPROTECT(1);
    return ret;
}




