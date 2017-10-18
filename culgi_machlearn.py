#!/usr/bin/python3


def remove_novar(df):
    """
    Function able to remove the 0 variance columns (constant values)
    """

    from sklearn.feature_selection import VarianceThreshold

    sel = VarianceThreshold()
    sel.fit_transform(df)

    # pick indices picked by selector
    indxs_sel = sel.get_support(indices=True)

    return (df[df.columns[indxs_sel]])


def univar_feat_sel(df, col="Solubility (cal/cc)", kvec=10, mutual=True):
    """
    Function able to find the kvec highest correlated columns with the one
    specified by col.  Two methods are available (f_reg & mutual_info)
    """

    from sklearn.feature_selection import SelectKBest

    # Select matrix and reference vector to compare with
    mat = df.drop(col, axis=1)
    ref = df[col]

    # select method
    if mutual:
        from sklearn.feature_selection import mutual_info_regression
        sel = SelectKBest(mutual_info_regression, k=kvec)
    else:
        from sklearn.feature_selection import f_regression
        sel = SelectKBest(f_regression, k=kvec)

    sel.fit_transform(mat, ref)

    # pick indices picked by selector
    indxs_sel = sel.get_support(indices=True)

    return (df[df.columns[indxs_sel]].join(ref))


def lassocv_feat_sel(df, col="Solubility (cal/cc)", nfeat=6):
    """
    Find the nfeat highest correlated columsn with the one specified by col.
    Similar to previous but now using regression with lasso regularization
    """

    from sklearn.feature_selection import SelectFromModel
    from sklearn.linear_model import LassoCV, Lasso

    # Select matrix and reference vector to compare with
    mat = df.drop(col, axis=1)
    ref = df[col]

    # We use the base estimator LassoCV since the L1 norm promotes sparsity of features.
    clf = LassoCV()

    # Set a minimum threshold of 0.25
    factor = 1.0
    sfm = SelectFromModel(clf, threshold=str(factor)+"*mean")
    sfm.fit(mat, ref)
    n_features = sfm.transform(mat).shape[1]
    print("Threshold: ",sfm.threshold_)

    # Reset the threshold till the number of features equals to nfeat
    # Note that the attribute can be set directly instead of repeatedly
    # fitting the metatransformer.
    while n_features > nfeat:

        factor += 0.2
        sfm.threshold = str(factor)+"*mean"
        x_transform = sfm.transform(mat)
        n_features = x_transform.shape[1]
        print("Threshold: ",sfm.threshold_)

    # pick indices picked by selector
    indxs_sel = sfm.get_support(indices=True)

    return (df[df.columns[indxs_sel]].join(ref))

