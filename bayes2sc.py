def rowcol2pixel(ad_input, ad_referece):
    import numpy as np
    from sklearn.linear_model import LinearRegression
    
    cols = []
    rows = []
    for s in ad_referece.obs.index:
        _, col, row = s.split('_')
        cols.append(int(col))
        rows.append(int(row))

    cols = np.array(cols)
    rows = np.array(rows)
    
    label_x, label_y = ad_referece.obsm["spatial"][:, 0], ad_referece.obsm["spatial"][:, 1]

    input_cols = ad_input.obs["col"].values
    input_rows = ad_input.obs["row"].values
    
    lr = LinearRegression()
    lr.fit(rows.reshape(-1, 1), label_y)
    pred_y = lr.predict(input_rows.reshape(-1, 1))

    lr = LinearRegression()
    lr.fit(cols.reshape(-1, 1), label_x)
    pred_x = lr.predict(input_cols.reshape(-1, 1))

    ad_input.obsm['spatial'] = np.column_stack((pred_x, pred_y))
    
def adata_HDspot_to_cell(adata_spot, node_feat, margin=220):
    """ Assigning ST of segmented cells by the nearest spots within the margin.
    
    Parameters
    ----------
    adata_spot: AnnData
        The AnnData object containing the spot data.
    node_feat: pd.DataFrame
        The dataframe containing the cell features.
    margin: float
        The distance threshold to find the closest spots. 

    Returns
    -------
    adata_cell: AnnData
        The AnnData object containing the cell data.
    spot_to_cell: list
        The list of spot indices that are mapped to cells.
    """
    import pandas as pd
    from sklearn.neighbors import NearestNeighbors
    
    print(
        "The first two columns in the node_feat DataFrame need to be consistent with the spatial coordinates from obsm['spatial']."
    )

    spot_pos = pd.DataFrame(
        adata_spot.obsm['spatial'],
        columns=['x', 'y'],
        index=adata_spot.obs.index
    )

    nbrs = NearestNeighbors(n_neighbors=1).fit(spot_pos[['x', 'y']].values)
    distances, indices = nbrs.kneighbors(node_feat.iloc[:, :2].values)

    if margin > 0:
        selected = (distances < margin)
        selected = selected[:, 0]

        # those are spot indices
        indices = indices[selected]
        node_feat = node_feat[selected]

    nearest_spot = spot_pos.iloc[indices.T[0], :]
    spot_to_cell = nearest_spot.index.values

    adata_cell = adata_spot[spot_to_cell]
    adata_cell.obs['spot_barcodes'] = adata_cell.obs.index.tolist()
    adata_cell.obsm['spatial'] = node_feat.iloc[:, :2].values
    

    df = pd.concat(
        [
            adata_cell.obs.set_index(node_feat.index),
            node_feat
        ],
        axis=1
    )

    #adata_cell.obs = df.set_index(adata_cell.obs.index)
    adata_cell.obs = df

    return adata_cell, spot_to_cell
