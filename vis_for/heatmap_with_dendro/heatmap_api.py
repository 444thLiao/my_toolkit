import plotly
import plotly.figure_factory as ff
import plotly.graph_objs as go


def create_heatmap(left_dist,  main_matrix,up_dist=None, accessory_matrix=None,height=1500,width=None):
    fig = plotly.tools.make_subplots(2, 3,
                                     specs=[[None, None, {}],
                                            [{}, {}, {}]],
                                     shared_yaxes=True,
                                     shared_xaxes=True)
    # This is the format of your plot grid:
    #     (empty)          (empty)      [ (1,3) x3,y1 ]
    # [ (2,1) x1,y2 ]  [ (2,2) x2,y2 ]  [ (2,3) x3,y2 ]
    sub_df = main_matrix
    # create up side dendrogram
    if up_dist is not None:
        up_dendro = ff.create_dendrogram(up_dist.values,
                                         orientation='bottom',
                                         labels=up_dist.index)

        sub_df = main_matrix.reindex(columns=up_dendro.layout.xaxis.ticktext)
        # coordinate to up dend
        for i in up_dendro.data:
            i.showlegend = False
            fig.append_trace(i, 1, 3)
    # Create left Side Dendrogram
    side_dendro = ff.create_dendrogram(left_dist.values,
                                       orientation='right',
                                       labels=left_dist.index)
    for i in side_dendro.data:
        i.showlegend = False
        fig.append_trace(i, 2, 1)
    sub_df = sub_df.reindex(side_dendro.layout.yaxis.ticktext)
    # coordinate to left dend
    ############################################################
    # draw (2,3) heatmap
    heatmap = go.Heatmap(
        x=sub_df.index,
        y=sub_df.columns,
        z=sub_df.values,
        colorscale='Earth',
        reversescale=True,
        showscale=False
    )
    if up_dist is not None:
        heatmap['x'] = up_dendro.layout.xaxis.tickvals
    heatmap['y'] = side_dendro.layout.yaxis.tickvals
    fig.append_trace(heatmap, 2, 3)

    ############################################################
    if accessory_matrix:
        _sub_df = accessory_matrix.reindex(side_dendro.layout.yaxis.ticktext)

        heatmap = go.Heatmap(
            x=_sub_df.columns,
            y=_sub_df.index,
            z=_sub_df.values,
            text=_sub_df.values,
            hoverinfo='all',
            colorscale='Earth',
            reversescale=True,
            showscale=False
        )
        fig.append_trace(heatmap, 2, 2)
    ############################################################

    if accessory_matrix is not None:
        fig.layout.xaxis1.domain = [0, 0.05]
        fig.layout.xaxis2.domain = [0.05, 0.1]
        fig.layout.xaxis3.domain = [0.1, 1]
    else:
        fig.layout.xaxis1.domain = [0, 0.17]
        fig.layout.xaxis2.domain = [0.0, 0.0]
        fig.layout.xaxis3.domain = [0.2, 1]
    if up_dist is not None:
        fig.layout.yaxis1.domain = [0.9, 1]
        fig.layout.yaxis2.domain = [0.0, 0.9]
    else:
        fig.layout.yaxis2.domain = [0.0, 0]
        fig.layout.yaxis1.domain = [0.0, 1]
    fig.layout.yaxis2.ticktext = side_dendro.layout.yaxis.ticktext
    fig.layout.yaxis2.tickvals = side_dendro.layout.yaxis.tickvals
    if up_dist is not None:
        fig.layout.xaxis3.ticktext = up_dendro.layout.xaxis.ticktext
        fig.layout.xaxis3.tickvals = up_dendro.layout.xaxis.tickvals
    # fig.layout.xaxis.anchor = 'x2'
    # fig.layout.margin.b = 250
    fig.layout.width = main_matrix.shape[1] * 13 if width is None else width
    fig.layout.height = height
    fig.layout.xaxis3.tickangle = 30
    fig.layout.xaxis1.zeroline = fig.layout.xaxis2.zeroline = fig.layout.xaxis3.zeroline = False
    fig.layout.yaxis1.zeroline = fig.layout.yaxis2.zeroline = False
    fig.layout.xaxis1.showgrid = False
    fig.layout.yaxis2.showgrid = False
    fig.layout.xaxis3.showgrid = False
    fig.layout.yaxis1.showgrid = False
    fig.layout.xaxis1.showticklabels = False
    fig.layout.xaxis3.showticklabels = False
    fig.layout.xaxis1.showspikes = False
    fig.layout.xaxis3.showspikes = False
    fig.layout.xaxis2.showspikes = False
    fig.layout.yaxis1.showticklabels = False
    return fig


def main(idir, odir):
    pass


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    main(args.indir, args.outdir)
