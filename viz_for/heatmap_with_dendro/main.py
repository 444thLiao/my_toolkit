import plotly.plotly as py
import plotly.graph_objs as go
import plotly.figure_factory as ff
import plotly
import numpy as np
from scipy.spatial.distance import pdist, squareform
import pandas as pd

data = pd.read_csv('/home/liaoth/Desktop/tmp_heatmap.data',sep='\t',index_col=0)
# data = data.iloc[:100,:]

data_array = data.T.values

# Initialize figure by creating upper dendrogram
figure = ff.create_dendrogram(data_array, orientation='bottom', labels=data.columns)
for i in range(len(figure['data'])):
    figure['data'][i]['yaxis'] = 'y2'

# Create Side Dendrogram
dendro_side = ff.create_dendrogram(data.values, orientation='right',labels=data.index)
for i in range(len(dendro_side['data'])):
    dendro_side['data'][i]['xaxis'] = 'x2'

# Add Side Dendrogram Data to Figure
for _data in dendro_side['data']:
    figure.add_trace(_data)

# Create Heatmap
heatmap = [
    go.Heatmap(
        x = data.columns,
        y = data.index,
        z = data.values,
    )
]

heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

# Add Heatmap Data to Figure
for data in heatmap:
    figure.add_trace(data)

# Edit Layout
figure['layout'].update({
                         'showlegend':False, 'hovermode': 'closest',
                         })
# Edit xaxis
figure['layout']['xaxis'].update({'domain': [.15, 1],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'ticks':""})
# Edit xaxis2
figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})

# Edit yaxis
figure['layout']['yaxis'].update({'domain': [0, .85],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})
# Edit yaxis2
figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})

# Plot!
plotly.offline.plot(figure, filename='dendrogram_with_heatmap')