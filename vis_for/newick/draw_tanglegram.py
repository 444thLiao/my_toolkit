"""
This script is mainly written for tanglegram which connect leafs from two dendrograms.
Further, it also implement some useful transferring function from iTOL.
    for adding color to the leaves of the dendrogram. you could use the iTOL file.. or something simpler like removing the description

"""

from collections import defaultdict

import click
import os
import plotly
import plotly.graph_objs as go
from tqdm import tqdm

from Tree_for.api_IO import read_tree
from Tree_for.plotly_draw_newick_tree import main as get_plotly_data_from_newick


def get_preferred_scale(newick1, newick2):
    t1 = read_tree(newick1, format='auto')
    t2 = read_tree(newick2, format='auto')
    num_1 = len(t1.get_leaf_names())
    num_2 = len(t2.get_leaf_names())

    yscale = num_1 / num_2
    return 1 / yscale


def parse_color_scheme_files(file, extra_set=False):
    lines = open(file).read().split('\n')
    lines = [_ for _ in lines if _]
    field_labels = [_ for _ in lines if _.startswith("FIELD_LABELS")]
    field_colors = [_ for _ in lines if _.startswith("FIELD_COLORS")]
    name2color = {}
    sep_indicator = [_ for _ in lines if _.startswith("SEPARATOR")][0]
    sep_indicator = sep_indicator.split(' ')[-1]
    s2s = {"TAB": '\t', "COMMA": ',', "SPACE": ' '}
    sep = s2s.get(sep_indicator, ',')

    if field_labels and field_colors:
        colors = field_colors[0].split('\t')[1:]
        anno_names = field_labels[0].split('\t')[1:]
        _name2data = [_ for _ in lines[lines.index("DATA") + 1:] if _ and not _.startswith('#')]
        for line in _name2data:
            line.strip('\n')
            vs = line.split('\t')
            name = vs[0]
            color = [c for c, _ in zip(colors, vs[1:]) if str(_) != '0']
            if color:
                name2color[name] = color[0]  # code. would overlap.. be careful..
    else:
        lines = [_ for _ in lines[lines.index("DATA") + 1:] if _ and not _.startswith('#')]
        for line in lines:
            name = line.split(sep)[0]
            if "range" in line:
                color = line.split(sep)[2]
            else:
                color = line.split(sep)[1]
            name2color[name] = color

    if str(extra_set) == 'rename':
        new_name2color = {k.split('_')[-1].replace('.', 'v'): v
                          for k, v in name2color.items()}
        name2color = new_name2color.copy()
    return name2color


def main(newick1, newick2,
         color_file1, color_file2,
         l_legnth='max', sep='_', extra_set=False):
    yscale = get_preferred_scale(newick1, newick2)
    fig = plotly.subplots.make_subplots(rows=1, cols=3, shared_yaxes=True)

    # get dendrogram parts
    tqdm.write('drawing dendrograms')
    datas, labels, _, labels_draw_text, labels_x, labels_y = get_plotly_data_from_newick(newick1,
                                                                                         fixed_length=l_legnth,
                                                                                         yscale=yscale)
    datas2, labels2, _, labels_draw_text2, labels_x2, labels_y2 = get_plotly_data_from_newick(newick2,
                                                                                              fixed_length=l_legnth)
    # add colors or something else into above datas
    l2color = {_: '#000000' for _ in labels_draw_text}
    r2color = {_: '#000000' for _ in labels_draw_text2}
    l2color.update(parse_color_scheme_files(color_file1, extra_set=False))
    r2color.update(parse_color_scheme_files(color_file2, extra_set=extra_set))
    # add color into generated trace
    tqdm.write('adding color')
    name2trace = {_['name']: _ for _ in datas if _['name'] is not None}
    for l, color in l2color.items():
        if color != '#00000' and l in labels_draw_text:
            trace = name2trace[l]
            trace['line']['color'] = color
    name2trace = {_['name']: _ for _ in datas2 if _['name'] is not None}
    for r, color in r2color.items():
        if color != '#00000' and r in labels_draw_text2:
            trace = name2trace[r]
            trace['line']['color'] = color

    # add dendrogram parts into figure.
    # data1/newick1 would be the left, data2/newick2 would be the right and the leaves of it will point to left
    # for _ in datas:
    fig.add_traces(datas, rows=[1] * len(datas), cols=[1] * len(datas))
    # for _ in datas2:
    fig.add_traces(datas2, rows=[1] * len(datas2), cols=[3] * len(datas2))

    # init data of middle part
    # get the y-coordinate information from below. put them into two dict.
    left_data = dict(zip(labels_draw_text, labels_y))
    right_data = dict(zip(labels_draw_text2, labels_y2))

    # get mapping relationship, default is from left to the right.. one to multi
    # so, the leaf names from the left tree should be the part of right, separate with `underline` or `space`
    l2r = defaultdict(list)
    for leaf1 in left_data.keys():
        leaf2s = [r for r in right_data.keys() if leaf1.split(sep)[0] == r]
        l2r[leaf1] = leaf2s

    # init the data from above mapping dict
    c2data = defaultdict(lambda: ([], []))
    for l, color in l2color.items():
        _xs = c2data[color][0]
        _ys = c2data[color][1]
        rs = l2r.get(l, [])
        l_y = left_data.get(l, 0)
        for r in rs:
            r_y = right_data[r]
            _xs += [0, 1, None]
            _ys += [l_y, r_y, None]

    for color, (_xs, _ys) in c2data.items():
        trace = go.Scatter(x=_xs,
                           y=_ys,
                           mode='lines',
                           line=dict(color=color, ),
                           hoverinfo='none',
                           showlegend=True)
        fig.add_traces([trace], rows=[1], cols=[2])
    fig.layout.xaxis3.autorange = 'reversed'
    fig.layout.xaxis.showticklabels = False
    fig.layout.xaxis2.showticklabels = False
    fig.layout.xaxis3.showticklabels = False
    fig.layout.xaxis.zeroline = False
    fig.layout.xaxis2.zeroline = False
    fig.layout.xaxis3.zeroline = False
    fig.layout.yaxis.zeroline = False
    fig.layout.yaxis.showticklabels = False

    return fig


@click.command()
@click.option("-newick1", help="first tree file, would be placed on the left")
@click.option("-newick2", help="second tree file, would be placed on the right")
@click.option("-output", "output_file", help="the path of html output ")
@click.option("-cf1", help="the file for color annotation to first tree and the links")
@click.option("-cf2", help="the file for color annotation to second tree ")
@click.option("-length", default="max", help="the length of leaves you want to extend to. normally all leaves will extend to identical length which similar to the longest one")
@click.option("-sep", default="_", help="the separator which used to . ")
@click.option("-extra_set", "extra_set", is_flag=True, default=False)
def cli(newick1, newick2, output_file, cf1, cf2, length, sep, extra_set):
    fig = main(newick1=newick1,
               newick2=newick2,
               color_file1=cf1,
               color_file2=cf2,
               l_legnth=length,
               sep=sep,
               extra_set=extra_set)

    fig.layout.width = 1400
    fig.layout.height = 3000

    if not exists(dirname(output_file)):
        os.makedirs(dirname(output_file))
    fig.write_html(output_file)


if __name__ == '__main__':
    from os.path import dirname, join, exists

    example_dir = join(dirname(dirname(dirname(__file__))), 'example', 'tanglegram')

    newick1 = join(example_dir, 'nxrA.newick')
    newick2 = join(example_dir, 'species.txt')

    color_file1 = join(example_dir, 'gene_annotation.txt')
    color_file2 = join(example_dir, 'phylum_annotate.txt')
    fig = main(newick1, newick2,
               color_file1, color_file2,
               l_legnth='max', sep='_', extra_set='rename')
    fig.layout.width = 1400
    fig.layout.height = 3000
    fig.write_html(join(example_dir, 'p.html'))
