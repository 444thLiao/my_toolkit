from collections import Counter

import matplotlib.pyplot as plt
import numpy as np


# matplotlib version
def sunburst(nodes, total=np.pi * 2, offset=0, level=0, ax=None, fig=None, fontsize=12):
    fig = fig or plt.figure(figsize=(10, 10))
    ax = ax or plt.subplot(111, projection='polar')

    if level == 0 and len(nodes) == 1:
        label, value, subnodes = nodes[0]
        ax.bar([0], [0.5], [np.pi * 2], alpha=0.5)
        ax.text(0, 0, label, ha='center', va='center', fontsize=fontsize)
        sunburst(subnodes, total=value, level=level + 1, ax=ax, fig=fig, fontsize=fontsize)
    elif nodes:
        d = np.pi * 2 / total
        labels = []
        widths = []
        local_offset = offset
        for label, value, subnodes in nodes:
            labels.append(label)
            widths.append(value * d)
            sunburst(subnodes, total=total, offset=local_offset,
                     level=level + 1, ax=ax, fig=fig, fontsize=fontsize)
            local_offset += value
        values = np.cumsum([offset * d] + widths[:-1])
        heights = [1] * len(nodes)
        bottoms = np.zeros(len(nodes)) + level - 0.5
        rects = ax.bar(values, heights, widths, bottoms, linewidth=1,
                       edgecolor='white', align='edge', alpha=0.5)
        for rect, label in zip(rects, labels):
            x = rect.get_x() + rect.get_width() / 2
            y = rect.get_y() + rect.get_height() / 2
            rotation = (90 + (360 - np.degrees(x) % 180)) % 360
            ax.text(x, y, label, rotation=rotation, ha='center', va='center', fontsize=fontsize)

    if level == 0:
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location('N')
        ax.set_axis_off()


def recursive_parse(sub_df, level):
    # parse a table data into suitable hierarchy structure.
    tax_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
    data = []
    _counter = Counter(sub_df.loc[:, level])

    if tax_levels.index(level) != len(tax_levels) - 1:
        for _c, v in _counter.items():
            data.append((_c, v, recursive_parse(sub_df.loc[sub_df.loc[:, level] == _c, :],
                                                tax_levels[tax_levels.index(level) + 1])))
    else:
        for _c, v in _counter.items():
            data.append((_c, v, []))
    return data


if __name__ == '__main__':
    import pandas as pd
    data = pd.read_csv("/home/liaoth/data2/project/NR-Pfam/cluster2domains/NasN info new.tab",sep='\t',index_col=0)
    data = data.loc[data.loc[:,'type']=='NasN',['Phylum', 'Class', 'Order', 'Family', 'Genus']]
    fig = plt.figure(figsize=(10, 10))
    sunburst([('Bacteria', data.shape[0], recursive_parse(data, 'Phylum'))], fig=fig, fontsize=5)
    plt.savefig('/home/liaoth/data2/project/NR-Pfam/cluster2domains/distribution.pdf', format='pdf')
    plt.savefig('/home/liaoth/data2/project/NR-Pfam/cluster2domains/distribution.png', format='png',dpi=300)

#
# data = [
#     ('/', 100, [
#         ('home', 70, [
#             ('Images', 40, []),
#             ('Videos', 20, []),
#             ('Documents', 5, []),
#         ]),
#         ('usr', 15, [
#             ('src', 6, [
#                 ('linux-headers', 4, []),
#                 ('virtualbox', 1, []),
#
#             ]),
#             ('lib', 4, []),
#             ('share', 2, []),
#             ('bin', 1, []),
#             ('local', 1, []),
#             ('include', 1, []),
#         ]),
#     ]),
# ]
#
# sunburst(data)
