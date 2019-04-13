import webbrowser
import json
import os
from collections import Counter
import io
dir_path = os.path.dirname(__file__)
def plot_net(graph, filename):
    with open(os.path.join(dir_path,"plot.js"), "r") as f:
        data_func = f.read()
    with open(filename, "w") as f:
        f.write(''.join([
            '<html>',
            '<head><meta charset="utf-8" /></head>',
            "<script src='https://d3js.org/d3.v5.min.js'></script>",
            '<body>',
            '<div style="text-align: center">',
            '<svg></svg>',
            '</div>'
            "<script>",
            data_func,
            "read_data({});".format(graph),
            "</script>"
        ]))
    webbrowser.open_new_tab(filename)

def recursive_parse(sub_df, level):
    # parse a table data into suitable hierarchy structure.
    tax_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']

    data = {"name": ""}

    _counter = Counter(sub_df.loc[:, level])

    if tax_levels.index(level) != len(tax_levels) - 1:
        data["children"] = []
        for _c, v in _counter.items():
            if type(_c) == str:
                data["name"] = _c
                data["children"].append(recursive_parse(sub_df.loc[sub_df.loc[:, level] == _c, :],
                                                tax_levels[tax_levels.index(level) + 1]))
    else:

        for _c, v in _counter.items():
            data["name"] = _c
            data["value"] = v
    return data
def parse(sub_df,level):
    data = []
    for phy in set(sub_df.loc[:,level]):
        _sub = sub_df.loc[sub_df.loc[:,level]==phy,:]
        data.append(recursive_parse(_sub,level))
    return {"name":"Bacteria","children":data}

if __name__ == "__main__":
    # 在python中打开一个json文件
    import pandas as pd
    data = pd.read_csv("/home/liaoth/data2/project/NR-Pfam/cluster2domains/NasN info new.tab",sep='\t',index_col=0)
    data = data.loc[data.loc[:,'type']=='NasN',['Phylum', 'Class', 'Order', 'Family', 'Genus']]
    # 绘图并打开一个html文件
    plot_net(parse(data,'Phylum'), "net.html")
