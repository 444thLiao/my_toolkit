import os
from os.path import dirname,abspath,join
import pandas as pd
root_path = os.path.abspath(__file__)
templates_dir = join(dirname(root_path),"templates")

templates_file = join(templates_dir,'dataset_color_strip_template.txt')


"""
LEGEND_TITLE,{legend_title}
LEGEND_SHAPES,{legend_shape}
LEGEND_COLORS,{legend_colors}
LEGEND_LABELS,{legend_labels}
"""


############################################################
# for color
def get_annotated_text(data_file,
                       template_file,
                       id_column,
                       annotate_column,
                       output_file,
                       color_scheme="Set1",
                       dataset_label='label1',
                       ):
    template_t = open(template_file).read()
    import seaborn as sns
    colors = sns.color_palette(color_scheme, 10).as_hex()
    nuc_df = pd.read_csv(data_file, index_col=None)
    labels = list(map(str, set(nuc_df.loc[:, annotate_column])))
    labels = [_ for _ in labels if _.lower() != 'nan']
    l2c = dict(zip(labels, colors))
    legend_title = 'which %s' % annotate_column
    legend_shape = ','.join(['1'] * len(l2c))
    legend_colors = ','.join([l2c[_]
                              for _ in labels])
    legend_labels = ','.join(list(labels))
    legend_text = f"""
    LEGEND_TITLE,{legend_title}
    LEGEND_SHAPES,{legend_shape}
    LEGEND_COLORS,{legend_colors}
    LEGEND_LABELS,{legend_labels}
    """
    annotate_text = ''
    for rid, row in nuc_df.iterrows():
        if str(row[annotate_column]).lower() == "nan":
            continue
        annotate_text += ','.join(map(str, [row[id_column],
                                            l2c[row[annotate_column]]])) + '\n'
    with open(output_file, 'w') as f1:
        f1.write(template_t.format(legend_text=legend_text,
                                   dataset_label=dataset_label) + annotate_text)


get_annotated_text('concat_all/total_nuc.csv',
                   'concat_all/colorize/dataset_color_strip_template.txt',
                   "seqID",
                   "family",
                   "concat_all/colorize/color_family_nuc.txt",
                   dataset_label='taxon family with nuc')

get_annotated_text('concat_all/total_nuc.csv',
                   'concat_all/colorize/dataset_color_strip_template.txt',
                   "seqID",
                   "host",
                   "concat_all/colorize/color_host_nuc.txt",
                   dataset_label='host with nuc')

get_annotated_text('concat_all/total_pro.csv',
                   'concat_all/colorize/dataset_color_strip_template.txt',
                   "accession",
                   "family",
                   "concat_all/colorize/color_family_pro.txt",
                   dataset_label='taxon family with nuc')

get_annotated_text('concat_all/total_pro.csv',
                   'concat_all/colorize/dataset_color_strip_template.txt',
                   "accession",
                   "host",
                   "concat_all/colorize/color_host_pro.txt",
                   dataset_label='host with nuc')


############################################################
def relabels_text(data_file,
                  template_file,
                  id_column,
                  annotate_column,
                  output_file,
                  ):
    nuc_df = pd.read_csv(data_file, index_col=None)
    template_t = open(template_file).read()
    annotate_text = ''
    for rid, row in nuc_df.iterrows():
        if str(row[annotate_column]).lower() == "nan":
            continue
        annotate_text += ','.join(map(str, [row[id_column],
                                            row[annotate_column]])) + '\n'
    with open(output_file, 'w') as f1:
        f1.write(template_t + annotate_text)


relabels_text('concat_all/total_pro.csv',
              'concat_all/colorize/labels_template.txt',
              "accession",
              "sci_name",
              "concat_all/colorize/labels_name_pro.txt")

relabels_text('concat_all/total_nuc.csv',
              'concat_all/colorize/labels_template.txt',
              "seqID",
              "sci_name",
              "concat_all/colorize/labels_name_nuc.txt")

if __name__ == '__main__':
    pass