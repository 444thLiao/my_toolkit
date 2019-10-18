import click
from bioservices.kegg import KEGG
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
from os.path import *
import os
def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d: len(iter) + 1])
    return n_iter

def retry_get(id,max_try=10):
    count_ = 0
    info_str = None
    while count_ <= max_try:
        info_str = kegg.get(id)
        # if failed, it will return 400 or 403. still an int
        if isinstance(info_str, str):
            break
        count_ += 1
    return info_str
 
def get_KO_info(ID, max_try=10):
    # ID could use + to concate multiple ko
    info_str = 0
    
    return_dict = {}
    info_str = retry_get(ID,max_try)
    info_dict_list = [kegg.parse('ENTRY ' + each_str)
                    for each_str in info_str.split('\nENTRY ')
                    if each_str]
    # make the first one entry startwith ENTRY instead of original locus.
    for info_dict in info_dict_list:
        if not isinstance(info_dict, dict):
            #print(info_dict)
            continue
        entry = info_dict.get('ENTRY', 'unknown').split(' ')[0]
        # entry name, normally is the input ko
        # sometimes, it maybe no entry.... =^=
        if entry.startswith('ENTRY'):
            entry = [_
                    for _ in info_dict.get('ENTRY', 'unknown').split(' ')
                    if _][1]
            # remove ENTRY to get after name
            
        _cache = [ori for ori, _ in zip(ID.split('+'),
                                        ID.lower().split('+'))
                if entry.lower() in _]
        # map to original ID
        if len(_cache) == 0:
            print(entry, ID)
            # failed to get it name
            continue
        ko = _cache[0]
        gene_name = ';'.join(info_dict.get('NAME', ['']))
        definition = info_dict.get('DEFINITION', '')
        pathway_dict = info_dict.get('PATHWAY',{})
        reference_t = ''
        if "REFERENCE" in info_dict:
            reference_t = ';'.join([str(_dict.get('TITLE', ''))
                                    for _dict in info_dict.get('REFERENCE', {})])

        return_dict[ko] = dict(gene_name=gene_name,
                            definition=definition,
                            reference_t=reference_t,
                            pathway=pathway_dict)
        
    return return_dict


def main(in_file,ko_col_num = 0):
    if in_file.endswith('.xlsx'):
        df = pd.read_excel(in_file,index_col=None)
#    elif in_file.endswith('.csv'):
 #       df = pd.read_csv(in_file,index_col=None)
    else:
        df = pd.read_csv(in_file,index_col=None)
    ko_col = list(df.iloc[:,ko_col_num].unique())
    pack10_up = batch_iter(ko_col, 10)
    
    all_ko_dict = {}
    tqdm.write('Start to fetch relative information of K number.')
    for ko_list in tqdm(pack10_up):
        ko_dict = get_KO_info('+'.join(ko_list))
        
        all_ko_dict.update(ko_dict)
    
    df.loc[:,'num_pathway'] = 0
    df.loc[:,'pathway'] = ''
    df.loc[:,'pathway_description'] = ''
    # add new columns
    
    new_df = df.set_index(df.columns[ko_col_num])
    all_p_names = []
    for ko in new_df.index:
        p_dict = all_ko_dict.get(ko,{}).get('pathway',{})
        ref_t = all_ko_dict.get(ko,{}).get('reference_t','')
        names = list(sorted(p_dict.keys()))
        all_p_names+=names
        new_df.loc[ko,'num_pathway'] = len(p_dict)
        new_df.loc[ko,'pathway'] = ';'.join(names)
        new_df.loc[ko,'pathway_description'] = ';'.join([p_dict.get(_) for _ in names])
        new_df.loc[ko,'title of refence'] = ref_t
        
    pathway_df_dict = defaultdict(dict)
    pathway2ko = defaultdict(list)
    pathway2des = {}
    for ko,info_dict in all_ko_dict.items():
        pathways = info_dict.get('pathway',{}).keys()
        for p in pathways:
            pathway2ko[p].append(ko)
            pathway2des[p] = info_dict.get('pathway',{}).get(p,'')
            
    sub_df = new_df.drop_duplicates()
    for p,ko_list in pathway2ko.items():
        pathway_df_dict[p]['name'] = p
        pathway_df_dict[p]['description'] = pathway2des[p]
        pathway_df_dict[p]['appear_ko_list'] = '+'.join(ko_list)
        pathway_df_dict[p]['Each appear times'] = '+'.join(list(map(lambda x: str(sub_df.loc[ko,'appear times']),ko_list)))
    pathway_df = pd.DataFrame.from_dict(pathway_df_dict,orient='index')
    # old_columns = pathway_df.columns
    # new_columns = ['Each appear times'] + list(old_columns)[:-1]
    #pathway_df.reindex(columns= new_columns)
    return new_df,pathway_df

@click.command()
@click.option('-i','in_file',help='input file.')
@click.option('-o','output_dir',help='output directory.')
@click.option('-k','ko_col_num',help='the number of columns of ko, normally is the first columns',default=1)
@click.option('-sub','substite_org_name',help='accept a file to replace the org name',default=None)
def cli(in_file,output_dir,ko_col_num,substite_org_name):
    ko_col_num = int(ko_col_num) -1
    if not exists(output_dir):
        os.makedirs(output_dir)
    new_df,pathway_df = main(in_file,ko_col_num)
    old_columns = list(new_df.columns)
    old_columns.pop('appear times')
    new_columns = ['appear times'] + old_columns
    new_df = new_df.reindex(columns=new_columns)
    if substite_org_name:
        rows = open(substite_org_name,'r').read().split('\n')
        org2new_name = dict([row.split('\t') for row in rows])
        new_df.loc[:,'name'] = [org2new_name.get(_,'Unknown') for _ in new_df.loc[:,'name']]
    new_df.to_csv(join(output_dir,basename(in_file)),index=1,index_label='k number')
    pathway_df.to_csv(join(output_dir,'pathway_info.csv'),index=0)

if __name__ == "__main__":
    kegg = KEGG()
    cli()