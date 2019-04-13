"""
fastest way to get taxonomy id from a higher level such as mycobacterium.
"""

import requests
import pandas as pd
from bs4 import BeautifulSoup


base_url = "https://www.ncbi.nlm.nih.gov/"


def get_assembly(tid,page=200):
    mission = "assembly/?term=txid%s[Organism:exp]&" % tid
    cond = ''
    url = base_url + mission + cond
    r =requests.get(url.encode())
    content = r.content

    parser_info = BeautifulSoup(content,'lxml')

def get_info_taxonomy(tid,lvl=5,records=("gcassembly",
                                         "biosample",
                                         "genome",
                                         "has_linkout")):
    mission = "Taxonomy/Browser/wwwtax.cgi?"
    cond = 'id=%s&lvl=%s' % (tid,lvl)
    cond += '&p='.join(['']+records)
    url = base_url + mission + cond
    r = requests.get(url)
    content = r.content

    parser_info = BeautifulSoup(content,'lxml')
    return parser_info

def parse_leaf(info):
    list_info = info.findAll('a')
    result = {'name_href':list_info[0].attrs['href'],
              'tid':''}
    if result['name_href']:
        _cache = [_ for _ in result['name_href'].split('&') if _.startswith('id=')]
        if _cache:
            result['tid'] = _cache[0].split('id=')[1]
    for i in list_info[1:]:
        title = i.attrs['title']
        href = i.attrs['href']
        result[title] = i.text
        result[title+'_href'] = href
    return list_info[0].text,result


def parse_tax_tree(parser_info):
    first_record = parser_info.find('li',attrs={'type':'circle'})
    # todo:
    leaf_nodes = parser_info.findAll(name='li', attrs={'type': 'square'})
    sto_dict = {}
    for leaf in leaf_nodes:
        name,leaf_dict = parse_leaf(leaf)
        sto_dict[name] = leaf_dict
    return sto_dict

if __name__ == '__main__':
    tid = '77643' # Mycobacterium tuberculosis complex
    sto_dict = get_info_taxonomy(tid)
    sto_df = pd.DataFrame.from_dict(sto_dict, orient='index')

