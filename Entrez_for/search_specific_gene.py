#################################################################################
#### for searching a gene sequence from published paper
#### using `pubmed` and `Entrez` to get required sequence
#### 2019.08.04    For Luo Lab
#################################################################################
from Bio import Entrez

Entrez.email = "l0404th@gmail.com"
from tqdm import tqdm
import pandas as pd


def get_info(id):
    info = Entrez.read(Entrez.esummary(id=id, db='nuccore'))
    assert len(info) == 1
    title = info[0]["Title"]
    giid = info[0]["Gi"].real
    accession_id = info[0]["AccessionVersion"]
    createDate = info[0]["CreateDate"]
    updateDate = info[0]["UpdateDate"]
    taxid = info[0]["TaxId"].real
    length = info[0]["Length"]
    collect_dict = dict(title=title, gi=giid, accession=accession_id, CreateDate=createDate, UpdateDate=updateDate, taxid=taxid, length=length)
    return collect_dict


def get_paper_info(id):
    info = Entrez.read(Entrez.esummary(id=id, db='pubmed'))
    assert len(info) == 1
    paper_title = info[0]["Title"]
    paper_id = info[0]["Id"]
    PubDate = info[0]["PubDate"]
    collect_dict = dict(paper_title=paper_title, paper_id=paper_id, PubDate=PubDate)
    return collect_dict


def filter_null():
    # removing some record with self adjustment.
    pass


# input word need to search, please check it by manually searching.
term = 'psbA cyanophage'
gene_relative = Entrez.read(Entrez.esearch(term=term, db='pubmed'))

# init a dataframe for stodge fetched metadata of sequences.
collect_df = pd.DataFrame(columns=['title',
                                   'gi',
                                   'accession',
                                   'CreateDate',
                                   'UpdateDate',
                                   'taxid',
                                   'length',
                                   'paper_title',
                                   'paper_id',
                                   'PubDate'])
for pid in tqdm(gene_relative["IdList"] + ["18586962"]):
    get_nucc = Entrez.read(Entrez.elink(dbfrom='pubmed', db='nucleotide', id=pid))
    get_ids = set([link["Id"] for _ in get_nucc for nid in _["LinkSetDb"] for link in nid["Link"]])
    if get_ids:
        print(pid, len(get_ids))
        for each_id in get_ids:
            collect_df = collect_df.append(pd.DataFrame.from_dict(get_info(each_id), orient="index").T,sort=False)

# output this dataframe, and manually check it again
collect_df.to_csv(".\psbA_relative.csv")