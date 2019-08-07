#################################################################################
#### Normally, we will perform the blast searching with online version and then download the table.
#### But the table normally could not provide enough data which only contains a id.
#### this script is prepared to process it.
####
#################################################################################

import pandas as pd
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "l0404th@gmail.com"


def get_info(info):
    title = info[0]["Title"]
    seqID = info[0]["Id"]
    accession_id = info[0]["AccessionVersion"]
    createDate = info[0]["CreateDate"]
    updateDate = info[0]["UpdateDate"]
    taxid = info[0]["TaxId"].real
    sci_name = Entrez.read(Entrez.esummary(id=taxid, db="taxonomy"))[0]["ScientificName"]
    length = int(info[0].get("Length", -1))
    collect_dict = dict(title=title, seqID=seqID, accession=accession_id, CreateDate=createDate, UpdateDate=updateDate, taxid=taxid, sci_name=sci_name, length=length)
    return collect_dict


def get_parsed_df(file_path, db):
    # it depends on the mode of blast it used
    # here I use blastx which use nucleotide to search protein database.
    collect_null_ids = []
    new_df = pd.DataFrame()

    df = pd.read_csv(file_path, index_col=0, header=None)
    df = df.drop_duplicates(subset=1)
    for rid, row in tqdm(df.iterrows(),
                         total=df.shape[0]):
        id = row.values[0]
        try:
            infos = Entrez.read(Entrez.esummary(id=id, db=db))
        except RuntimeError as e:
            print(id, e.args)
            collect_null_ids.append(id)
            continue
        # [DictElement({'Item': [], 'Id': '498530462', 'Caption': 'AGL46827', 'Title': 'photosystem II D1 core reaction center protein, partial [Synechococcus T7-like virus S-TIP37]', 'Extra': 'gi|498530462|gb|AGL46827.1|[498530462]', 'Gi': IntegerElement(498530462, attributes={}), 'CreateDate': '2013/05/18', 'UpdateDate': '2013/05/18', 'Flags': IntegerElement(0, attributes={}), 'TaxId': IntegerElement(1332145, attributes={}), 'Length': IntegerElement(227, attributes={}), 'Status': 'live', 'ReplacedBy': '', 'Comment': '  ', 'AccessionVersion': 'AGL46827.1'}, attributes={})]
        if infos:
            collect_dict = get_info(infos)
            new_df = new_df.append(pd.DataFrame.from_dict(collect_dict, orient="index").T)
    return new_df, collect_null_ids


def get_nuccore_df(new_df: pd.DataFrame):
    "from protein info dataframe to fetch corresponding nucc dataframe"
    nucc_df = pd.DataFrame(columns=["accession"])
    for rid, row in tqdm(new_df.iterrows(), total=new_df.shape[0]):
        gid = row["accession"]
        result = Entrez.read(Entrez.elink(id=gid, dbfrom='protein', db="nuccore"))
        assert len(result) == 1
        result = result[0]["LinkSetDb"]
        if not len(result) == 1:
            result = [linkid for diff_link in result for linkid in diff_link["Link"]]
        else:
            result = result[0]["Link"]
        nids = list(set([_["Id"] for _ in result]))
        for nid in nids:
            # if nid in nucc_df.loc[:, "accession"]:
            #     continue
            nucc_info = Entrez.read(Entrez.esummary(id=nid, db="nuccore"))
            sub_dict = get_info(nucc_info)
            sub_dict["protein id"] = gid
            sub_df = pd.DataFrame.from_dict(sub_dict, orient="index").T
            nucc_df = nucc_df.append(sub_df, sort=False)
    return nucc_df


if __name__ == '__main__':
    file_path = "file:///D:/Desktop/psbA_nr_search/MCHYU3F1015-Alignment-HitTable.csv"
    db = "protein"
    protein_df, null_ids = get_parsed_df(file_path, db)
    protein_df.to_csv("D://Desktop//psbA_nr_search//parsed_result.csv", index=False)
