from Bio import Entrez
from bs4 import BeautifulSoup

Entrez.email = 'l0404th@gmail.com'


def org2genome_def(org):
    genome_def = ''
    genome_record = Entrez.read(Entrez.esearch(term=org, db='genome'))
    if genome_record['IdList']:
        gid = genome_record['IdList'][0]
        genome_info = Entrez.read(Entrez.esummary(id=gid, db='genome'))
        genome_def = genome_info[0]['DefLine']
    return genome_def
def org2bioproj_des(org):
    genome_def = ''
    ass_record = Entrez.read(Entrez.esearch(term=org, db='bioproject'))
    #all_des = []
    if ass_record['IdList']:
        for aid in ass_record['IdList']:
            ass_info = Entrez.read(Entrez.esummary(id=aid, db='bioproject'))
            #proj_des = ass_info['DocumentSummarySet']['DocumentSummary'][0].get('Project_Description','')
            proj_dict = ass_info['DocumentSummarySet']['DocumentSummary'][0]
            return proj_dict
            all_des.append(proj_des)
    return proj_dict
def org2bioproj_dict(org):
    genome_def = ''
    ass_record = Entrez.read(Entrez.esearch(term=org, db='bioproject'))
    #all_des = []
    if ass_record['IdList']:
        for aid in ass_record['IdList']:
            ass_info = Entrez.read(Entrez.esummary(id=aid, db='bioproject'))
            #proj_des = ass_info['DocumentSummarySet']['DocumentSummary'][0].get('Project_Description','')
            proj_dict = ass_info['DocumentSummarySet']['DocumentSummary'][0]
            return proj_dict

from tqdm import tqdm
org2proj_dict = {}
for org in tqdm(t):
    if org not in org2proj_dict:
        try:
            org2proj_dict[org] = org2bioproj_dict(org)
        except:
            org2proj_dict[org] = {}
            
            
# from tid to biosample
def fetch_metadata(org):
    record_dict = {}
    raw_dict = {}
    genome_record = Entrez.read(Entrez.esearch(term=org, db='genome'))
    if genome_record['IdList']:
        gid = genome_record['IdList'][0]
        
        genome_info = Entrez.read(Entrez.esummary(id=gid, db='genome'))
        genome_def = genome_info[0]['DefLine']
        
    ass_record = Entrez.read(Entrez.esearch(term=org, db='assembly'))
    if record.get('IdList', []):
        record = Entrez.read(Entrez.esummary(id=','.join([str(_) for _ in record.get('IdList', [])]), db='biosample'))
        re_text = record['DocumentSummarySet']['DocumentSummary']
        for each in re_text:
            sample_data = each.get('SampleData', '')
            if sample_data:
                soup = BeautifulSoup(sample_data)
                all_attr = soup.findAll('attribute')
                for attr in all_attr:
                    if 'env' in attr.attrs.get('attribute_name', '') or 'iso' in attr.attrs.get('attribute_name', ''):
                        n_d = attr.attrs
                        n_d['text'] = attr.text

                        record_dict[tid].append(n_d)
                n_d = {}
                n_d['attribute_name'] = 'basic_info'
                n_d['BioSample'] = each.get('Accession', '')
                n_d['Organism'] = each.get('Organism', '')
                n_d['Taxonomy'] = each.get('Taxonomy', '')
                n_d['Identifiers'] = each.get('Identifiers', '')
                record_dict[tid].append(n_d)
                raw_dict[tid].append(all_attr)
    return record_dict, raw_dict