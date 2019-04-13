from Bio import Entrez
from bs4 import BeautifulSoup

Entrez.email = 'l0404th@gmail.com'


# from tid to biosample
def fetch_metadata(tid):
    record_dict = {}
    raw_dict = {}
    record = Entrez.read(Entrez.esearch(term='txid' + tid, db='biosample'))
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
