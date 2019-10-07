from Bio import Entrez
from ete3 import NCBITaxa
ncbi = NCBITaxa()
Entrez.email = 'l0404th@gmail.com'

def pid2taxid(pid):
    taxid_r = Entrez.read(Entrez.elink(dbfrom='protein',id=pid,db='taxonomy'))
    try:
        tid = taxid_r[0]['LinkSetDb'][0]['Link'][0]['Id']
    except:
        return
    
    return tid
    
def print_tax(pid):
    tid = pid2taxid(pid)
    if tid is None:
        print('')
        return
    lineage = ncbi.get_lineage(int(tid))
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    rank2name = {v:names[t] for t,v in ranks.items()}

    if rank2name.get('phylum','') == 'Proteobacteria':
        print(rank2name['class'])
    else:
        phy = rank2name.get('phylum','')
        if not phy:
            print(tid)
        else:
            print(phy)