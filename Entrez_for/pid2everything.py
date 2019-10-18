from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
import io 

ncbi = NCBITaxa()
Entrez.email = 'l0404th@gmail.com'

def pid2taxid(pid):
    taxid_r = Entrez.read(Entrez.elink(dbfrom='protein',id=pid,db='taxonomy'))
    try:
        tid = taxid_r[0]['LinkSetDb'][0]['Link'][0]['Id']
    except:
        return
    
    return tid

def pid2org(pid):
    p_result = Entrez.read(Entrez.esummary(db='protein',id=pid))
    try:
        org_name = p_result[0]['Title'].split('[')[-1].split(']')[0]
        return org_name
    except:
        return ''
    
def pid2seq(pid):
    p_fa = Entrez.efetch(db='protein', id=pid, retmode='text', rettype='fasta').read()
    try:
        record = SeqIO.read(io.StringIO(p_fa),format='fasta')
        return str(record.seq)
    except:
        return ''


def pid2accession(pid):
    p_fa = Entrez.efetch(db='protein', id=pid, retmode='text', rettype='fasta').read()
    try:
        record = SeqIO.read(io.StringIO(p_fa),format='fasta')
        return str(record.id)
    except:
        return ''
    
def pid2tax(pid):
    tid = pid2taxid(pid)
    if tid is None:
        #print('')
        return ''
    lineage = ncbi.get_lineage(int(tid))
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    rank2name = {v:names[t] for t,v in ranks.items()}
    if '48479' in map(str,lineage):
        return 'ENV'
    if rank2name.get('phylum','') == 'Proteobacteria':
        return(rank2name['class'])
    else:
        phy = rank2name.get('phylum','')
        if not phy:
            return(tid)
        else:
            return(phy)