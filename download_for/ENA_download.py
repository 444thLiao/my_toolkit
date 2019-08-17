from ftplib import FTP
from tqdm import tqdm
ftp = FTP("ftp.sra.ebi.ac.uk")
ftp.login()

def get_filename(idlist):
    base_dir = "/vol1"
    ftp.cwd(base_dir)

    id2filenames = dict()
    for each_id in tqdm(idlist):
        dirname = each_id[:6]
        sub_dir = each_id
        ftp.cwd('/'.join([base_dir,dirname,sub_dir]))
        file_names = ftp.nlst()
        id2filenames[each_id] = file_names