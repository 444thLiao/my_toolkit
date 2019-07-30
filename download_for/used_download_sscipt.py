#################################################################################
####  for download genome info from refseq of NCBI
####
####
#################################################################################
import pandas as pd
import os
from ftplib import FTP
from os.path import join, basename
from tqdm import tqdm
import click
import wget


def get_dirname_multiple(path, n=1):
    _count = 0
    path = os.path.realpath(os.path.abspath(path))
    while _count < n:
        path = os.path.dirname(path)
        _count += 1
    return path


def get_constant():
    refseq_link = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    domain_refseq_link = {}
    for domain in ["archaea",
                   "bacteria",
                   "fungi",
                   "plant",
                   "viral",
                   "invertebrate",
                   "vertebrate_mammalian",
                   "vertebrate_other"]:
        domain_refseq_link[domain] = "ftp://ftp.ncbi.nih.gov/genomes/refseq/%s/assembly_summary.txt" % domain
    collect_dict = dict(refseq=refseq_link,
                        domain_refseq=domain_refseq_link)
    return collect_dict


def rewrite_md5(infile):
    dirname = get_dirname_multiple(infile, 2)
    list_dir = [os.path.basename(_).lower()
                for _ in os.listdir(dirname)]
    result = ''
    with open(infile) as f1:
        for row in f1:
            suffix = row.split('.')[-2]
            if suffix in list_dir:
                result += suffix
    with open(infile, 'w') as f1:
        f1.write(result)


@click.command()
@click.option("--odir", '-o', required=True)
@click.option("--target", '-t', required=False, default='all')
@click.option("--infile", '-i', required=False, default=None)
@click.option("--only_infile", is_flag=True, )
def cli(odir, target, infile):
    odir = os.path.abspath(odir)
    collect_dict = get_constant()
    if infile is None:
        if target == 'all':
            infile = collect_dict['refseq']
        elif target.lower() in collect_dict["domain_refseq"]:
            infile = collect_dict["domain_refseq"][target.lower()]
        else:
            raise SyntaxError("target is not inside required list, infile also not given.")
        wget.download(infile, out=odir)
        infile = join(odir, basename(infile))

    if (not os.path.exists(infile)):
        raise SyntaxError("infile doesn't exist.")
    else:
        main(odir, infile, target)


def main(odir, infile, target, max_try=5):
    assembly_summary_df = pd.read_csv(infile, sep='\t', header=None, comment='#', low_memory=False)
    # download from assembly_summary.txt

    odir = os.path.join(odir, target)
    os.makedirs(odir, exist_ok=True)
    ftp = FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()

    # make dirs
    os.makedirs(join(odir, "MD5"), exist_ok=True)
    type_ = ['fna', 'gbff']
    for t in type_:
        os.makedirs(join(odir, t), exist_ok=True)

    for idx, row in tqdm(assembly_summary_df.iterrows(),
                         total=assembly_summary_df.shape[0]):
        path = row[19]
        if pd.isna(path):
            # in case there is missing path
            continue
        cwd_path = '/' + '/'.join(path.split('/')[3:])
        name = path.split('/')[-1]
        type_ = ['fna', 'gbff']

        if ftp.pwd() != cwd_path:
            ftp.cwd(cwd_path)

        for type_file in type_:
            download_target = name + '_genomic.%s.gz' % type_file
            ofile = join(odir, type_file, download_target)
            if os.path.exists(ofile):
                continue
            while not os.path.exists(ofile):
                count = 1
                try:
                    count += 1
                    ftp.retrbinary('RETR %s' % download_target,
                                   open(ofile,
                                        'wb').write)
                except:
                    if count >= max_try:
                        break

        download_target = "md5checksums.txt"
        ofile = join(odir, "MD5", "%s.md5" % name)
        if os.path.exists(ofile):
            continue
        if ftp.pwd() != cwd_path:
            ftp.cwd(cwd_path)
        ftp.retrbinary('RETR %s' % download_target,
                       open(ofile,
                            'wb').write)
        rewrite_md5(ofile)
        # back to the root path
        ftp.cwd('/')
