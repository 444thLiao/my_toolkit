from glob import glob
from tqdm import tqdm
from subprocess import check_call
import os

prokka_p = "/usr/local/bin/prokka"


def main(idir, odir):
    for seq in tqdm(glob(os.path.join(idir, '*.fa'))):
        # sample_name = os.path.basename(seq).split('.')[0]
        sample_name = '_'.join(os.path.basename(seq).split('_')[:2])
        os.makedirs(f"{odir}/{sample_name}", exist_ok=True)
        cmd = "{prokka_p} {infile} --outdir {odir}/{sample_name} --force --cpus 0 "

        check_call(cmd.format(prokka_p=prokka_p,
                              infile=seq,
                              odir=odir,
                              sample_name=sample_name,),
                   executable="/home-user/thliao/anaconda3/bin/zsh",
                   shell=True,
                   stdout=open(os.path.join(odir,'prokka_o.log'),'w'),
                   stderr=open(os.path.join(odir,'prokka_o.log'),'w'))

if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    main(args.indir, args.outdir)
