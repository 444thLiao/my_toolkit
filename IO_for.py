import csv
import pandas as pd

def data_parser(path, ft='csv', verbose=1, **kwargs):
    if type(path) != str and ft != 'metadatas':
        df = path.copy()
        if type(df) != pd.DataFrame:
            df = pd.DataFrame(df)
    else:
        if ft == 'csv':
            sniffer = csv.Sniffer()
            sniffer.preferred = [',', '\t', '|']
            dialect = sniffer.sniff(open(path, 'r').readline().strip('\n'))
            df = pd.read_csv(path, sep=dialect.delimiter, index_col=False, header=0, low_memory =False,**kwargs)
            # low_memory is important for sample name look like float and str. it may mistaken the sample name into some float.
        elif ft == 'metadatas': # indicate there may be multiple input
            datas = [data_parser(_,ft='csv',verbose=0) for _ in path]
            assert len(set([tuple(_.index) for _ in datas])) == 1  # all datas must have same samples
            merged_data = pd.concat(datas, axis=1)
            col_dict = {p: list(data.columns) for p, data in zip(path,datas)}
            return merged_data, col_dict
        else:
            df = pd.read_excel(path, index_col=0, header=0, **kwargs)

    return df