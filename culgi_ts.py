#!/usr/bin/python3


def first_instance(rawstr, filename):
    """takes a string and return the first instance line number

    :rawstr: raw string
    :filename: filename
    :returns: line number

    """

    import re

    regex = re.compile(rawstr)

    with open(filename, 'r') as f:
        n = 0
        for line in f:
            n += 1
            if regex.search(line):
                return(n)
        raise ValueError("String not found")


def read_culgi_ts(ctf_file):
    """function able to read the timeseries file obtained from culgi

    :ctf_file: file containing the timeseries (ctf format)
    :returns: a pandas dataframe containing all ts.

    """

    import numpy as np
    import pandas as pd

    fi = first_instance(r'^Data', ctf_file)

# dropna(1) removes the created (empty) last column
    df = pd.read_csv(
            ctf_file, index_col="#Index", sep='\t',
            skiprows=list(np.append(np.arange(fi), fi+1))).dropna(1)

    return(df)


def read_culgi_tstime(ctf_file, nmol=None, dropna=True):
    """read the timeseries file obtained from culgi and converts the index
    into time format
    :ctf_file: file containing the timeseries (ctf format)
    :returns: a pandas dataframe containing all ts and index in datetime
    format
    """

    df = read_culgi_ts(ctf_file)
    df = df.set_index(df.index.to_datetime())

# if nmol is passed the Solubility is calculated
    if nmol:
        import numpy as np
        nAv = 6.022e23
        convCtte = -1.0 * float(nmol) * 1.0e3 / (nAv * 1e-24)

        df['Solubility (cal/cc)'] = np.sqrt(convCtte * (
            df['Electrostatics Energy (kcal/mol)'] +
            df['VdW Energy (kcal/mol)'] +
            df['Hydrogen bonding Energy (kcal/mol)']) /
            (df['X (A)'] * df['Y (A)'] * df['Z (A)']))

    if dropna:
        return(df.dropna())

    return(df)


def culgi_ts_panel(chainf, nmol=None, dropna=True):
    """creates a pandas panel from a list of ctf file paths
    :chainf: files to be included (obtained as chainf = !ls */*.ctf)
    :nmol: number or molecs to calculate solubility [optional]
    :returns: a pandas panel containing each file under the name of the
    previous directory:
            path/to/nice/file.ctf  ---> panel['nice'] = DF_from_file.ctf
    """
    import os.path as op
    import pandas as pd
    import pprint

    data = {}

    print("Including %d files to panel" % len(chainf), end="")

    count = 0
    err_list = []
    for fp in chainf:
        # base = op.basename(fp)
        # print("reading ", fp)
        path = op.dirname(fp)
        last = path.split("/")[-1]

        temp_df = read_culgi_tstime(fp, nmol, dropna)
        if not temp_df.empty:
            data[last] = temp_df
        else:
            count += 1
            err_list.append(last)

    if count > 0:
        print(" (%d rejections)" % count)
        print("\t", end="")
        pprint.pprint(err_list)
    else:
        print()

    return(pd.Panel(data))


if __name__ == '__main__':
    pass
