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


def read_culgi_hist_ang(ctf_file):
    """function able to read the angle histogram

    :ctf_file: file containing the histogram
    :returns: a pandas dataframe containing histogram.

    """

    import numpy as np
    import pandas as pd

    fi = first_instance(r'^Data', ctf_file)

# dropna(1) removes the created (empty) last column
    df = pd.read_csv(
            ctf_file, index_col="#Angle (degrees)", sep='\t',
            skiprows=list(np.append(np.arange(fi), fi+1))).dropna(1)

    return(df)


def read_culgi_hist_bond(ctf_file):
    """function able to read the bond histogram

    :ctf_file: file containing the histogram
    :returns: a pandas dataframe containing histogram.

    """

    import numpy as np
    import pandas as pd

    fi = first_instance(r'^Data', ctf_file)

# dropna(1) removes the created (empty) last column
    df = pd.read_csv(
            ctf_file, index_col="#Length", sep='\t',
            skiprows=list(np.append(np.arange(fi), fi+1))).dropna(1)

    return(df)

if __name__ == '__main__':
    pass
