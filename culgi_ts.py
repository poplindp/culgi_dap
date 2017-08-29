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

    import pandas as pd
    df = read_culgi_ts(ctf_file)
# df = df.set_index(df.index.to_datetime())
    df = df.set_index(pd.to_datetime(df.index))

# if nmol is passed the Solubility is calculated
    if nmol:
        import numpy as np
        nAv = 6.022e23
        convCtte = -1.0 * float(nmol) * 1.0e3 / (nAv * 1e-24)

# adapting for old versions
        prefix_lbl = ""
        if sum(df.columns.str.contains(r'Intermolecular')):
            prefix_lbl = "Intermolecular "
        el_lbl = prefix_lbl + 'Electrostatics Energy (kcal/mol)'
        vdw_label = prefix_lbl + 'VdW Energy (kcal/mol)'
        hb_lbl = prefix_lbl + 'Hydrogen bonding Energy (kcal/mol)'

        df['Solubility (cal/cc)'] = np.sqrt(convCtte * (
            df[el_lbl] + df[vdw_label] + df[hb_lbl]) /
            (df['X (A)'] * df['Y (A)'] * df['Z (A)']))

    if dropna:
        return(df.dropna())

    return(df)


def read_culgi_descriptors(out_file):
    """function able to read the descriptors file obtained from culgi
    :out_file: file containing the timeseries (ctf format)
    :returns: a pandas series containing all descriptors.
    print("descriptors")
    """

    import numpy as np
    import pandas as pd
    fi = first_instance(r'^------', out_file)+1

# dropna(2,1) removes the 2nd (empty) column
    df = pd.read_csv(
            out_file, header=None, index_col=0, sep='[=|\t]', engine='python',
            skiprows=list(np.arange(fi))).drop(2, 1)
    df = df.drop(df.index[len(df)-1])  # removing last row
# converting to numbers
    df[1] = pd.to_numeric(df[1], errors='coerce', downcast='signed')
    return(df)


def build_dfmi(system, sample, files="mm06*.ctf", nmol=1,
               axis=0, verb=False):
    """function building a multiindex dataframe containing the averages
    of culgi time series
    if descriptor files are present, add it to table
    :system: = 'P*' name of system
    :sample: = '*' name of sample
    :files: = name of ctf files
    :nmol: # number of molecules per box  (for solubility calculation)
    :axis: # if multiindex are in rows (0, uglier but useful to use
    with .loc or .iloc) or columns (1, just nicer)
    :verb: prints which file is taking
    :returns: a pandas dataframe multindex containing all <ts> and descript.
    """

    import pandas as pd
    import os.path
    import glob
    import re

    dict1mF = {}
    system_list = [file for file in glob.glob(system, recursive=True)]

    for sys in system_list:
        dict2mF = {}
        sys_name = os.path.basename(sys)
        if verb:
            print(sys_name)
        sample_list = [file for file in
                       glob.glob(sys + "/" + sample, recursive=True)]

        for sam in sample_list:
            dict3mF = {}
            sam_name = os.path.basename(sam)
            if verb:
                print("\t" + sam_name)
            files_list = [file for file in
                          glob.glob(sam + "/" + files, recursive=True)]

            for fil in files_list:
                fil_name = os.path.basename(fil)
                if verb:
                    print("\t\t" + fil_name)
                ctf = read_culgi_tstime(fil, nmol)
#                dict3[fil_name.split('_')[0]] = ctf
#                dict3m[fil_name.split('_')[0]] = ctf.mean()

                fil_descrp = re.sub(r'_Inst.*.ctf', r'.cof_descriptors.out',
                                    fil)
                fil_exists = os.path.exists(fil_descrp)
                print(fil_descrp)
#                if descrip:
                if fil_exists:
                    print("\t\t--> " + fil_descrp)
                    cdes = read_culgi_descriptors(fil_descrp)
                    dict3mF[fil_name.split('_')[0]] = pd.concat([ctf.mean(),
                                                                cdes[1]],
                                                                axis=0)
                else:
                    dict3mF[fil_name.split('_')[0]] = ctf.mean()
            dict2mF[sam_name.split('-')[-1]] = pd.concat(dict3mF, axis=1)
        dict1mF[sys_name] = pd.concat(dict2mF, axis=1)
    dfmi_full = pd.concat(dict1mF, axis=1).transpose()

    return(dfmi_full)


# This function is not useful anymore and a multindex dataframe should
# be created in situ using a dictionry and concatenating it
# dict_n[index_level_n] = cts.read_culgi_tstime(fil, nmol)
# dict_n-1[index_level_n-1] = pd.concat(dict_n, axis=axis)
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
