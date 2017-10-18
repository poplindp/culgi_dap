#!/usr/bin/python3



def get_line_number(rawstr, filename):
    """takes a string (maybe a regex) and return the first instance line number
    :rawstr: raw string
    :filename: filename
    :returns: line number
    """

    import re
    regex = re.compile(rawstr)

    with open(filename, 'r') as myfile:
        for num, line in enumerate(myfile, 1):
            if regex.search(line):
                return(num)
        raise ValueError("String not found")


#
# Functions for ctf df creation
#

def read_ctf(ctf_file, time_step=0.5):
    """function able to read the timeseries file obtained from culgi
    :ctf_file: file containing the timeseries (ctf format)
    :time_step: simulation time step
    :returns: a pandas dataframe containing the ts.
    """

    import pandas as pd

# dropna(1) removes the created (empty) last column
    ctf = pd.read_csv(
        ctf_file, index_col="#Index", sep='\t', parse_dates=True,
        skiprows=get_line_number('^Data:', ctf_file))

# remove Unnamed empty column (read_cvs bug), rename and rescale index column
    ctf = ctf.loc[:, ~ctf.columns.str.contains('^Unnamed')]
    ctf.index.rename("Time (ns)", inplace=True)
    ctf.index = ctf.index * time_step

    ctf.filename = ctf_file

    return(ctf)


def ctf_intra(ctf, intra_suffix="_IntraEne"):
    """
    takes a ctf file and a ctf intraenergy file to build and "old format a"
    """
    import sys

    if any("Intermolecular" in word for word in ctf.columns):
        print("Intermolecular terms already present.  Not doing anything.")
        return(ctf)

    ctf_intrafilename = ctf.filename.replace("_Inst", intra_suffix)
    try:
        ctf_intra = read_ctf(ctf_intrafilename)
    except IOError:
        print("Cannot open Intramolecular energy file: '{}'".format(ctf_intrafilename)
            , file=sys.stderr)
        sys.exit(1)

# adding Intramolecular and units strings
    ctf_intra.columns = ['Intramolecular ' + str(col) + ' (kcal/mol)'
            for col in ctf_intra.columns]

    ctf = ctf.join(ctf_intra)

# Substracting intramolecular contributions

    el_lbl = "Electrostatics Energy (kcal/mol)"
    vdw_lbl = "VdW Energy (kcal/mol)"
    hb_lbl = "Hydrogen bonding Energy (kcal/mol)"
    hb_lbl2 = "Hydrogen Bonding Energy (kcal/mol)"

    ctf[el_lbl] = ctf[el_lbl] - ctf["Intramolecular "+ el_lbl]
    ctf[vdw_lbl] = ctf[vdw_lbl] - ctf["Intramolecular "+ vdw_lbl]
    ctf[hb_lbl] = ctf[hb_lbl] - ctf["Intramolecular " + hb_lbl2]

    if (ctf["Intramolecular " + hb_lbl2] != 0).any():
        print("Intramolecular " + hb_lbl2 + "different to zero!")

    ctf = ctf.rename(columns={el_lbl: 'Intermolecular ' + el_lbl})
    ctf = ctf.rename(columns={vdw_lbl: 'Intermolecular ' + vdw_lbl})
    ctf = ctf.rename(columns={hb_lbl: 'Intermolecular ' + hb_lbl})

    return(ctf)


def solubility(ctf, nmol=1, dropna=False):
    """read the timeseries file obtained from culgi and converts the index
    into time format
    :ctf_file: file containing the timeseries (ctf format)
    :returns: a pandas dataframe containing all ts and index in datetime
    format
    """

    import pandas as pd

    if not any("Intermolecular" in word for word in ctf.columns):
        print("Intermolecular terms not present. Calculating them.")
        ctf = ctf_intra(ctf)

# FIXME is this convertion really necessary?
    # df = df.set_index(pd.to_datetime(df.index))

    import numpy as np
    n_av = 6.022e23
    conv_ctte = -1.0 * float(nmol) * 1.0e3 / (n_av * 1e-24)

    prefix_lbl = "Intermolecular "
    el_lbl = prefix_lbl + 'Electrostatics Energy (kcal/mol)'
    vdw_label = prefix_lbl + 'VdW Energy (kcal/mol)'
    hb_lbl = prefix_lbl + 'Hydrogen bonding Energy (kcal/mol)'

    ctf['Solubility (cal/cc)'] = np.sqrt(conv_ctte *
            (ctf[el_lbl] + ctf[vdw_label] + ctf[hb_lbl]) /
            (ctf['X (A)'] * ctf['Y (A)'] * ctf['Z (A)']))

    if dropna:
        return(ctf.dropna())

    return(ctf)


# This function should only be called (mostly) inside build_dfmi
def read_culgi_descriptors(out_file):
    """function able to read the descriptors file obtained from culgi
    :out_file: file containing the timeseries (ctf format)
    :returns: a pandas series containing all descriptors.
    print("descriptors")
    """

    import pandas as pd

    df = pd.read_csv(out_file, header=None, index_col=0, sep='[=|\t]',
            engine='python', skiprows=get_line_number('^------', out_file),
            usecols=[0,1], skipfooter=2)
# converting to numbers
    df[1] = pd.to_numeric(df[1], errors='coerce', downcast='signed')

    return(df)


#
# Functions dfmi creation
#

# FIXME nmol should be read from the input file
#       in fact input files should be added

def build_dfmi(system, sample, files="mm06*.ctf", nmol=1,
               axis=0, verb=False, descrbool=True):
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

    It requires a given directory format:
        ./system/"runs/"sample/files
    """

    import pandas as pd
    import os.path
    import glob
    import re

    dict1mf = {}
    system_list = [file for file in glob.glob(system, recursive=True)]

    for sys in system_list:
        dict2mf = {}
        sys_name = os.path.basename(sys)
        if verb:
            print(sys_name)
        sample_list = [file for file in
                       glob.glob(sys + "/runs/" + sample, recursive=True)]

        for sam in sample_list:
            dict3mf = {}
            sam_name = os.path.basename(sam)
            if verb:
                print("\t" + sam_name)
            files_list = [file for file in
                          glob.glob(sam + "/" + files+"_Inst.ctf", recursive=True)]

            for fil in files_list:
                fil_name = os.path.basename(fil)
                if verb:
                    print("\t\t" + fil_name)

                ctf = read_ctf(fil)
                ctf = solubility(ctf, nmol)


# FIXME Here also input file should be added.
                fil_descrp = re.sub(r'_Inst.*.ctf', r'.cof_descriptors.out',
                                    fil)
                fil_exists = os.path.exists(fil_descrp)
#                if verb:
#                    print(fil_descrp)
#                if descrip :
                if fil_exists and descrbool:
                    if verb:
                        print("\t\t--> " + fil_descrp)
                    cdes = read_culgi_descriptors(fil_descrp)
                    dict3mf[fil_name.split('_')[0]] = pd.concat([ctf.mean(),
                                                                cdes[1]],
                                                                axis=0)
                else:
                    dict3mf[fil_name.split('_')[0]] = ctf.mean()
            dict2mf[sam_name.split('-')[-1]] = pd.concat(dict3mf, axis=1)
        dict1mf[sys_name] = pd.concat(dict2mf, axis=1)
    dfmi_full = pd.concat(dict1mf, axis=1).transpose()

    return(dfmi_full)
