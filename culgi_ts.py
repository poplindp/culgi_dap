#!/usr/bin/python3
#
# Functions for ctf df creation
#

def read_cts(cts_csv, solub=True, intra_corr=True, intra_suffix="_IntraEne",
        time_step=0.5, nmol=1, debug=True):
    """function able to read the timeseries file obtained from culgi and
    correcting, if present, by intramolecular terms

    :cts_csv: CSV culgi time series
    :solub: add solubility column if possible
    :intra_corr: if Intramolecular Energy correction is to be done
    :intra_suffix: suffix of CSV containing intramolecular energies
    :time_step: simulation time step
    :nmol: number of polymer molecules (for solubility calculation)
    :debug: print information
    :returns: a pandas dataframe containing the time series.
    """

    import sys
    import pandas as pd

    #df = pd.read_csv(cts_csv, sep='\t', index_col="Index", parse_dates=True)
    df = pd.read_csv(cts_csv, index_col="Index", parse_dates=True)

    if solub or intra_corr:
        cts_intra_csv = cts_csv.replace("_Inst", intra_suffix)  # change name
        df = intra_energy(df, cts_intra_csv)

    if solub:
        solubility(df, nmol=nmol)

# rename and rescale index column
    df.index.rename("Time(ns)", inplace=True)
    df.index = df.index * time_step

# FIXME: this was useful if internal analysis was done with
# assign filename value (to keep track where this was coming from)
    df.filename = cts_csv


    return(df)



def intra_energy(df, cts_intra_csv, debug=True):
    """ Given a culgi time series df with total energies and a file with intra
    enegies, modify df to add Inter and Intramolecular Energy cols

    :df: culgi time series dataframe containing total energies
    :cts_intra_csv: CSV culgi time series containing intra energies
    :debug: print information
    :returns: a pandas dataframe containing the time series with intra ener
    """

    import sys
    import pandas as pd

    if not any("Intermolecular" in word for word in df.columns):

        try:
            df_intra = pd.read_csv(cts_intra_csv, index_col="Index",
                    parse_dates=True)
        except IOError:
            if debug:
                print("Cannot open Intramolecular energy file: '{}'"
                        .format(cts_intra_csv) , file=sys.stderr)
            sys.exit(1)

    # adding Intramolecular and units strings
        df_intra.columns = ['Intramolecular' + str(col) + '(kcal/mol)'
                for col in df_intra.columns]

    # NOTE: here join forces to return a value since df is a new df, is a copy
        df = df.join(df_intra)

    # Substracting intramolecular contributions

        el_lbl = "ElectrostaticsEnergy(kcal/mol)"
        vdw_lbl = "VdWEnergy(kcal/mol)"
        hb_lbl = "HydrogenbondingEnergy(kcal/mol)"
        hb_lbl2 = "HydrogenBondingEnergy(kcal/mol)"

        df[el_lbl] = df[el_lbl] - df["Intramolecular"+ el_lbl]
        df[vdw_lbl] = df[vdw_lbl] - df["Intramolecular"+ vdw_lbl]
        df[hb_lbl] = df[hb_lbl] - df["Intramolecular" + hb_lbl2]

        if (df["Intramolecular" + hb_lbl2] != 0).any():
            if debug:
                print("Intramolecular " + hb_lbl2 + "different to zero!")

        df = df.rename(columns={el_lbl: 'Intermolecular' + el_lbl})
        df = df.rename(columns={vdw_lbl: 'Intermolecular' + vdw_lbl})
        df = df.rename(columns={hb_lbl: 'Intermolecular' + hb_lbl})

    return(df)



def solubility(df, nmol=1, debug=True):
    """ Given a df containing intramolecular and intermolecular energies
    calculate and add a new solubility column
    :df: culgi dataframe containing the energies
    """

    import sys
    import pandas as pd
    import numpy as np

    if any("Solubility" in word for word in df.columns):
        if debug:
            print("  Solubility column already present")
#        return(df)

    n_av = 6.022e23
    conv_ctte = -1.0 * float(nmol) * 1.0e3 / (n_av * 1e-24)

    prefix_lbl = "Intermolecular"
    el_lbl = prefix_lbl + 'ElectrostaticsEnergy(kcal/mol)'
    vdw_label = prefix_lbl + 'VdWEnergy(kcal/mol)'
    hb_lbl = prefix_lbl + 'HydrogenbondingEnergy(kcal/mol)'

    try:
        df['Solubility(cal/cc)'] = np.sqrt(conv_ctte *
                (df[el_lbl] + df[vdw_label] + df[hb_lbl]) /
                (df['X(A)'] * df['Y(A)'] * df['Z(A)']))
    except KeyError:
        print("DataFrame not containing Inter or Intermolecular Energy terms")
        sys.exit(1)

#   in case to need return a df, instead of calculating df['Solub...'] we could
#   have used a foo['Solub...'] and do df=df.join(foo['Solub...'])
#    return(df)


def read_culgi_descriptors(cvs_file, singlerow=True):
    """Takes a csv file containing descriptors for all molecules
    it returns a dataframe with the mean of all descriptors for everymolecule
    type "mol_name" (contains filename used by culgi) column is removed and
    error strings are converted into nan
    :cvs_file: file containing the timeseries (ctf format)
    :single_row: returns all descriptors in a single row df
    :returns: a pandas df containing all descriptors.
    print("descriptors")
    """

    import pandas as pd
    import numpy as np

    # read and remove "mol_name" column
    df = pd.read_csv(cvs_file).drop("mol_name", axis=1)

    # FIXME: this line is prone error, there is no a good way to determine
    # which name the polymer is supposed to have.
    # replace word with "-" usually polymer should contain that
    df['molecule'] = df['molecule'].replace(r'\b\w+-\w+\b', 'polymer', regex=True)

    # set index to molecule
    df = df.set_index('molecule')

    # replace string cells (usually errors) into nan
    df.replace(r'\s+', np.nan, regex=True, inplace=True)

    # convert everything to numeric
    df = df.apply(pd.to_numeric)

    # caculate mean (separating by index)
    df = df.mean(level='molecule')

    if singlerow:
        # fuse rows adding index name to columns
        df = df.stack(dropna=False).to_frame().T
        df.columns = ['{}_{}'.format(*c) for c in df.columns]

    return(df)



#
# Functions dfmi creation
#

# FIXME nmol should be read from the input file
#       in fact input files should be added

def build_dfmi(files="*_Inst.csv", nmol=1, multindex=False, axis=0,
        verb=False, descrbool=True, condit=True):
    """function building a multiindex dataframe containing the averages
    of culgi time series
    if descriptor files are present, add it to table
    :files: = name of csv energy files (*Inst.csv avoids problem with
    IntraEne.csv
    :nmol: # number of molecules per box  (for solubility calculation)
    :axis: # if multiindex are in rows (0, uglier but useful to use
    with .loc or .iloc) or columns (1, just nicer)
    :condit: includes initial contidions
    :verb: prints which file is taking
    :returns: a pandas dataframe multindex containing all <ts> and descript.

    It requires a given directory format:
        ./system/"runs/"sample/files
    """

    import pandas as pd
    import os.path as osp
    import glob
    import re


    search_path = "**/" + files

    result = pd.DataFrame()

    for fileraw in glob.iglob(search_path, recursive=True):

        lst = []

        # Skip internal energy files
        if "IntraEne" in fileraw:
            continue

        file = osp.normcase(osp.normpath(fileraw))
        if verb:
            print("Processing:", file)

        # Processing path to extract polymer, system, condition, replica
        filename = osp.basename(file)
        dirname = osp.dirname(file)
        (polsys, rep) = osp.split(dirname)
        (polymer_full, syscond) = osp.split(polsys)
        polymer = osp.basename(polymer_full)
        (system, condition)=syscond.split("-")

        try:
            replica = int(re.search('.*-(\d+)$', rep)[1])
        except:
            print("** Problem parsing replica in:", fileraw)
            break


# add the columns that identify the system
        ini = pd.DataFrame({'Polymer': polymer, 'Replica': replica,
            'Filename': filename}, index=[0])

        lst.append(ini)

# add conditions and systems

        if condit:
            cond = pd.read_csv(dirname+"/data/Conditions.csv")
            conditions = cond[cond['Condition'].str.contains(condition)]
            conditions.index = [0] # Need to reset this index to 0
            lst.append(conditions)

            syst = pd.read_csv(dirname+"/data/Systems.csv")
            systems = syst[syst['System'].str.contains(system)]
            systems.index = [0] # Need to reset this index to 0
            nmol = systems['nChains']
            lst.append(systems)

        else:
            ini['Condition'] = condition
            ini['System'] = system


# read energy time series
        cts_ener = read_cts(fileraw, solub=True, nmol=nmol, intra_corr=True,
               intra_suffix="_IntraEne", time_step=0.5, debug=True)
        ener_mean = cts_ener.mean().to_frame().T

        lst.append(ener_mean)

# FIXME: calculate autocorr

# add descriptors
        fil_descrp = re.sub(r'_Inst', r'_desc', fileraw)
        fil_exists = osp.exists(fil_descrp)
        if fil_exists and descrbool:
            if verb:
                print("  Descriptors --> " + fil_descrp)
            cdes = read_culgi_descriptors(fil_descrp)
            lst.append(cdes)

        if verb:
            print()

        result = result.append(pd.concat(lst, axis=1))

    if multindex:
        result = result.set_index(['Polymer', 'System', 'Condition', 'Replica', 'Filename'])

    return(result)




# --- From here older versions

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

def read_ctf_old(ctf_file, time_step=0.5):
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


def ctf_intra_old(ctf, intra_suffix="_IntraEne"):
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


def solubility_old(ctf, nmol=1, dropna=False):
    """read the timeseries file obtained from culgi and converts the index
    into time format
    :ctf_file: file containing the timeseries (ctf format)
    :returns: a pandas dataframe containing all ts and index in datetime
    format
    """

    import pandas as pd

    if not any("Intermolecular" in word for word in ctf.columns):
        print("  Intermolecular terms not present. Calculating them.")
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
def read_culgi_descriptors_old(out_file):
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


def build_dfmi_OLD(system, sample, files="mm06*.ctf", nmol=1,
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
    #FIXME: would be great to only list directories... 
    system_list = [file for file in glob.glob(system, recursive=True)]

    for sys in system_list:
        dict2mf = {}
        sys_name = os.path.basename(sys)
        if verb:
            print(sys_name)
    #FIXME: same here it would be great to only list directories... 
        sample_list = [file for file in
                       glob.glob(sys + "/" + sample, recursive=True)]
                      # glob.glob(sys + "/runs/" + sample, recursive=True)]

        for sam in sample_list:
            dict3mf = {}
            sam_name = os.path.basename(sam)
            if verb:
                print("\t" + sam_name)
            files_list = [file for file in
                          glob.glob(sam + "/" + files+"_Inst.ctf", recursive=True)]

            print(files_list)
            for fil in files_list:
                fil_name = os.path.basename(fil)
                if verb:
                    print("\t\t" + fil_name)

                ctf = read_ctf(fil)
                ctf = solubility(ctf, nmol)


# FIXME Here also input file should be added.
                #fil_descrp = re.sub(r'_Inst.*.ctf', r'.cof_descriptors.out',
                fil_descrp = re.sub(r'_Inst.ctf', r'_desc.out',
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
