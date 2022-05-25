import pandas as pd
import subprocess
import os
from pymol import cmd, stored
import platform


def form(pdbstr, x):
    '''
    Part of lambda function to format dataframe to Pymol compatible form:
    from "A:183:LEU:O" to "/2akr/A/A/LEU`183/O"
    :param pdbstr: pdb code of handled structure
    :param x: entry within dataframe
    '''
    chain, resID, amiaci, atom = x.split(':')
    x = f'(/{pdbstr}//{chain}/{resID}/{atom})'
    return x


def pymol_display(df):
    '''
    '''
    zilis = list(zip(df['ACC'].tolist(), df['DONO'].tolist()))
    # print(zilis)
    for i in zilis:
        cmd.distance(i[0], i[1])


def startHBsearch():

    # Setting environment variable
    os.environ['PSE_FILE'] = 'period-table-info.txt'
    # Determine operation system
    osys = platform.system()
    # Executing hb_search
    hbs = subprocess.run(f"./{osys}/hb-search -hb hb-define.txt 2akr.pdb",
                         capture_output=True, shell=True, check=True, text=True).stdout
    return hbs


def readInHBS(pdbstr, hbsfile):
    hbs_columns = [i for i in hbsfile.split('\n')]
    hbs_split = [i.split() for i in hbs_columns]

    HEAD_LST = ['IDENT', 'ACC', 'sep1', 'DONO',
                ':', 'x', 'y', 'z', 'sep2', 'a', 'b']

    df = pd.DataFrame(hbs_split, columns=HEAD_LST)
    df = df[df["IDENT"] == "HBOND"]
    df = df[['ACC', 'DONO']]
    df['ACC'] = df['ACC'].map(lambda x: form(pdbstr, x))
    df['DONO'] = df['DONO'].map(lambda x: form(pdbstr, x))
    return df


def hbsearch(pdbstr: str) -> pd.DataFrame():
    '''
    Executing hb_search with set parameters and extract HBOND-entries from output
    :return df_hbond: Dataframe with all HBOND entries from hb_search output 
    '''

    hbs = startHBsearch()
    df = readInHBS(pdbstr, hbs)

    pymol_display(df)

    return df


cmd.extend('hbsearch', hbsearch)
