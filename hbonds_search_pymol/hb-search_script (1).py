#!/usr/bin/env python
# coding: utf-8

# # GUI hbond

# ### Initialization

# In[1]:


import pandas as pd
import subprocess
import os
from pymol import cmd, stored
import platform


# ### Functions

# In[33]:



"""def form(pdbstr, x):
    '''
    Part of lambda function to format dataframe to Pymol compatible form:
    from "A:183:LEU:O" to "/2akr/A/A/LEU`183/O"
    :param pdbstr: pdb code of handled structure
    :param x: entry within dataframe
    '''
    temp = x.split(':')
    x = f'/{pdbstr}//{temp[0]}/{temp[2]}/{temp[3]}'
    return x

def pymol_display(df):
    '''
    '''
    zilis = list(zip(df['ACC'].tolist(), df['DONO'].tolist()))
    for i in zilis:
        cmd.distance( i[0] , i[1])
    
def hbsearch(pdbstr:str) -> pd.DataFrame():
    '''
    Executing hb_search with set parameters and extract HBOND-entries from output
    :return df_hbond: Dataframe with all HBOND entries from hb_search output 
    '''
    
    HEAD_LST = ['IDENT','ACC','sep1','DONO',':','x','y','z','sep2','a','b']
    
    # Setting environment variable
    os.environ['PSE_FILE'] = 'period-table-info.txt'
    
    # Executing hb_search
    hbs = subprocess.run(f"./hb-search -hb hb-define.txt {pdbstr}.pdb", 
                         stdout=subprocess.PIPE, shell=True, check=True, text=True)
    
    # Decode and format hb_search output
    hbs_hbb = [i for i in hbs.stdout.split('\n') if i[0:5] == "HBOND"]
    hbs_splt = [i.split() for i in hbs_hbb]

    # Return dataframe
    df_hbond = pd.DataFrame(hbs_splt, columns = HEAD_LST)
    
    df_hbond = df_hbond[['ACC', 'DONO']]
    df_hbond['ACC'] = df_hbond['ACC'].map(lambda x: form(pdbstr,x) )
    df_hbond['DONO'] = df_hbond['DONO'].map(lambda x: form(pdbstr,x) )
    
    pymol_display(df_hbond)
    
    return df_hbond

cmd.extend('hbsearch', hbsearch)"""


# In[43]:


def changeDirectory(programDirectory = "."):
    
    cmd.cd(programDirectory)
    os.chdir(os.path.normpath(programDirectory))
    # Usefull if we want to give an error! To tell the person in which directory they are located!
    cwd = os.getcwd()
changeDirectory()


# In[55]:


def startHBsearch():

    # Setting environment variable
    os.environ['PSE_FILE'] = 'period-table-info.txt'
    # Determine operation system
    system = platform.system()
    # Executing hb_search
    hbs = subprocess.run(os.path.normpath(f"./{system}/hb-search -hb hb-define.txt 4awn.pdb"), capture_output=True, shell=True, check = True, text = True).stdout
    return hbs


# In[62]:


def readInHBS(hbsfile):
    hbs_columns = [i for i in hbsfile.split('\n')]
    hbs_split = [i.split() for i in hbs_columns]
    
    HEAD_LST = ['IDENT','ACC','sep1','DONO',':','x','y','z','sep2','a','b']
    
    df = pd.DataFrame(hbs_split, columns = HEAD_LST)
    df = df[df["IDENT"] == "HBOND"]
    return df
readInHBS(startHBsearch())


# In[69]:


def prepareLists(dataframe):
    
    acceptor_pre = list(dataframe["ACC"])
    donor_pre = list(dataframe ["DONO"])
    
    acceptor = []
    donor = []
    
    for i in range(len(acceptor_pre)):
        acceptor.append(tuple(acceptor_pre[i].split(":")))
        
    for j in range(len(donor_pre)):
        donor.append(tuple(donor_pre[j].split(":")))


    
    return acceptor, donor
acceptor, donor = prepareLists(readInHBS(startHBsearch()))


# In[74]:


def displayDistances(acceptor, donor, obj = "4awn"):
    
    bondList = []
    
    for i in range(len(acceptor)):
        print(acceptor[i], donor[i])
        cmd.distance(f"HydrogenBond{i}", 
                     f"{obj}//{acceptor[i][0]}/{acceptor[i][1]}/{acceptor[i][3]}", 
                     f"{obj}//{donor[i][0]}/{donor[i][1]}/{donor[i][3]}", )

        bondList.append(f"HydrogenBond{i}")
    cmd.group("HydrogenBonds", " ".join(bondList))
    cmd.hide("labels", "HydrogenBonds")


# In[75]:


def main():
    changeDirectory()
    hbs_output = startHBsearch()
    hbs_dataframe = readInHBS(hbs_output)
    acceptor, donor = prepareLists(hbs_dataframe)
    displayDistances(acceptor, donor)


# In[76]:


cmd.extend("hbsearch", main)


# In[ ]:





# ### Main body

# In[34]:





# In[35]:




