#!/usr/bin/env python
# coding: utf-8

# # GUI hbond

# ### Initialization

# In[73]:


import pandas as pd
import subprocess
import os
from pymol import cmd, stored
import platform
from typing import List


# ### Functions

# In[20]:


"""
def form(pdbstr, x):
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

#cmd.extend('hbsearch', hbsearch)"""


# In[74]:


def changeDirectory(programDirectory: str = "."):
    
    cmd.cd(os.path.normpath(programDirectory))
    #os.chdir(os.path.normpath(programDirectory)) #can be deleted?
    # Usefull if we want to give an error! To tell the person in which directory they are located!
    #cwd = os.getcwd()


# In[75]:


def useObject(input_molecule: str):
    
    cmd.save(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"), input_molecule)


# In[76]:


def removeObject(input_molecule: str):
    
    os.remove(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"))


# In[77]:


def fetchPDB(pdbID: str, object_name: str = ""):
    if object_name == "":
        object_name = pdbID
    #setting fetch_path to desired folder
    cmd.set("fetch_path", os.path.normpath("./pdb_files/"))
    #fetching pdb file if not in folder
    cmd.fetch(pdbID, name = object_name, type = "pdb")


# In[78]:


def startHBsearch(molecule: str, hb_file: str, solvent_key:str, pse_file:str, connections: str):

    # Setting environment variable
    os.environ['PSE_FILE'] = pse_file
    # Determine operation system
    system = platform.system()
    # Executing hb_search
    hbs_output = subprocess.run(os.path.normpath(f"./{system}/hb-search -hb {hb_file} -solv {solvent_key} -con {connections} ./pdb_files/{molecule}.pdb"), capture_output=True, shell=True, check = True, text = True).stdout
    return hbs_output


# In[79]:


def readInHBS(hbs_output: str):
    hbs_rows = [i for i in hbs_output.split('\n')]
    hbs_split = [i.split() for i in hbs_rows]
    
    HEAD_LST = ['IDENT','ACC','sep1','DONO',':','x','y','z','sep2','a','b']
    
    df = pd.DataFrame(hbs_split, columns = HEAD_LST)
    df = df[df["IDENT"] == "HBOND"]
    return df


# In[80]:


def prepareLists(dataframe: pd.DataFrame):
    
    acceptor_pre = list(dataframe["ACC"])
    donor_pre = list(dataframe ["DONO"])
    
    acceptor = []
    donor = []
    
    for i in range(len(acceptor_pre)):#modify with zip
        acceptor.append(tuple(acceptor_pre[i].split(":")))
        
    for j in range(len(donor_pre)):
        donor.append(tuple(donor_pre[j].split(":")))


    
    return acceptor, donor


# In[81]:


def displayDistances(acceptor: List, donor: List, object_name: str):
    
    bondList = []
    
    for i in range(len(acceptor)): #modify with zip
        
        cmd.distance(f"{object_name}_hydrogenBond_{i}", 
                     f"{object_name}//{acceptor[i][0]}/{acceptor[i][1]}/{acceptor[i][3]}", 
                     f"{object_name}//{donor[i][0]}/{donor[i][1]}/{donor[i][3]}", )

        bondList.append(f"{object_name}_hydrogenBond_{i}")
    cmd.group(f"{object_name}_hydrogenBonds", " ".join(bondList))
    cmd.hide("labels", f"{object_name}_HydrogenBonds")


# In[90]:


def showSticks(acceptor: List,donor: List, object_name: str):
    
    stickList = []
    
    for i in range(len(acceptor)):
        stickList.append(f"/{object_name}//{acceptor[i][0]}/{acceptor[i][1]}")
        stickList.append(f"/{object_name}//{donor[i][0]}/{donor[i][1]}") 
    print(stickList)
    cmd.select(f"Connections_Sticks_{object_name}", " ".join(stickList))
    cmd.show("sticks", f"Connections_Sticks_{object_name}")
    cmd.deselect()


# In[86]:


def main(molecule:str, molecule_name = "", directory:str = ".", 
         use_object: str = "0", remove_object = "1", hb_file: str = "hb-define.txt", 
         solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", connections: str = "0"):
    
    changeDirectory(directory)
    
    print(molecule_name)
    if use_object == "0":
        fetchPDB(molecule, molecule_name)
    elif use_object == "1":
        useObject(molecule)
        
    hbs_output = startHBsearch(molecule, hb_file, solvent_key, pse_file, connections)
    hbs_dataframe = readInHBS(hbs_output)
    acceptor, donor = prepareLists(hbs_dataframe)
    
    if molecule_name == "":
        displayDistances(acceptor, donor, molecule)
        showSticks(acceptor,donor, molecule)
    else:
        displayDistances(acceptor, donor, molecule_name)
        showSticks(acceptor,donor, molecule_name)
    
    if remove_object == "1":
        removeObject(molecule)


# In[87]:


cmd.extend("hbsearch", main)


# In[72]:





# ### Main body

# In[34]:





# In[35]:




