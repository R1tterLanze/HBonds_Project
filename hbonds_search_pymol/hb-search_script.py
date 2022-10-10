#!/usr/bin/env python
# coding: utf-8

# # GUI hbond

# ### Initialization

# In[16]:


import pandas as pd
import subprocess
import os
from pymol import cmd, stored
import platform
from typing import List, Dict
from datetime import datetime


# ### Functions

# In[17]:


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


# # HB-Search

# In[18]:


def changeDirectory(programDirectory: str = "."):
    """
    Changes working directory to the folder containing the needed programs. 
    Path to program folder needs to be entered if started from other location
    :param programDirectory: Path to the program folder. Standard: Set to current working directroy
    """
    
    #Changes working directory to the entered path adjusted to the operating system used.
    cmd.cd(os.path.normpath(programDirectory))
    #os.chdir(os.path.normpath(programDirectory)) #can be deleted?
    # Usefull if we want to give an error! To tell the person in which directory they are located!
    #cwd = os.getcwd()


# In[19]:


def useObject(input_molecule: str):
    """
    Enables to use a chosen object of the pymol session for HBsearch/HBnetwork.
    Saves the chosen object as pdb-file in the pdb_files folder of the program folder
    under the name of the object.
    :param input_molecule: pymol object name of the structure that should be used for the HBsearch/HBnetwork run
    """
    
    #Saves pymol object by its name in the pdb_files folder.
    cmd.save(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"), input_molecule)


# In[20]:


def removeObject(input_molecule: str):
    """
    Removes pdb-file of chosen object name. Used after HBsearch run for deleting unwanted data.
    :param input_molecule: Name of the molecule to be removed.
    """
    
    #Removes pdb-file from pdb_files ordner by entered name.
    os.remove(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"))


# In[21]:


def fetchPDB(pdbID: str, object_name: str = ""):
    """
    Fetched PDB-file by PDB-ID to pymol session and sets object name in pymol session.
    PDB-file is saved under the PDB-ID in pdb_files folder.
    If PDB-file is already existant with same name. No new file will be fetched.
    Internet conncetion required.
    :param pdbID: PDB-ID of protein structure to be fetched and used for HBsearch/HBnetwork run.
    :param object_name: Name the pdb file shoudl be displayed as in the pymol session.
    """
    
    #Checks if a specific name is chosen. If no name is chosen, file will be named after the PDB-ID.
    if object_name == "":
        object_name = pdbID
        
    #setting fetch_path to desired folder
    cmd.set("fetch_path", os.path.normpath("./pdb_files/"))
    
    #fetching pdb file if not in folder
    cmd.fetch(pdbID, name = object_name, type = "pdb")


# In[22]:


def startHBsearch(molecule: str, hb_file: str, solvent_key:str, pse_file:str, connections: str) -> str:
    """
    Starts the HBsearch run. HBsearch finds hydrogen bonds based on a set of parameters based on the chemical nature
    of the atoms of the protein/cofactor/nucleic acid (maximal Van-der-Waals radii and maximal covalent radii) and their distances to each other.
    Possible pairs forming hydrogen bonds are listed in the hb-define file. Chemical nature of atoms are listed in period-table.info file as standard.
    :param molecule: Strucutre used for HBsearch run
    :param hb_file: HB-file sued to define possible hydrogen bond interactions. Standard set to hb-define.txt file
    :param solvent_key: If hydrogen bond bridges with solvent should be considered. Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.
    :param pse_file: File containing the chemical nature of the atoms. Standard set to period-table.info. To create own one see standard file for structure.
    :param connections: If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, but could be possible due to rotation of the residues. If connections = 1: Special conncetions will be taken into account. Standard set to 0: No special connectoins will be taken into account.
    :return hbs_output: Returns string with the output of the HBsearch run.
    """

    # Setting environment variable. Sets chemical nature of atoms.
    os.environ['PSE_FILE'] = pse_file
    # Determine operation system to start correct HBsearch application.
    system = platform.system()
    # Executing hb_search with chosen parameters
    hbs_output = subprocess.run(os.path.normpath(f"./{system}/hb-search -hb {hb_file} -solv {solvent_key} -con {connections} ./pdb_files/{molecule}.pdb"), capture_output=True, shell=True, check = True, text = True).stdout
    
    return hbs_output


# In[23]:


def readInHBS(hbs_output: str) -> pd.DataFrame():
    """
    Reads in HBsearch output and converts it to pandas dataframe for further processing.
    :param hbs_output: HBsearch output as string:
    :return df: Returns dataframe containing sorted and filtered HBsearch results
    """
    
    #Creates a list with each line of the HBsearch output being assigned to an entry in the list.
    hbs_rows = [i for i in hbs_output.split('\n')]
    #Splits the entries in the row-list further, creating an array with the first order list containing the rows and the second order list separating by column.
    hbs_split = [i.split() for i in hbs_rows]
    
    #Column names for the later created dataframe containing the HBsearch output
    HEAD_LST = ['IDENT','ACC','sep1','DONO',':','x','y','z','sep2','a','b']
    
    #Creating pandas dataframe by reading in the primitive dataframe based on lists from the HBsearch output and setting column names.
    df = pd.DataFrame(hbs_split, columns = HEAD_LST)
    #Filters dataframe, so only entries containing the partners of a hydrogen bond are listed in the dataframe.
    df = df[df["IDENT"] == "HBOND"]
    
    return df


# In[24]:


def prepareLists(dataframe: pd.DataFrame) -> List:
    """
    Extracts the acceptors and donors from the HBsearch/HBnetwork output into separate lists.
    :param dataframe: Dataframe containing the HBsearch/HBnetwork output splitted into columns.
    :return acceptor: Returns list with acceptor atoms in hydrogen bond with complementary indexes to their partner in the donor list.
    :return donor: Returns list with donor atoms in hydrogen bond with complemenary indexes to their partner in the acceptor list.
    """
    
    #Extracting important acceptor and donor entries from dataframe containing HBsearch results as separate lists.
    acceptor_pre = list(dataframe["ACC"])
    donor_pre = list(dataframe ["DONO"])
    
    #Creating lists empty lists for further use. Lists should contain later entries from acceptor and donor entries, respectively, separated by the chain, residue, residue ID, and atom. 
    acceptor = []
    donor = []
    
    #Separates entries for acceptors in HBsearch output by chain, residue, residue ID, and atom for faciliated reorganisation for PyMol input.
    for i in range(len(acceptor_pre)):#only one list can be used for iteration of list index, since each acceptor needs a donor to form a hydrogen bond. So acceptor and donor list need to have same length.
        acceptor.append(tuple(acceptor_pre[i].split(":")))
    #Separates entries for donors in HBsearch output by chain, residue, residue ID, and atom for faciliated reorganisation for PyMol input.    
    for j in range(len(donor_pre)):
        donor.append(tuple(donor_pre[j].split(":")))

    return acceptor, donor


# In[25]:


def displayDistances(acceptor: List, donor: List, object_name: str, run_information:str) -> None:
    """
    Displays the hydrogen bonds of the hydrogen bond acceptors and their respective 
    donors found by the HBsearch run as distances without labeling in PyMol.
    :param acceptor: Lists of acceptor atoms with complementary index to their respective donors in donor list. Entry tuples contain chain, residue, residue ID and atom.
    :param donor:  Lists of donor atoms with complementary index to their respective donors in acceptor list. Entry tuples contain chain, residue, residue ID and atom.
    :param object_name: Name of program this function is used in.
    :param run_information: Displays in the selection of sticks, which program was used. In case of HBnetwork, additionally the search query is shown.
    """
    
    #Creating list of Hbonds for grouping. Faciliates displaying
    bondList = []
    
    #Displays connection between the acceptor atom and the respective donor atom in pymol. Therfore, creates distance line object
    for i in range(len(acceptor)):        
        cmd.distance(f"{run_information}_{object_name}_hydrogenBond_{i}", #Name of the distance line object
                     f"{object_name}//{acceptor[i][0]}/{acceptor[i][1]}/{acceptor[i][3]}", #Acceptor molecule. Tuple entries of acceptor list ordered by PyMol selection format.
                     f"{object_name}//{donor[i][0]}/{donor[i][1]}/{donor[i][3]}")        #Donor molecule. Tuple entries of acceptor list ordered by PyMol selection format.
        
        #Creates a list containing each distance object. Faciliates further adjustments of the distance line objects. 
        bondList.append(f"{run_information}_{object_name}_hydrogenBond_{i}")
        
    #Groups all distance line objects into one cluster for better clarity in PyMol window.
    cmd.group(f"{run_information}_{object_name}_hydrogenBonds", " ".join(bondList))
    #Hides distance labels (distance in Angstrom) for better clarity.
    cmd.hide("labels", f"{run_information}_{object_name}_HydrogenBonds")


# In[26]:


def showSticks(acceptor: List, donor: List, object_name: str, run_information: str):
    """
    Displays the residues participating in hydrogen bonds as sticks in PyMol.
    :param acceptor: Lists of acceptor atoms with complementary index to their respective donors in donor list. Entry tuples contain chain, residue, residue ID and atom.
    :param donor: Lists of donor atoms with complementary index to their respective donors in acceptor list. Entry tuples contain chain, residue, residue ID and atom.
    :param object_name: Name of the object the HBsearch/HBnetwork run was performed on.
    :param run_information: Displays in the selection of sticks, which program was used. In case of HBnetwork, additionally the search query is shown.
    """
    
    #Creating list of residues that should be displayed as sticks. Improves speed of the program by grouping and showing sticks of group.
    stickList = []
    
    #Creating list with all residues participating in hydrogen bonds. Entries are converted to PyMol selection format.
    for i in range(len(acceptor)):##maybe possible to modify with zip(); Iterating over all entries
        stickList.append(f"/{object_name}//{acceptor[i][0]}/{acceptor[i][1]}") #Appending to stick-list all acceptor residues. Residue information is converted to PyMol selection format.
        stickList.append(f"/{object_name}//{donor[i][0]}/{donor[i][1]}") #Appending to stick-list all donor residues. Residue information is converted to PyMol selection format.
    
    #Selecting all sticks by selecting all entries in the sticklist
    cmd.select(f"{run_information}_Connections_Sticks_{object_name}", " ".join(stickList)) #For selection list entries need to be converted to string, with entries separeated by space character.
    cmd.show("sticks", f"{run_information}_Connections_Sticks_{object_name}") #Selection is displayed as sticks.
    cmd.deselect() #Selection is deselected for better clarity and to spare the user deselecting selection by him-/herself.


# In[27]:


def hbsearch(molecule:str, molecule_name: str = "", directory:str = ".", 
         use_object: str = "0", remove_object: str = "1", hb_file: str = "hb-define.txt", 
         solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", connections: str = "0"):
    """
    Runs HBsearch using a given biomolecule (either fetched by the PDB-ID or using an object of the PyMol session)
    and displays the hydrogen bonds occuring in the molecule.
    :param molecule: PDB-ID of molecule HBsearch run should be performed on or in case an object in the PyMol session has to be used: The name of the object.
    :param molecule_name: If object fetched by PDB-ID. Standard set to blank.
    :param directory: Directory of the program folder. Directes working directory in PyMol to entered directory. Obsolete when script is started from program folder and HBsearch is started without changing the directory. Standard set to current working directory (".").
    :param use_object: If an object in the PyMol session has to be used for the HBsearch run. 0: No; 1: Yes. Standard set to 0. If set to 0. Molecule will be fetched by PDB-ID.
    :param remove_object: If created PDB-file in pdb_files folder should be deleted after HBsearch run. 0: No; 1: Yes. Standard set to 0.
    :param hb_file: HB-file sued to define possible hydrogen bond interactions. Standard set to hb-define.txt file
    :param solvent_key: If hydrogen bond bridges with solvent should be considered. Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.
    :param pse_file: File containing the chemical nature of the atoms. Standard set to period-table.info. To create own one see standard file for structure.
    :param connections: If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, but could be possible due to rotation of the residues. If connections = 1: Special conncetions will be taken into account. Standard set to 0: No special connectoins will be taken into account.
    """
    
    changeDirectory(directory) #Change directory to program folder directory. Needed PyMol directory is not set to program folder for HBsearch run
    
    #Checks if user wants to use own object in PyMol session or wants to fetch a protein structure from the PDB
    if use_object == "0": #User wants to fetch a protein from the PDB using a PDB-ID
        fetchPDB(molecule, molecule_name) #Uses fetch command to fetch PDB_ID and allows user to name the fetched object in PyMol session.
    elif use_object == "1": #User wants to use own object.
        useObject(molecule) #Object with entered name is saved as pdb-file in the pdb_files folder and is used for the following HBsearch run
    
    
    hbs_output = startHBsearch(molecule, hb_file, solvent_key, pse_file, connections) #Starts an HBsearch run with given parameters. Output is stored
    hbs_dataframe = readInHBS(hbs_output) #Converts HBsearch run to dataframe
    acceptor, donor = prepareLists(hbs_dataframe) #Prepares acceptor and donor lists for displaying in PyMol
    
    run_information = "HBsearch" #Which program is run. Needed for displayDistances() and showSticks()
    
    if molecule_name == "": #Checking if custom molecule name was entered: If no name was entered, PDB-ID or saved object name is used.
        displayDistances(acceptor, donor, molecule, run_information) #Displays distances of hydrogen bond acceptors with their respective donors as distance objects without label in PyMol and groups distance objects according to the strucutre object they are based on.
        showSticks(acceptor,donor, molecule, run_information) #Shows residues participating in hydrogen bonds as sticks
    else: #When molecule name was entered, PyMol strucutre object posses molecule name. So this is used for following PyMol selection based commands.
        displayDistances(acceptor, donor, molecule_name, run_information) #Displays distances of hydrogen bond acceptors with their respective donors as distance objects without label in PyMol and groups distance objects according to the strucutre object they are based on.
        showSticks(acceptor,donor, molecule_name, run_information) #Shows residues participating in hydrogen bonds as sticks
    
    if remove_object == "1": #Checks if parameter remove_object is set to 1. If yes: created PDB-file in pdb_files folder is deleted.
        removeObject(molecule)


# In[28]:


#Creation of command in PyMol.
cmd.extend("hbsearch", hbsearch) #When read in in PyMol the script creates a command in PyMol which can be started via the PyMol command line.


# # HB-Network - Initialization

# In[29]:


class saveBot:
    def __init_():
        directory_name_save:str = None
        dictionary_save: Dict = None
        molecule_name_save: str = None

save = saveBot()


# In[30]:


def createDirectory(molecule_name:str) -> str:
    """
    Creates new directory, containing HBnetwork files later on. directory is created in program folder in ./HB_network/.
    Directory name: MoleculeName_Date_Time.
    :param molecule_name: Name of the molecule HBnetwork is run on. 
    :return directory_name: Returns name of created folder in HB_network folder. Needed for guiding HBnetwork.
    """
    #Get current time and format it
    current_time = datetime.now()
    formated_time = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    
    #Set folder name in HB_network directory
    directory_name = f"{molecule_name}_{formated_time}"
    
    #Create directory with set folder name in HB_network folder
    os.mkdir(os.path.normpath(f"./HB_network/{directory_name}"))
    
    #Informs user, where cluster files will be saved.
    print(f"HBnetwork run will be saved in", os.path.normpath(f"HB_network/{directory_name}"))
    
    return directory_name


# In[31]:


def createHBnetwork(molecule: str, directory_name: str, hb_file: str = "hb-define.txt", 
                    solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", 
                    connections: str = "0"):
    """
    Runs HBnetwork on the entered molecule. Files containing information about the hydrogen bond
    clusters are created in the HB_network folder in the program folder in the directory.
    Final directory is termed: MoleculeName_Date_Time. Returns summary of all hydrogen bond clusters
    :param molecule: Strucutre used for HBsearch run. HBnetwork is run on the same molecule using the output of HBsearch.
    :param directory_name: Saving folder for the HBnetwork run in HB_network subfolder in the program folder. 
    :param hb_file: HB-file sued to define possible hydrogen bond interactions. Standard set to hb-define.txt file
    :param solvent_key: If hydrogen bond bridges with solvent should be considered. Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.
    :param pse_file: File containing the chemical nature of the atoms. Standard set to period-table.info. To create own one see standard file for structure.
    :param connections: If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, but could be possible due to rotation of the residues. If connections = 1: Special conncetions will be taken into account. Standard set to 0: No special connectoins will be taken into account.
    :return hbn_output: Summary of all hydrogen bond clusters. Used later for indexing.
    """
    
    #Returns operating system the user utilizes. Needed to guide HBsearch/HBnetwork to the correct folders containing appropriate executables.
    system = platform.system()
    
    #Creating a .hb file for the HBnetwork program based. HBsearch run is needed as parameter for HBnetwork run.
    hbs_output = startHBsearch(molecule, hb_file, solvent_key, pse_file, connections) #Gathering HBsearch output

    with open(os.path.normpath(f"./HB_network/{directory_name}/{molecule}.hb"),"w") as fh: #writing HBsearch output in textfile in correct directory
        fh.write(hbs_output)
        
    #Runs HBnetwork. Creates output files in HB_network/molecule_date_time/ 
    cluster_dir = os.path.normpath(f"./HB_network/{directory_name}")
    hbnetwork_dir = os.path.normpath(f"../../{system}/hb-network {molecule}.hb") #Program_dir hb_file_dir
    hbn_output = subprocess.run(hbnetwork_dir, cwd = cluster_dir, capture_output=True, shell=True, check = True, text = True).stdout #Runs HBnetwork. Cluster files are created. Summary is given as output
    
    return hbn_output


# In[32]:


def cleanHBnetwork(directory_name: str):
    """
    Cleans HBnetwork CLUSTER directory of the run, by deleting all files with a size of 0 bytes.
    :param directory_name: Name of the directory the HBnetwork cluster files are saved in.
    """
    
    ##Deletes not needed files without any network
    ##Give later out: If file is not existant than: "No network found"
    ##Giving networks with aminoacids and with atoms
    
    directory = os.path.normpath(f"./HB_network/{directory_name}/CLUSTER") #By HBnetwork created CLUSTER directory is targeted
    
    #Files with a size of 0 bytes are deleted since they dont contain clusters.
    for file in os.listdir(directory): #Iterating over the directory 
        
        file_dir = os.path.normpath(f"{directory}/{file}")
        
        if os.path.getsize(file_dir) == 0: #Files with a size of 0 bytes are removed
            os.remove(file_dir)          
    


# In[33]:


def indexHbnetwork(hbn_output: str):
    """
    Creates a dictionary with each atom participating in hydrogen bonds. Atoms are set as keys, the respective cluster
    is set as value. Needed for later searching of the correct hydrogen bond cluster containg the hydrogen bond network
    the selected atom participates. 
    :param hbn_output: Output of HBnetwork run. Contains summary of all clusters, containing each participating atom uniquly in a cluster.
    :return cluster_dict: Returns a dictionary containg the clusters as values and the atoms in it as unique keys. Used as index dictionary for search of correct file later on.
    """
    #Only using the last part of the HBnetwork output: The cluster summary
    starting_index = hbn_output.find("Cluster") #Finding the start of the summary
    hbn_sorted = hbn_output[starting_index:].split("\n") #Slicing the output string and splitting the resulting summary by lines
    
    #Start index for sorting the splitted output by cluster entries
    start = 0 
    
    #Creating list and dictionary for indexed HBnetwork summary
    cluster_list = []
    cluster_dict = {}
    
    #Creating a list containing each cluster as an separate entry
    for i in range(len(hbn_sorted)):
        hbn_sorted[i] = hbn_sorted[i].strip() #Stripping the entries, since otherwise \t might interefere later on
        if hbn_sorted[i] == "END": #END as marker that respective cluster entry ended
            cluster_list.append(hbn_sorted[start:i]) #Appending the cluster entry without the "END"
            start = i + 1 #Setting new start index one later than the last "END"
            
    #Creating dictionary with the cluster containing the selected atom as a value and the atom as the unique key.
    for j in range(len(cluster_list)): #Iterating over all cluster entries
        for k in range(1,len(cluster_list[j])): #Iterating over all atoms in that entry. Index 0 is spared out, since it contains the cluster number
            
            #Removing residue name from entry for faciliated data selection and processing in further steps
            cluster_atom_list = cluster_list[j][k].split(":") #Splitting the string to a list
            cluster_atom_list.pop(2) #Removing the residue three letter code from the entry
            cluster_atom_key = "/".join(cluster_atom_list) #Joining it again with format as it will be used later on, oriented on PyMol selection format
            
            #Creating variable for dictionary value containing the cluster in which atom (key) participates.
            #Processes cluster name to correct format, faciliating later use.
            #" " is replaced with "-". "Cluster" is lowered to "cluster". ":" at the end is cut off
            #Final format: "cluster-1"
            cluster_number = cluster_list[j][0].replace(" ","-").lower()[:-1]
            
            #Filling dictionary with atom as key and cluster as value
            cluster_dict[cluster_atom_key] = cluster_number #Setting the cluster as value to the atom as the key.
    
    return cluster_dict


# In[34]:


def initiateHBnetwork(molecule:str, molecule_name = "", directory:str = ".", 
         use_object: str = "0", remove_object = "0", hb_file: str = "hb-define.txt", 
         solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", connections: str = "0"):
    """
    
    """
    
    changeDirectory(directory) #Change directory to program folder directory. Needed PyMol directory is not set to program folder for HBnetwork run
    
    #Checks if user wants to use own object in PyMol session or wants to fetch a protein structure from the PDB
    if use_object == "0": #User wants to fetch a protein from the PDB using a PDB-ID
        fetchPDB(molecule, molecule_name) #Uses fetch command to fetch PDB_ID and allows user to name the fetched object in PyMol session.
    elif use_object == "1": #User wants to use own object.
        useObject(molecule) #Object with entered name is saved as pdb-file in the pdb_files folder and is used for the following HBsearch run
    
    #Creating directory for HB_network run, where cluster files are going to be saved.
    #Directory_name variable is saved in object "save", since it will be needed later on and program does not have to be rerun every time.
    directory_name = createDirectory(molecule)
    
    save.directory_name_save = directory_name #Saves the directory name in the object "save" for further use in the next steps.
    
    #Running HBsearch, since output is needed for HBnetwork run.
    hbn_output = createHBnetwork(molecule, directory_name)
    
    #########Cleans HB_network folder from empty cluster files. Improving search efficiency later on.
    #Could also be deleted, since not necessairy for function.
    #cleanHBnetwork(directory_name)
    
    #Created cluster dictionary serving as an indexing dictionary later on
    #Most important step for initalization.
    hbn_cluster_dict = indexHbnetwork(hbn_output)
    
    save.dictionary_save = hbn_cluster_dict #Indexing dictionary is saved to object "save" for furhter use later on.
                                            #True ouput of function, since this function should be called only once for each HBnetwork.
                                            #Searches should be possible without reoccuring initation runs. Therefore saving output in object.
    
    #
    if molecule_name == "": #Checking if custom molecule name was entered: If no name was entered, PDB-ID or saved object name is used.
        save.molecule_name_save = molecule #
    else: #When molecule name was entered, PyMol strucutre object posses molecule name. So this is used for following PyMol selection based commands.
        save.molecule_name_save = molecule_name #
    
    
    
    
    #Checks if parameter remove_object is set to 1. If yes: created PDB-file in pdb_files folder is deleted.
    if remove_object == "1":
        removeObject(molecule)  


# In[35]:


initiateHBnetwork("4akr")


# In[36]:


cmd.extend("initiateHBnetwork", initiateHBnetwork)


# # HB-Network PyMol-Display

# In[37]:


def readoutHBnetwork(query: str, checkFor: str = "RESIDUE" ) -> str:
    
    #Checks if only hydrogen network for atom or whole residue should be shown
    #Reads out dictionary, which Cluster files to check
    #How should the user input b Chain:ResIDAtom. Best would be PyMol input:
    #Last "/" important for closing arguement. otherwise other atoms could be chosen. E.g. 3 or 31 if A/3 and not A/3/
    #If List empty: Resdue does not participate in any cluster
    
    #Sets the index dictionary from the object "save" from the previous initiation run as current index dictionary for the function.
    dictionary = save.dictionary_save
    
    #Creates an empty list, which contains later on the names of the clusters containing the atom/residue.
    destination_cluster_list = []
    
    #Checks if the clusters of single atom should be checked or if the clusters the whole residue participates should be displayed.
    if checkFor == "ATOM":
        #Appends the cluster name with the searched atom to the cluster list.
        #Appends cluster with a matching atom in the format "chain/residue ID/atom"
        destination_cluster_list.append(dictionary[query])
        
        return destination_cluster_list #Cluster list is passed to next function for searching respective cluster files.
                                        #Returns blank list, if no cluster, the residue participates, were found.
    
    #Checks if all clusters, where any atom of the residue participate, should be processed.
    elif checkFor =="RESIDUE":
        #Appends all cluster matching the query input of "chain/residue/"
        [destination_cluster_list.append(value) for key, value in dictionary.items() if query in key]
        
        return destination_cluster_list #Cluster list is passed to next function for searching respective cluster files.
                                        #Returns blank list, if no cluster, the residue participates, were found.


# In[38]:


def readInHBN(cluster_file_output: str) -> pd.DataFrame:
    """
    
    """
    
    
    #Creates a list with each line of the HBnetwork cluster file output being assigned to an entry in the list.
    cluster_rows = [i for i in cluster_file_output.split('\n')]
    #Splits the entries in the row-list further, creating an array with the first order list containing the rows and the second order list separating by column.
    cluster_split = [i.split() for i in cluster_rows]
    
    #Column names for the later created dataframe containing the HBnetwork cluster files output.
    HEAD_LST_CLUSTER = ['IDENT','ACC','sep1','DONO']
    
    #Creating pandas dataframe by reading in the primitive dataframe based on lists from the HBsearch output and setting column names.
    df_cluster = pd.DataFrame(cluster_split, columns = HEAD_LST_CLUSTER)
    
    return df_cluster


# In[39]:


def prepareDataFrameHBnetwork(destination_cluster_list: str) -> pd.DataFrame:

    #Get directory from global variable
    #Reads out files with given file names
    #Creates dataframe for each file -->readInHBN
    #Append dataframe acceptor and donor to repective lists
    #Return lists for displaying in Pymol as sticks an connected with distances
    #Create readInHBN. Orient on readInHBS
    #Iteration through files see cleanHBnetwork
    #Create acceptor and donor lists so that prepareLists, showDistances and showSticks work on it. Less work
    #--> Therfore append everything to an acceptor and donor list returning to prepareList
    #--> When no dataframe is empyt --> Atom /Residue does not participate in any hydrogen bonds
    
    directory = save.directory_name_save
    
    #Directory, where cluster files for initiation run are saved
    cluster_directory = os.path.normpath(f"./HB_network/{directory}/CLUSTER") #By HBnetwork created CLUSTER directory is targeted
    
    #Creating empty DataFrame. All dataframes created by reading out the generated cluster files is appended to this dataframe
    hBond_cluster_dataframe = pd.DataFrame()
    
    #
    for file in destination_cluster_list: #Iterating over given clusters.
        
        #Save path to cluster file for reading in.
        cluster_dir = os.path.normpath(f"{cluster_directory}/{file}.hb-ntw")
        
        #Reading in cluster file.
        with open(cluster_dir,"r") as cluster_file:
            cluster_output = cluster_file.read()
            
        cluster_output = cluster_output.strip()
        
        #Processing cluster file output to dataframe
        if cluster_output:
            new_cluster_df = readInHBN(cluster_output)
        else:
            continue
            
        hBond_cluster_dataframe = pd.concat([hBond_cluster_dataframe, new_cluster_df], ignore_index = True) 
    return hBond_cluster_dataframe


# In[42]:


def showNetwork(query: str, checkFor = "RESIDUE"):
    pass
    ####### Learn objects. Save data like dictionary, directory name and co. as objects -->Classes, Objects
    destination_cluster_list = readoutHBnetwork(query, checkFor)
    
    hbn_dataframe = prepareDataFrameHBnetwork(destination_cluster_list)
    
    if hbn_dataframe.empty:
        print(f"This {checkFor} does not participate in any hydrogen bonds")
        return None
    
    else:
        acceptor, donor = prepareLists(hbn_dataframe)
        print(acceptor)
        print(donor)
    #
    molecule_selection = save.molecule_name_save
        
    #
    run_information = f"HB_network"
    print(run_information)
    displayDistances(acceptor, donor, molecule_selection, run_information)
    showSticks(acceptor, donor, molecule_selection, run_information)


showNetwork("A/9/")


# In[ ]:


cmd.extend("showNetwork", showNetwork)


# In[ ]:


print(os.getcwd())


# In[ ]:


"""<!--       _
       .__(.)< (MEOW)
        \___)   
 ~~~~~~~~~~~~~~~~~~-->"""


# In[ ]:





# In[ ]:




