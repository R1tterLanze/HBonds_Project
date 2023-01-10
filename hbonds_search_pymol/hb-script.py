#!/usr/bin/env python
# coding: utf-8

# # GUI hbond

# ### Initialization

# In[38]:


import pandas as pd
import subprocess
import os
from pymol import cmd, stored
import platform
from typing import List, Dict
from datetime import datetime


# ### Functions

# # HB-search

# ## Error Classes

# In[39]:


class DirectoryError(Exception):
    
    def __init__(self, *args):

        if args:
            self.message = args[0]
        else:
            self.message = None
    
    def __str__(self):

        if self.message:
            return self.message
        else:
            return "DirectoryError has been raised"
        


# In[40]:


class SelectionError(Exception):
    
    def __init__(self, *args):

        if args:
            self.message = args[0]
        else:
            self.message = None
    
    def __str__(self):

        if self.message:
            return self.message
        else:
            return "SelectionError has been raised"
        


# In[41]:


class ProcessError(Exception):
    
    def __init__(self, *args):

        if args:
            self.message = args[0]
        else:
            self.message = None
    
    def __str__(self):

        if self.message:
            return self.message
        else:
            return "ProcessError has been raised"
        


# In[42]:


class InputError(Exception):
    
    def __init__(self, *args):

        if args:
            self.message = args[0]
        else:
            self.message = None
    
    def __str__(self):

        if self.message:
            return self.message
        else:
            return "Input has been raised"
        


# In[43]:


def changeDirectory(programDirectory: str = "."):
    #Check if directory is the right one. Like check if the needed programs are there and the needed folders:
    #Needed folders: Darwin, HBnetwork, Linux,pdb_files, Windows, hb-define.txt, hb-search_script.py, 
    #period-table-info.txt, README.md
    #If not in correct folder, say to direct into correct folder
    #Tell user in which directory he is
    """
    Changes working directory to the folder containing the needed programs. 
    Path to program folder needs to be entered if started from other location
    :param programDirectory: Path to the program folder. Standard: Set to current working directroy
    """
    
    #Checking if all files/folder for HBsearch and HBnetwork run are present
    #Here: Defining which files are of utmost necessity
    NEEDED_FOLDERS = ["Darwin", "HB_network", "Linux", "pdb_files", "Windows"]
    NEEDED_FILES = ["hb-define.txt", "period-table-info.txt"]
       
    path = os.path.normpath(programDirectory) #Path converted suitable to current operating system
    
    if not os.path.exists(path): #Checks if entered path does exist.
        raise DirectoryError("Entered path does not exists")
    
    #List to be iterated over for file/directory existance checking
    required_items = NEEDED_FILES + NEEDED_FOLDERS
    #Checking if all required folders/files are present.
    #######Program files like hb-search will be checked for their existance later on.
    for i in required_items:
        if not os.path.exists(os.path.normpath(f"{path}/{i}")):
            raise DirectoryError(f"Wrong directory or folders/files are missing/renamed. You are currently in {os.getcwd()}.\nPlease change to correct directory by using the \"directory\" attribute when running the program.\nExample: hbsearch 4awn, directory = Path")
    
   
    #Changes working directory to the entered path adjusted to the operating system used.
    cmd.cd(path)


# In[44]:


def useObject(input_molecule: str):
    #Say when the object user wants to use does not exist.
    """
    Enables to use a chosen object of the pymol session for HBsearch/HBnetwork.
    Saves the chosen object as pdb-file in the pdb_files folder of the program folder
    under the name of the object.
    :param input_molecule: pymol object name of the structure that should be used for the HBsearch/HBnetwork run
    """
    #Get list of all objects of the current PyMol session
    object_list = cmd.get_object_list("(all)")
    #Check if the object to be used exists.
    if input_molecule not in object_list:
        raise SelectionError(f"The object \"{input_molecule}\" does not exist.")
        
    #Saves pymol object by its name in the pdb_files folder.
    cmd.save(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"), input_molecule)


# In[45]:


def removeObject(input_molecule: str):
    #Say when there is no such object to be removed.
    """
    Removes pdb-file of chosen object name. Used after HBsearch run for deleting unwanted data.
    :param input_molecule: Name of the molecule to be removed.
    """
    #Dont really need the following. Otherwise the whole program wouldnt work. But maybe someone manages to get this far so...
    #Get list of all objects of the current PyMol session
    object_list = cmd.get_object_list("(all)")
    #Check if the object to be used exists.
    if input_molecule not in object_list:
        raise SelectionError(f"The object \"{input_molecule}\" does not exist and cant therefore be removed")
    
    #Removes pdb-file from pdb_files ordner by entered name.
    os.remove(os.path.normpath(f"./pdb_files/{input_molecule}.pdb"))


# In[46]:


def fetchPDB(pdbID: str, object_name: str = ""):
    #Error, when PDB-ID is nonexistant.
    #Error when objectname already exists.
    """
    Fetched PDB-file by PDB-ID to pymol session and sets object name in pymol session.
    PDB-file is saved under the PDB-ID in pdb_files folder.
    If PDB-file is already existant with same name. No new file will be fetched.
    Internet conncetion required.
    :param pdbID: PDB-ID of protein structure to be fetched and used for HBsearch/HBnetwork run.
    :param object_name: Name the pdb file shoudl be displayed as in the pymol session.
    """
    
    
    #Get list of all objects of the current PyMol session
    object_list = cmd.get_object_list("(all)")
    #Check if the object name to be used already exists.
    if object_name in object_list:
        raise SelectionError(f"The object \"{object_name}\" already exists.")
    
    #Checks if a specific name is chosen. If no name is chosen, file will be named after the PDB-ID.
    #If object name with PDB ID as name already exists PyMol handles it. 
    if object_name == "":
        object_name = pdbID
        
    #setting fetch_path to desired folder
    cmd.set("fetch_path", os.path.normpath("./pdb_files/"))
    
    #fetching pdb file if not in folder
    cmd.fetch(pdbID, name = object_name, type = "pdb")


# In[47]:


def startHBsearch(molecule: str, hb_file: str, solvent_key:str, pse_file:str, connections: str) -> str:
    #Raise Error, when gives out zeroNullError or something like that. Most likely program cant handle that molecule.
    #Raise error, when file could not be found --> Probably in wrong directory.
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
    
    #Paths to input files
    path_PDB = os.path.normpath(f"./pdb_files/{molecule}.pdb")
    path_hb_file = os.path.normpath(f"./{hb_file}")
    path_pse_file = os.path.normpath(f"./{pse_file}")
    #Putting all paths into one list for iterating
    paths = [path_PDB, path_hb_file, path_pse_file]
    for i in paths:
        if not os.path.exists(i):
            raise FileNotFoundError(f"File \"{i}\" does not exist")
    
    try:
        hbs_output = subprocess.run(os.path.normpath(f"./{system}/hb-search -hb {hb_file} -solv {solvent_key} -con {connections} ./pdb_files/{molecule}.pdb"), capture_output=True, shell=True, check = True, text = True).stdout
    except subprocess.CalledProcessError:
        raise ProcessError("hb-search could not be run. Try with PDB-ID \"4AWN\" to make sure program works properly.")
    
    return hbs_output


# In[48]:


def readInHBS(hbs_output: str) -> pd.DataFrame():
    #Here could also be errors, but dont know which. Wait for the testing. 
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


# In[49]:


def prepareLists(dataframe: pd.DataFrame, water_connections: str) -> List:
    #Dont know what could go wrong here.
    """
    Extracts the acceptors and donors from the HBsearch/HBnetwork output into separate lists.
    :param dataframe: Dataframe containing the HBsearch/HBnetwork output splitted into columns.
    :param water_connections: Define if hydrogen bridges with water molecules involved should be considered. "0" for NO, "1" for YES.
    :return acceptor: Returns list with acceptor atoms in hydrogen bond with complementary indexes to their partner in the donor list.
    :return donor: Returns list with donor atoms in hydrogen bond with complemenary indexes to their partner in the acceptor list.
    """
    
    #Extracting important acceptor and donor entries from dataframe containing HBsearch results as separate lists.
    acceptor_pre = list(dataframe["ACC"])
    donor_pre = list(dataframe ["DONO"])

    #Creating lists empty lists for further use. Lists should contain later entries from acceptor and donor entries, respectively, separated by the chain, residue, residue ID, and atom. 
    acceptor = []
    donor = []
    
    #Entries that can occur when solvent key is not set to "NONE". Have to be filtered out, since no specific atom is participating in Hbond.
    solvent_connections_list = ["SOLV", "POS", "NEG"]
    
    #Separates entries for acceptors in HBsearch output by chain, residue, residue ID, and atom for faciliated reorganisation for PyMol input.
    for i in range(len(acceptor_pre)):#only one list can be used for iteration of list index, since each acceptor needs a donor to form a hydrogen bond. So acceptor and donor list need to have same length.
        
        #When special solvent keys other than "NONE" are used it may occur to hbonds between non defined atoms. Here they are filtered out
        if acceptor_pre[i] in solvent_connections_list:
            continue
        
        #Appending processed atom selection to acceptor list
        acceptor.append(tuple(acceptor_pre[i].split(":")))
        
    #Separates entries for donors in HBsearch output by chain, residue, residue ID, and atom for faciliated reorganisation for PyMol input.    
    for j in range(len(donor_pre)):
        
        #Separates entries for donors in HBsearch output by chain, residue, residue ID, and atom for faciliated reorganisation for PyMol input.
        if donor_pre[j] in solvent_connections_list:
            continue
        
        #Appending processed atom selection to donor list    
        donor.append(tuple(donor_pre[j].split(":")))
    
    #Deleting entries with hydrogen bonds to water molecules.
    if water_connections == "0":
        #Creating list with indexes to be popped
        pop_list = []
        
        #Append entry indexes containing water to pop_list
        for k in range(len(acceptor)):
            if acceptor[k][2] == "HOH":
                pop_list.append(k)
                
            if donor[k][2] == "HOH":
                pop_list.append(k)
        
        #Removing duplicates from pop_list
        pop_list = list(dict.fromkeys(pop_list)) #Remove duplicates from list
        
        #Deleting entries containing water.
        for index in sorted(pop_list, reverse = True): #Reverse list is used to obtain order of indexes.
            acceptor.pop(index)
            donor.pop(index)
        
    return acceptor, donor


# In[50]:


def displayDistances(acceptor: List, donor: List, object_name: str, run_information:str) -> None:
    #Maybe object already exists error or something like that. should be all working. Dont know what could possibly go wrong.
    """
    Displays the hydrogen bonds of the hydrogen bond acceptors and their respective 
    donors found by the HBsearch run as distances without labeling in PyMol.
    :param acceptor: Lists of acceptor atoms with complementary index to their respective donors in donor list. Entry tuples contain chain, residue, residue ID and atom.
    :param donor:  Lists of donor atoms with complementary index to their respective donors in acceptor list. Entry tuples contain chain, residue, residue ID and atom.
    :param object_name: Name of object this function is used on.
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


# In[51]:


def showSticks(acceptor: List, donor: List, object_name: str, run_information: str):
    #Should also be working without exceptions. lets wait for the testing.
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


# In[52]:


def hbsearch(molecule:str, molecule_name: str = "", directory:str = ".", 
         use_object: str = "0", remove_object: str = "0", water_connections:str = "0", hb_file: str = "hb-define.txt", 
         solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", connections: str = "0"):
    #use string.casefold() or string.casefold().upper() to ensure correct inputs.
    #Check for correct inputs and form. Guide the user by the hand. 
    """
    Runs HBsearch using a given biomolecule (either fetched by the PDB-ID or using an object of the PyMol session)
    and displays the hydrogen bonds occuring in the molecule.
    :param molecule: PDB-ID of molecule HBsearch run should be performed on or in case an object in the PyMol session has to be used: The name of the object.
    :param molecule_name: If object fetched by PDB-ID. Standard set to blank.
    :param directory: Directory of the program folder. Directes working directory in PyMol to entered directory. Obsolete when script is started from program folder and HBsearch is started without changing the directory. Standard set to current working directory (".").
    :param use_object: If an object in the PyMol session has to be used for the HBsearch run. 0: No; 1: Yes. Standard set to 0. If set to 0. Molecule will be fetched by PDB-ID.
    :param remove_object: If created PDB-file in pdb_files folder should be deleted after HBsearch run. 0: No; 1: Yes. Standard set to 0.
    :param water_connections: Define if hydrogen bridges with water molecules involved should be considered. "0" for NO, "1" for YES.
    :param hb_file: HB-file sued to define possible hydrogen bond interactions. Standard set to hb-define.txt file
    :param solvent_key: If hydrogen bond bridges with solvent should be considered. Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.
    :param pse_file: File containing the chemical nature of the atoms. Standard set to period-table.info. To create own one see standard file for structure.
    :param connections: If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, but could be possible due to rotation of the residues. If connections = 1: Special conncetions will be taken into account. Standard set to 0: No special connectoins will be taken into account.
    """
    
    #Allowed inputs for various variables.
    allowed_use_object = ["0", "1"]
    allowed_remove_object = ["0", "1"]
    allowed_water_connections = ["0", "1"]
    allowed_solvent_key = ["NONE", "HOMO","MEMB"]
    allowed_connections = ["0", "1"]
    
    #Checking if user inputs are allowed.
    if use_object not in allowed_use_object:
        raise InputError(f"{use_object} is not allowed as input for\"use_object\". Only \"0\" and \"1\" are possible inputs")
    if remove_object not in allowed_remove_object:
        raise InputError(f"{remove_object} is not allowed as input for \"remove_object\". Only \"0\" and \"1\" are possible inputs")    
    if water_connections not in allowed_water_connections:
        raise InputError(f"{water_connections} is not allowed as input for\"water_connections\". Only \"0\" and \"1\" are possible inputs")
    if solvent_key not in allowed_solvent_key:
        raise InputError(f"{solvent_key} is not allowed as input for \"solvent_key\". Only \"NONE\", \"HOMO\" and \"MEMB\" are possible inputs")
    if connections not in allowed_connections:
        raise InputError(f"{connections} is not allowed as input for \"connections\". Only \"0\" and \"1\" are possible inputs")
    
    changeDirectory(directory) #Change directory to program folder directory. Needed PyMol directory is not set to program folder for HBsearch run
    
    #Checks if user wants to use own object in PyMol session or wants to fetch a protein structure from the PDB
    if use_object == "0": #User wants to fetch a protein from the PDB using a PDB-ID
        fetchPDB(molecule, molecule_name) #Uses fetch command to fetch PDB_ID and allows user to name the fetched object in PyMol session.
    elif use_object == "1": #User wants to use own object.
        useObject(molecule) #Object with entered name is saved as pdb-file in the pdb_files folder and is used for the following HBsearch run
    
    hbs_output = startHBsearch(molecule, hb_file, solvent_key, pse_file, connections) #Starts an HBsearch run with given parameters. Output is stored
    hbs_dataframe = readInHBS(hbs_output) #Converts HBsearch run to dataframe
    acceptor, donor = prepareLists(hbs_dataframe, water_connections) #Prepares acceptor and donor lists for displaying in PyMol
    
    run_information = "HBsearch" #Which program is run. Needed for displayDistances() and showSticks()
    
    if molecule_name == "": #Checking if custom molecule name was entered: If no name was entered, PDB-ID or saved object name is used.
        displayDistances(acceptor, donor, molecule, run_information) #Displays distances of hydrogen bond acceptors with their respective donors as distance objects without label in PyMol and groups distance objects according to the strucutre object they are based on.
        showSticks(acceptor,donor, molecule, run_information) #Shows residues participating in hydrogen bonds as sticks
    else: #When molecule name was entered, PyMol strucutre object posses molecule name. So this is used for following PyMol selection based commands.
        displayDistances(acceptor, donor, molecule_name, run_information) #Displays distances of hydrogen bond acceptors with their respective donors as distance objects without label in PyMol and groups distance objects according to the strucutre object they are based on.
        showSticks(acceptor,donor, molecule_name, run_information) #Shows residues participating in hydrogen bonds as sticks
    
    if remove_object == "1": #Checks if parameter remove_object is set to 1. If yes: created PDB-file in pdb_files folder is deleted.
        removeObject(molecule)


# In[53]:


#Creation of command in PyMol.
cmd.extend("hbsearch", hbsearch) #When read in in PyMol the script creates a command in PyMol which can be started via the PyMol command line.


# # HB-Network - Initialization

# In[54]:


class saveBot:
    def __init_():
        directory_name_save:str = None
        dictionary_save: Dict = None
        molecule_name_save: str = None
        water_connections_save: str = None

save = saveBot()


# In[55]:


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
    os.mkdir(os.path.normpath(f"./HB_network/{directory_name}")) #Error will be handled by os.mkdir() intern error report.
    
    #Informs user, where cluster files will be saved.
    print(f"HBnetwork run will be saved in", os.path.normpath(f"HB_network/{directory_name}"))
    
    return directory_name


# In[56]:


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

    try:
        hbn_output = subprocess.run(hbnetwork_dir, cwd = cluster_dir, capture_output=True, shell=True, check = True, text = True).stdout #Runs HBnetwork. Cluster files are created. Summary is given as output
    except subprocess.CalledProcessError:
        raise ProcessError("hb-network could not be run. Try with PDB-ID \"4AKR\" to make sure program works properly. It is possible that hb-network does not work with your biomolecule of choice.")

    return hbn_output


# In[57]:


def cleanHBnetwork(directory_name: str):
    pass
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
    


# In[58]:


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
                                                #Returns -1 if the word "Cluster" isn't found --> No Clusters existant or program not working properly
    
    if starting_index == -1:
        raise ProcessError("Either the molecule does not posses hydrogen bond clusters or something went wrong")
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


# In[59]:


def initiateHBnetwork(molecule:str, molecule_name = "", directory:str = ".", water_connections: str = "0", 
         use_object: str = "0", remove_object = "0", hb_file: str = "hb-define.txt", 
         solvent_key:str = "NONE", pse_file:str ="period-table-info.txt", connections: str = "0") -> None:
    """
    Initiates the HBnetwork run. Creates index directory containing all atoms participating in hydrogen bonds with their
    respective cluster. Needed for further showNetwork runs.
    Index dictionary, molecule name, as well as Cluster directory is saved in object "save" of class saveBot.
    :param molecule: Strucutre used for HBsearch run. HBnetwork is run on the same molecule using the output of HBsearch.
    :param directory_name: Saving folder for the HBnetwork run in HB_network subfolder in the program folder. 
    :param hb_file: HB-file sued to define possible hydrogen bond interactions. Standard set to hb-define.txt file
    :param solvent_key: If hydrogen bond bridges with solvent should be considered. Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.
    :param pse_file: File containing the chemical nature of the atoms. Standard set to period-table.info. To create own one see standard file for structure.
    :param connections: If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, but could be possible due to rotation of the residues. If connections = 1: Special conncetions will be taken into account. Standard set to 0: No special connectoins will be taken into account.
    """
    
    #Allowed inputs for various variables.
    allowed_use_object = ["0", "1"]
    allowed_remove_object = ["0", "1"]
    allowed_water_connections = ["0", "1"]
    allowed_solvent_key = ["NONE", "HOMO","MEMB"]
    allowed_connections = ["0", "1"]
    
    #Checking if user inputs are allowed.
    if use_object not in allowed_use_object:
        raise InputError(f"{use_object} is not allowed as input for\"use_object\". Only \"0\" and \"1\" are possible inputs")
    if remove_object not in allowed_remove_object:
        raise InputError(f"{remove_object} is not allowed as input for \"remove_object\". Only \"0\" and \"1\" are possible inputs")    
    if water_connections not in allowed_water_connections:
        raise InputError(f"{water_connections} is not allowed as input for\"water_connections\". Only \"0\" and \"1\" are possible inputs")
    if solvent_key not in allowed_solvent_key:
        raise InputError(f"{solvent_key} is not allowed as input for \"solvent_key\". Only \"NONE\", \"HOMO\" and \"MEMB\" are possible inputs")
    if connections not in allowed_connections:
        raise InputError(f"{connections} is not allowed as input for \"connections\". Only \"0\" and \"1\" are possible inputs")

    save.water_connections_save = water_connections
    
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

    #Created cluster dictionary serving as an indexing dictionary later on. Dictionary is saved in object "save".
    #Most important step for initalization.
    hbn_cluster_dict = indexHbnetwork(hbn_output)
    
    save.dictionary_save = hbn_cluster_dict #Indexing dictionary is saved to object "save" for furhter use later on.
                                            #True ouput of function, since this function should be called only once for each HBnetwork.
                                            #Searches should be possible without reoccuring initation runs. Therefore saving output in object.
    
    #Checking which molecule name should be saved for later use in the showNetwork() function.
    if molecule_name == "": #Checking if custom molecule name was entered: If no name was entered, PDB-ID or saved object name is used.
        save.molecule_name_save = molecule #saving PDB-ID of fetched molecule as molecule name or name of saved object is use_object = 1.
    else: #When molecule name was entered, PyMol strucutre object posses molecule name. So this is used for following PyMol selection based commands.
        save.molecule_name_save = molecule_name #saving molecule_name as the object name used by showNetwork(). Applies, if fetched strucutre is named.
    
    
    
    
    #Checks if parameter remove_object is set to 1. If yes: created PDB-file in pdb_files folder is deleted.
    if remove_object == "1":
        removeObject(molecule)  


# In[60]:





# In[61]:


#Creates a pymol command starting initiateHBnetwork().
cmd.extend("initiateHBnetwork", initiateHBnetwork)


# # HB-Network PyMol-Display

# In[62]:


def readoutHBnetwork(query: str, checkFor: str = "RESIDUE" ) -> str:
    """
    Checks the index dictionary for hydrogen bond clusters, the query atom/residue participates and returns them in a list.
    :param query: The residue/atom of which the hydrogen cluster should be displayed in PyMol. Format: Chain/Residue ID/ for displaying clusters the whole residue participates. Chain/Residue ID/atom for displaying clusters the chosen atom participates.
    :param checkFor: Viable inputs are: ATOM and RESIDUE. Input determines if clusters where just the atom participates should be shown or the cluster where the whole residue participates.
    :return destination_cluster_list: List with all hydrogen bond cluster the residue/atom participates in.
    """
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
        
    
    #Checks if all clusters, where any atom of the residue participate, should be processed.
    elif checkFor =="RESIDUE":
        
        #Query statement for residues needs to be closes with "/". Otherwise multiple clusters starting with the entered query are chosen for displaying.
        if not query[-1] == "/":
            raise SelectionError("Please close your query statement for whole residues with \"/\" so the selection is unambiguous.")
        
        #Appends all cluster matching the query input of "chain/residue/"
        [destination_cluster_list.append(value) for key, value in dictionary.items() if query in key]
    
    if not destination_cluster_list:
        raise SelectionError("The selected atom/residue could not be found or does not participate in hydrogen bonds")
    
    #Returns cluster list with all clusters the residue/atom participates. 
    return destination_cluster_list #Cluster list is passed to next function for searching respective cluster files.
                                    #Returns blank list, if no cluster, the residue participates, were found.


# In[63]:


print(readoutHBnetwork("A/3/"))


# In[64]:


def readInHBN(cluster_file_output: str) -> pd.DataFrame:
    """
    Given cluster files output genereted by HBnetwork is converted to a dataframe for furhter processing in following functions.
    :param cluster_file_output: Output of cluster file that should be conferted to pandas dataframe.
    :return df_cluster: Returns converted cluster file output as pandas dataframe.
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


# In[65]:


def prepareDataFrameHBnetwork(destination_cluster_list: str) -> pd.DataFrame:
    """
    Iterates through all cluster files indicated in cluster list, containing all clusters a residue/atom participates,
    processes them via readInHBN() to dataframes and appends all output dataframes to one, which is used for displaying
    connections and showing participating residues as sticks later on.
    :param destination_cluster_list: List containing all clusters a residue/atom participates.
    :return hBond_cluster_dataframe: Returns dataframe with all acceptor/donor relationships of all clusters a residue/atom participates.
    """
    
    #Extracts directory, where HBnetwork run is saved from "save" object
    directory = save.directory_name_save
    
    #Directory, where cluster files for initiation run are saved
    cluster_directory = os.path.normpath(f"./HB_network/{directory}/CLUSTER") #By HBnetwork created CLUSTER directory is targeted
    
    #Creating empty DataFrame. All dataframes created by reading out the generated cluster files is appended to this dataframe
    hBond_cluster_dataframe = pd.DataFrame()
    
    #Iterates through all files indicated in the cluster list containing all clusters a residue/atom participates.
    for file in destination_cluster_list: #Iterating over given clusters.
        
        #Save path to cluster file for reading in.
        cluster_dir = os.path.normpath(f"{cluster_directory}/{file}.hb-ntw")
        
        #Reading in cluster file.
        with open(cluster_dir,"r") as cluster_file:
            cluster_output = cluster_file.read()
            
        cluster_output = cluster_output.strip() #Cluster file is stripped, since last blanc line is interfering with creation of dataframe.
        
        #Processing cluster file output to dataframe
        if cluster_output: #Checks if cluster contains any hydrogen bonds.
            new_cluster_df = readInHBN(cluster_output) #Processes cluster file output to dataframe
        else: #If cluster is empty -> Continueing with next cluster of list.
            continue
        
        #Appending all dataframes from output of cluster files to one dataframe containing all acceptor/donor relationships
        #of all clusters the residue/atom participates.
        hBond_cluster_dataframe = pd.concat([hBond_cluster_dataframe, new_cluster_df], ignore_index = True) 
    
    return hBond_cluster_dataframe


# In[66]:


def showNetwork(query: str, checkFor = "RESIDUE") -> None:
    """
    Displays the hydrogen bridges between accetors and donors of all clusters the query residue/atom participates in PyMol.
    Before started initializeHBnetwork needs to be performed with the object of choice.
    :param query: The residue/atom of which the hydrogen bond network should be displayed in PyMol. Format: Chain/Residue ID/ for displaying clusters the whole residue participates. Chain/Residue ID/atom for displaying clusters the chosen atom participates.
    :param checkFor: Determines if hydrogen networks of the whole residue or just the atom should be displayed. Viable inputs: ATOM or RESIDUE.
    """
    
    
    #Prepares the list of all clusters the residue/atom participates.
    destination_cluster_list = readoutHBnetwork(query, checkFor)
    #Processes the cluster file outputs of respective clusters to dataframes and appends them to one single dataframe.
    hbn_dataframe = prepareDataFrameHBnetwork(destination_cluster_list)
    
    #Prepares acceptor/donor lists using the dataframe.
    if hbn_dataframe.empty: #Check if dataframe is empty. In case it is empty --> Stop function. Show notification.
        print(f"This {checkFor} {query} does not participate in any hydrogen bonds")
        return None
    
    else:
        water_connections = save.water_connections_save
        acceptor, donor = prepareLists(hbn_dataframe, water_connections) #Preparing acceptor/donor lists of hydrogen network

    #Extracting the object name of the object initializeHBnetwork was performed on from saving object "save".
    molecule_selection = save.molecule_name_save
        
    #Preparing the query name for better differentiation in PyMol. 
    query_adj = query.replace("/", "\\") #Need to change direction, since PyMol otherwise recognizes the selection name as the first distance. Dont know why.
    run_information = f"HB_network_{query_adj}" #Displaying in distance and stick representation names that HBnetwork was used and which residue/atom was targeted.
    
    #Displays the hydrogen bonds between acceptors and donors of a hydrogen network of a chosen residue/atom.
    displayDistances(acceptor, donor, molecule_selection, run_information)
    showSticks(acceptor, donor, molecule_selection, run_information)


# In[68]:


showNetwork("A/3/")


# In[31]:


#Creates PyMol command for showNetwork().
cmd.extend("showNetwork", showNetwork)


# In[32]:


print(os.getcwd())


# In[ ]:


"""<!--       _
       .__(.)< (MEOW)
        \___)   
 ~~~~~~~~~~~~~~~~~~-->"""


# In[ ]:





# In[ ]:




