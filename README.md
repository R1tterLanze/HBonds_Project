# HBonds_Project
This program is a PyMol script calculating which atoms in a protein interact via hydrogen bonds and displays the results in PyMol.

##Installation
To use this program, you will need to have the latest version of Python (best would be Anaconda) and PyMol, as well as the Pandas library, installed on your computer.


To install the latest version of Python, follow these steps:

Download the Python installer from the official website and run it to begin the installation process.
Follow the prompts to complete the installation.


To install PyMol, follow these steps:

Download the PyMol installer from the official website and run it to begin the installation process.

Or:
Open the Anaconda command prompt and run the following command to install the open souce version of PyMOL:
conda install -c conda-forge pymol-open-source


To install Pandas, follow these steps:

Open a command prompt or terminal and run the following command to install Pandas using pip:
pip install pandas
or:
conda install pandas (if using Anaconda)

Once you have Python, PyMol, and Pandas installed, you can download the entire program folder from this repository and place it in a convenient location on your computer. The program folder contains the script file (hb-script.py) as well as the necessary structures and executables


##Usage
To use the program, open PyMol and load the script. Either by running the script via command prompt (run path/hb-script.py) or by navigating through File -> Run script and selecting the script in the program folder to run.

This will execute the script and load the three functions (hbsearch, initiateHBnetwork, searchNetwork) but will not automatically execute them. 

To run the hydrogen bond calculation and display the hydrogen bonds, you will need to call the hbsearch function in the PyMol command prompt.

hbsearch molecule[, molecule_name[, directory[, 
         	use_object[, remove_object[, water_connections[,
		hb_file[, solvent_key[, 
		pse_file[, connections]]]]]]]]]

ARGUMENTS:

molecule = 			PDB-ID of molecule HBsearch run should be performed on or
				in case usage of an object in the PyMol session: The name of the object. 
				(Herefore, use_object has to be set to 1)
				If no existing object is used. PyMOL will try to fetch and display the protein with the entered PDB-ID.
				Further calculations will be performed this fetched structure.

molecule_name = 		If structure was fetched via PDB-ID using HBsearch, molecule_name can be used to name the object in the fetching process.
				Standard: ""
				If argument is empty the PDB-ID of the fetched protein will be used as object name.

directory =			Directory of the program folder. Directs working directory in PyMol to entered directory. 
				Obsolete when script is started via directing File -> Run script -> .../hb-script.py from program folder and HBsearch is started without changing the directory after loading the script. 
				Standard: Set to current working directory (".").

use_object =		Determines if an object in the PyMol session has to be used for the HBsearch run. 0: No; 1: Yes. 
				Standard set to 0. If set to 0: Molecule will be fetched by PDB-ID.

remove_object = 		Determines if created PDB-file in pdb_files folder should be deleted after HBsearch run. 0: No; 1: Yes. 
				Standard set to 0.

water_connections = 	Define if hydrogen bridges with water molecules involved should be considered. "0" for NO, "1" for YES. 
				Standard: 0

hb_file =			HB-file used to define possible hydrogen bond interactions. 
				Standard set to hb-define.txt file

solvent_key = 		Determines if hydrogen bond bridges with solvent should be considered. 
				Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.

pse_file = 			File containing the chemical nature of the atoms. 
				Standard set to period-table-info.txt To create own one see standard file for structure.

connections = 		If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, 
				but could be possible due to rotation of the residues. 
				If connections = 1: Special conncetions will be taken into account. 
				Standard set to 0: No special connectoins will be taken into account.


To calculate and display hydrogen bond clusters either by participating residue or participating atom, HBnetwork is used.
HBnetwork constis of two commands:
"initiateHBnetwork" calculates hydrogen bonds between atoms and residues and orders connected hydrogen bonds into clusters.
After initialization of the hydrogen bond network using "initiate HBnetwork, the cluster/s the atom/residue of choice participates can be displayed in the structure via the command "showNetwork".

initiateHBnetwork molecule[, molecule_name[, directory[, 
         	use_object[, remove_object[, water_connections[,
		hb_file[, solvent_key[, 
		pse_file[, connections]]]]]]]]]

HBsearch output and cluster files listing each atom participating in respective cluster are saved in the program folder under "./HB_network/molecule_name_date_time/

ARGUMENTS:

molecule = 			PDB-ID of molecule HBsearch run should be performed on or
				in case usage of an object in the PyMol session: The name of the object. 
				(Herefore, use_object has to be set to 1)
				If no existing object is used. PyMOL will try to fetch and display the protein with the entered PDB-ID.
				Further calculations will be performed this fetched structure.

molecule_name = 		If structure was fetched via PDB-ID using initializeHBnetwork, molecule_name can be used to name the object in the fetching process.
				Standard: ""
				If argument is empty the PDB-ID of the fetched protein will be used as object name.

directory =			Directory of the program folder. Directs working directory in PyMol to entered directory. 
				Obsolete when script is started via directing File -> Run script -> .../hb-script.py from program folder and HBsearch is started without changing the directory after loading the script. 
				Standard: Set to current working directory (".").

use_object =		Determines if an object in the PyMol session has to be used for the initialization of HBnetwork. 0: No; 1: Yes. 
				Standard set to 0. If set to 0: Molecule will be fetched by PDB-ID.

remove_object = 		Determines if created PDB-file in pdb_files folder should be deleted after intitialization. 0: No; 1: Yes. 
				Standard set to 0.

water_connections = 	Define if hydrogen bridges with water molecules involved should be considered. "0" for NO, "1" for YES. 
				Standard: 0

hb_file =			HB-file used to define possible hydrogen bond interactions. 
				Standard set to hb-define.txt file

solvent_key = 		Determines if hydrogen bond bridges with solvent should be considered. 
				Standard NONE: No solvent H-Bonds. Further possible: HOMO: Homogenous solvent; MEMB: Membrane environment.

pse_file = 			File containing the chemical nature of the atoms. 
				Standard set to period-table-info.txt To create own one see standard file for structure.

connections = 		If special connections should be taken into account. Here: Hydrogen bonds that would not be recognized by parameters, 
				but could be possible due to rotation of the residues. 
				If connections = 1: Special conncetions will be taken into account. 
				Standard set to 0: No special connectoins will be taken into account.



showNetwork query[, checkFor]

ARGUMENTS:

query = 			The atom or residue of which the hydrogen bond network cluster should be displayed
				Form for atoms as query: Chain/Residue ID/Atom e.g. A/6/OG
				Form for residues: Chain/Residue ID/ e.g. A/38/ (last "/" important here to close query)

checkFor =			Determines if the clusters of the whole residue should be displayed or only for the atom
				Possible inputs: 	"RESIDUE" for displaying the clusters a residue participates.
							"ATOM" for displaying the cluster an atom participates.

Limitations
This program is intended for use with protein structures in the PDB format. It may not work with other types of structures or file formats.

In addition, the program relies on a set of predefined criteria for identifying hydrogen bonds, which may not be suitable for all protein structures. As such, the results produced by the program should be interpreted with caution.

License
This program is provided under the terms of the MIT License. See the LICENSE file for details.