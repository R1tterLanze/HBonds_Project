# HBonds_Project
PyMol Hydrogen Bonding Calculator
This program is a PyMol script that calculates which atoms in a protein interact via hydrogen bonds and displays the results in PyMol.

nstallation
To use this program, you will need to have the latest version of Python and PyMol, as well as the Pandas library, installed on your computer.

To install the latest version of Python, follow these steps:

Download the Python installer from the official website and run it to begin the installation process.
Follow the prompts to complete the installation.
To install PyMol, follow these steps:

Download the PyMol installer from the official website and run it to begin the installation process.
Follow the prompts to complete the installation.
To install Pandas, follow these steps:

Open a command prompt or terminal and run the following command to install Pandas using pip:
pip install pandas

Once you have Python, PyMol, and Pandas installed, you can download the entire program folder from this repository and place it in a convenient location on your computer. The program folder contains the script file (hydrogen_bonding.py) as well as the necessary structures and executables
Usage
To use the program, open PyMol and load the protein structure you want to analyze. Then, in the PyMol command prompt, navigate to the folder where you saved the program folder and run the following command:

Copy code
run hydrogen_bonding.py
This will execute the script and load the three functions (hbsearch, initiateHBnetwork, searchNetwork) but will not automatically execute them. To run the hydrogen bonding calculation, 
you will need to call the hbsearch function in the PyMol command prompt, like this:

Copy code
hbsearch(molecule, molecule_name, directory, use_object, remove_object, water_connections, solvent_key, hb_file, pse_file, connections)
This will run the hbsearch function and perform the hydrogen bonding calculation using the specified parameters. The molecule parameter should be the PDB ID or object name of the molecule you want to analyze,
and the molecule_name parameter allows you to specify a custom name for the PDB file that will be fetched or used for the calculation. 
The directory parameter allows you to specify the directory where the program files are located,
and the use_object and remove_object parameters determine whether an existing PyMol object should be used for the calculation and whether the PDB file should be deleted after the calculation is complete.

The water_connections parameter determines whether hydrogen bonds to water molecules should be displayed, and the solvent_key parameter allows you to specify the solvent conditions to be considered. 
The hb_file and pse_file parameters allow you to specify the files that define the possible hydrogen bond interactions and the chemical nature of the atoms, respectively. 
Finally, the connections parameter determines whether special connections should be taken into account in the calculation.

Limitations
This program is intended for use with protein structures in the PDB format. It may not work with other types of structures or file formats.

In addition, the program relies on a set of predefined criteria for identifying hydrogen bonds, which may not be suitable for all protein structures. As such, the results produced by the program should be interpreted with caution.

License
This program is provided under the terms of the MIT License. See the LICENSE file for details.