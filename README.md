# Intelligent Design of High-Performance MOFs Based on Material Fingerprints and Transfer Learning

This repository contains the transfer learning model and the automated design submission script codes for hypothetical MOFs, along with the corresponding datasets described in our paper.

## Project Description
This project focuses on the automated design of high-performance MOFs (Metal-Organic Frameworks) through material fingerprints and transfer learning. The repository includes code for running transfer learning models, automating the MOF design process, and managing related datasets. This project aims to accelerate the discovery and optimization of hypothetical MOFs by predicting performance and facilitating efficient molecular simulation workflows.

## Project Organization
-----------------------
├── README.md                <- The top-level README for developers using this project.
├── main_auto.sh                       <- Main execution script
├── required_path.sh                       <- Stores paths and concurrency settings
├── __pycache__           <- Python cache files
├── all_nodes             <- All nodes data
├── all_topologies          <- All topologies data
├── best_edges_cif        <- Selected best substructure CIF files
├── best_edges_mol        <- Selected best substructure MOL files
├── best_edges_png        <- Selected best substructure PNG image files
├── gcmc_rest_mofs        <- Folder of MOFs that do not require RASPA calculation
├── gcmc_selected_mofs    <- Folder of MOFs selected for RASPA calculations
│   ├── DatatoCif.py             <- Script to convert LAMMPS output .DATA files to .CIF files
│   ├── error_list_energy_mini.txt             <- List of MOFs that failed during energy minimization
│   ├── lammps_input.in          <- LAMMPS script for energy minimization of MOFs
│   └── MOFs_gcmc.csv         <- List of MOFs for MC simulation
├── Linker_summary        <- Folder containing linker summaries
├── LINUX                <- Files for running in a Linux environment
│   ├── EQeq          <- Charged mol files required for RASPA calculations
│   ├── gcmc_required             <- Required files for running RASPA, including forcefield and molecular model files
│   ├── prepared_mofs                       <- MOF files ready for RASPA MC simulation
│   ├── MOFs_gcmc.csv         <- List of MOFs for MC simulation
│   └── simulation.input         <- Input script for running RASPA
├── new_designed_mofs     <- Folder of all newly designed MOFs
│   └── All_designed_mofs.xlsx           <- Data of all newly designed MOFs, including unit cell sizes, the one-hot encoding of their metal centers and topology structures, MACCS fingerprints  of organic linkers (255-bit encoding), and performance predictions from the pre-trained deep neural network model
├── New_MOF_summary       <- Folder for summaries of new MOFs
├── optimal_linker        <- Folder of optimal linkers
├── TL    <- Folder containing files related to transfer learning data
│   ├── __pycache__           <- Python cache files
│   ├── MOF_verify           <- Folder of MOFs used for verifying DNN model performance
│   │   ├── DatatoCif.py             <- Script to convert LAMMPS output .DATA files to .CIF files
│   │   ├── lammps_input.in          <- LAMMPS script for energy minimization of MOFs
│   │   └── MOFs_verify_model.csv         <- List of MOFs for MC simulation to verify DNN model
│   ├── calculate_R2.py           <- Script to calculate R2 of MOF performance predicted by fine-tuned pre-trained DNN model
│   ├── calculate_R2_original_DNN.py           <- Script to calculate R2 of MOF performance predicted by the original pre-trained DNN model
│   ├── counter.json           <- JSON file used to record the current loop count
│   ├── Data_output.py           <- Script for generating the Output.png
│   ├── Data_summaize.py           <- Script for summarizing the data into final_data.xlsx at the end of each loop
│   ├── final_data.xlsx           <- Summary of structure encodings (255-bit) and performance of all batches of MOFs after program completion
│   ├── Pretrained_model.ckpt           <- Pre-trained DNN model
│   ├── MOF_verify_model.xlsx           <- List of MOFs used for verifying the DNN model
│   ├── Output.png           <- Summary graph of the program"s final results, including the average performance, highest performance, and total number of MOFs in each batch
│   ├── TL_data_target_task_test.xlsx           <- List of MOFs used to validate the fine-tuned DNN model
│   ├── TL_data_target_task_train.xlsx           <- List of MOFs used to fine-tune the pre-trained DNN model
│   ├── Transferlearning_finetune.py           <- Script for fine-tuning the original pre-trained DNN model with some new MOFs, followed by transfer learning on the remaining MOFs
│   ├── Transferlearning_originaldnn.py           <- Script for performing transfer learning using the original pre-trained DNN model to predict the performance of the remaining MOFs
│   ├── TSN_cal.py           <- Script for calculating MOF performance based on MC simulation data
│   └── TSN_cal_verify.py           <- Script for calculating MOF performance based on MC simulation data, used to verify the DNN model
├── tobacco_1.0-master    <- Folder containing Tobacco_1.0 files used for MOF synthesis
├── 500000.sdf                  <- Large substructure database
├── Adsorption_simulation.sh                       <- GCMC simulation for gas adsorption
├── CiftoMol.py                 <- Script for converting CIF files to MOL files
├── colourMol.py                <- Script for coloring PNG image files
├── DefineXAtoms.py             <- Script for defining X atoms
├── EncodeMOFs.py               <- Script for encoding MOFs
├── Energy_minimization.sh                       <- Geometry optimization for MOF structures
├── Energy_minimization_verify.sh                       <- Optimization for verification MOFs only
├── EQeq_calculation.sh                       <- EQeq-based partial charge assignment
├── error_list_energy_mini.txt           <- List of MOFs that failed during energy minimization
├── error_list_eqeq.txt         <- List of MOFs that failed during charge calculation (EQeq)
├── error_list_obtain_data.txt         <- List of MOFs that failed to generate DATA file via lammps-interface
├── Linkers_summary.py          <- Script for summarizing the selected linkers in each cycle
├── MOF_design_part1.py         <- Part 1 of MOF design code
├── MOF_design_part2.py         <- Part 2 of MOF design code
├── MOFs_summary.py             <- Script for summarizing synthesized MOFs in each cycle
├── mol_with_atom_index.py      <- Script for processing MOL files with atom index
├── MolFormatConversion.py      <- Script for converting MOF file formats
├── MoltoCif.py                 <- Script for converting MOL files to CIF files
├── netcode.py                  <- Script for encoding nodes and topology structures
├── OptimalFingerprint.py       <- Script for selecting the optimal molecular fingerprint
└── ScreenLinkers.py            <- Script for screening linkers

## Requirements
1. A Python environment with version 3.7 or higher is required.
2. This script should be run in a Linux environment.
3. Make sure your Python environment contains the following packages: "RDKit", "torch", "scikit-learn", "openpyxl", "pytorch-ignite", "ase" and "openmpi", "pymatgen", "paramiko", "lammps-interface".
4. Make sure you have installed LAMMPS. [https://www.lammps.org/#gsc.tab=0]
5. Ensure that the MOLECULE and EXTRA-MOLECULE packages are enabled in LAMMPS.
6. Make sure you have installed RASPA. [https://iraspa.org/raspa/]
7. Make sure that the required molecular files (.def) and forcefield files are placed in "/LINUX/gcmc_required". Moreover, Raspa input files (.input) is placed in "/LINUX".
8. Make sure you have installed parallel. [http://ftpmirror.gnu.org/parallel/]

## Usage
1. Fill in the LAMMPS path in the "LAMMPS_DIR" variable of the "required_path.sh" file, ensuring that the "lmp_mpi" can be found in the "${LAMMPS_DIR}/src" directory.
2. Fill in the RASPA path in the "RASPA_DIR" variable of the "required_path.sh" file, ensuring that the "${RASPA_DIR}/bin", "${RASPA_DIR}/lib", and "${RASPA_DIR}/share" directories exist and contain the required files.
3. Fill in the number of concurrent tasks in the "required_path.sh" file, please ensure the server has enough CPU cores to support the maximum number of concurrent tasks.
4. Place the MOL files of the top-performing MOFs' linkers into the "/optimal_linker" directory.
5. Place the CIF files for all the nodes (.cif) and the TEMPLATE files for all the topologies (.template) into the "/all_nodes" directory and "/all_topologies" directory, respectively.
6. Also, copy the node files (.cif) and topology files (.template) into the "/tobacco_1.0-master/nodes_bb" directory and "/tobacco_1.0-master/templates" directory, respectively.
7. A batch of MOF data with computed performance needs to be added to "TL/final_data.xlsx". This data should include the MACCS fingerprints derived from the linkers, one-hot encodings for the nodes and topologies (the number of bits in the one-hot encoding should match the number of files in "/all_nodes" and "/all_topologies" directories), and the unit cell sizes (Length A, B, C) of the MOFs. The data should be organized in the following order: MOF name, MACCS fingerprint, node one-hot encoding, topology one-hot encoding, crystal size, and MOFs" performance. This worksheet should be named "Cycle 0".
8. Place a pre-trained model file named "Pretrained_model.ckpt" in the "/TL" directory. Ensure that the input layer size of "Pretrained_model.ckpt" matches the total number of MACCS fingerprints, as well as the number of files in the "/all_nodes" and "/all_topologies" directories.
9. Run the following commands in the Linux terminal:

   # 1. Activate the conda environment
   # Replace ${YOUR_ENV_NAME} with the name of your actual conda environment.
   # You can find the available environments by running: conda env list
   conda activate ${YOUR_ENV_NAME}
   
   # Example:
   # conda activate dl-env

   # 2. Grant the script execution permission
   # This command grants the "main_auto.sh" script the necessary execute permissions.
   chmod +x main_auto.sh
   
   # Example:
   # chmod 777 main_auto.sh

   # 3. Run the main script in the background using nohup
   # "nohup" allows the script to continue running even after you close the terminal.
   # The output will be saved in "output.log".
   nohup sh main_auto.sh > output.log 2>&1 &

## Example
After running the following command:

    conda activate dl-env
    chmod 777 main_auto.sh
    nohup sh main_auto.sh > output.log 2>&1 &

The expected output includes a set of MOF performance results stored in "TL/final_data.xlsx". Additionally, the program will generate a graph ("TL/Output.png") summarizing key performance metrics such as average performance, highest performance, and the number of MOFs for each batch. All generated MOFs will be saved in the "New_MOF_summary" folder, while the corresponding linkers will be stored in the "Linker_summary" folder.

## References
[LAMMPS for energy minimization] A. P. Thompson, H. M. Aktulga, R. Berger, D. S. Bolintineanu, W. M. Brown, P. S. Crozier, P. J. in "t Veld, A. Kohlmeyer, S. G. Moore, T. D. Nguyen, R. Shan, M. J. Stevens, J. Tranchida, C. Trott, S. J. Plimpton, LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales, Comp. Phys. Comm., 271, 10817 (2022). [https://doi.org/10.1016/j.cpc.2021.108171]

[LAMMPS-interface for generating LAMMPS files] P. G. Boyd, S. M. Moosavi, M. Witman & B. Smit, Force-Field Prediction of Materials Properties in Metal-Organic Frameworks. J. Phys. Chem. Lett. 8, 357-363 (2017). [https://dx.doi.org/10.1021/acs.jpclett.6b02532]

[RASPA for calculating MOFs' adsorption performance] D. Dubbeldam, S. Calero, D.E. Ellis, and R.Q. Snurr, RASPA: Molecular Simulation Software for Adsorption and Diffusion in Flexible Nanoporous Materials, Mol. Simulat., 42(2), 81-101 (2015). [http://dx.doi.org/10.1080/08927022.2015.1010082]

[ToBaCCo_1.0 for MOFs synthesis] Y. J. Colón, D. A. Gómez-Gualdrón, R. Q. Snurr, Topologically guided, automated construction of metal–organic frameworks and their evaluation for energy-related applications, Cryst. Growth Des. 17, 5801-5810 (2017). [https://doi.org/10.1021/acs.cgd.7b00848]

[EQeq for charge equilibration] C. E. Wilmer, K. C. Kim, R. Q. Snurr, An Extended Charge Equilibration Method, J. Phys. Chem. Lett., 3(17), 2506-2511 (2012). [http://doi.org/10.1021/jz3008485]

[RDKit] RDKit: Open-source cheminformatics. (https://www.rdkit.org). [https://doi.org/10.5281/zenodo.591637]

## License
This project is licensed under the [Creative Commons Zero v1.0 Universal License (CC0)] [https://choosealicense.com/licenses/cc0-1.0/]. See the LICENSE file for more details.

## Credits
This work is supported by Guangzhou Key Laboratory New Energy and Green Catalysis, Joint Institute of Guangzhou University & Institute of Corrosion Science and Technology, Guangzhou University and Center for Research Computing at the University of Notre Dame for computational resources. Special thanks to the project collaborators for their contributions and insights.

## Citation
If you use this code or dataset in your research, please cite the following publication:"# ADTL_v1" 
"# ADTL_v1" 
