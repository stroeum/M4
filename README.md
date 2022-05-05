# M4
Last updated M4 by Kellen Sappington (ksappington2018@my.fit.edu) on 5/4/2022.


## To Install
You will need to run the following commands:
	sudo apt-get update
	sudo apt-get install make
	sudo apt-get install build-essential
	sudo apt-get install petsc-dev

For ssh to blueshark, you must type
	chmod 600 ~/.ssh/id_rsa (This file is your private key)


## Before running
Make sure that in your ~/.bashrc you are exporting the proper PETSC_DIR

	Locally:
		PETSC_DIR ~= /usr/lib/petsc  or  /usr/lib/petsc-<version>

	Blueshark:
		PETSC_DIR ~= export/spack/linux-centos7-x86_64/gcc-9.1.0/petsc-<version>


makefile must be tailored to the specific computer being used (primarily for directories, and hyperthreading locally).

For hyperthreading... (if you don't know about this, this may make the program run significantly faster on your computer)
the --oversubscribe should be enough, but if it isn't you'll likely need to create what's called a "rankfile".
For Windows users, you can check if hyperthreading is working by going to Task Manager -> Performance -> CPU
Here, "Logical processors" shows the number of threads/tasks you should be able to achieve via hyperthreading.
Right click on the display and "Change graph to -> Logical Processors". If when the program is running you see all at 100%
then you're all set.

## Commands to run

	1.  make clean_data
	2.  make all
	3.  make run
	4.  Optional: Can save data with
	    make n=<name> save
	    
	    Note: this will create a directory named SaveData which will hold all future saved data (This folder can be used for Matlab when analyzing data)

	Or, if on blueshark
	1. sbatch start.job


## Debugging information

M4/output directory must exist and must be manually created if it does not already exist!
(This should be changed in the future, such that the program automatically creates the directory if it does not already exist)

Seg faults:
	Most likely due to hardcoded directory pointing to the wrong place.
	The output, bin, and build directories must exist prior to running the mpiexec
	The input main.in file, myctx.c, and the makefile contain hardcoded directory pointers, double-check these.
	This also depends where the program is being run from. Should be /M4, not /M4/bin

btl_tcp error on BlueShark:
	This occurs when the program is attempting to run on >1 node. You need to add command-line arguments to mpiexec that look similar to
	-mca btl_tcp_if_include eno1 -mca btl tcp,self
	The eno1 can be located using the ifconfig command. Beyond this, I don't know much about the error which is why I hopefully say "look similar to"

Problem during Linking as part of "make all"
	May need to link a library using "sudo ln -s <library you do have> <library it's looking for and failing to find>
	Ex: sudo ln -s /usr/lib/python3 /usr/lib/python

	May also not have installed build-essential with "sudo apt-get install build-essential" (See "To Install" section near top)


## Introductory file system information

bin/: (stands for binary) This contains the binary executable, nothing to see, read, or edit here. This is created and destroyed by "make all" and "make clean_data" respectively.

build/: Only contains binary object files, same story as bin/.

include/: this is the folder containing what are called header files which establish all of the constants, functions, data structures, and file dependencies that the associated .c file needs

input/:   this holds all data/input that the simulation requires. There are 3 primary files, main.in, Profiles.in, and Partition.in.
	main.in: This is where the user sets all desired flags/parameters. This is like your "settings" for the simulation
	Profiles.in: This is effectively your first "slice" of data that the simulation will build off of overtime. Each line represents an altitude, there is 1 line for every 1 km from 400 km down to 100 km.
		within each line, the magnetic field geometry & magnitude, neutral gas density, electron density and temperature, ion and neutral temperature, and pressures are all defined at each km value.
	Parititon.in: This is formatted the same way as Profiles.in however this just includes the altitude, and the fraction of ion species. Each row should add to 1 to represent 100%

output/:  This is where all of the simulation output is temporarily stored. There are 3 main files, all .bin files, out.txt, and t.out.
	*.bin: These are the raw binary files containing all of the data. The number represents what "timestep" the file represents. X vs Y represents whether it was part of the 9 additional "ExtraDiagnostics" variables
	out.txt: This contains some basic run-time information such as the domain dimensions (mx, my, mz), some constants, your toggles (chemistry, collisions, etc.), and then documents information on each time the simulation
		saves data to *.bin files.
	t.out: This file just contains the time in seconds corresponding to every *.bin file
	***Note, any *.bin.info files are not necessary for anything and can be deleted

src/: (stands for source, as in source code) This is where all of the actual code is written (outside of the short header files in include/). Currently there are 7 files that make up the M4 program.
	MyChemistry.c: This just contains functions for all accounted for chemistry collisions and defines the reaction rates for each collision.
	MyCollisions.c: Same as above, except for collisions
	MyCtx.c: (Ctx = context) This contains many utility functions. The primary function is InitCtx which reads in all "settings" from input/main.in, and then sets up various parts of the simulation like the size of the domain
	MyMonitor.c: This deals with formatting and outputting the data calculated by the simulation. This creates all the output.
	MyPDEs.c: This is where 99% of the calculations and physics are done. You will also likely spend most your time here.
	MyProfiles.c: Just reads in the profiles and partition input files, and also defines some handy interpolation functions
	main.c: First file ran by the program on startup, this creates the PETSC and parallel-C environment.

viz_dir/: (stands for visualization_directory) This is where all of the figures and animations created by the viz_fun/M4_binary_viewer.mlx MATLAB script is stored.

viz_fun/: (visualization_functions) This is where all the suppporting MATLAB scripts are stored.

SaveData/: By using the "make n=<filename> save" commands, you can save whatever is currently in the output/ folder here, so that you can load it into the MATLAB scripts. Each folder within SaveData/ contains all of the binary files
	as well as the associated main.in, output.out, t.out, start.job, and input profile & partition. This helps with reproducing old simulation results

makefile: A makefile is a file that contains instructions for running and compiling a large program with multiple files. This sets up a lot of the "command-line arguments" that are necessary, and allows the user to simply use
	"make run" to begin the simulation. This also includes very handy functions like "make all", "make clean_data", and "make n=<filename> save".

start.job: This file is exclusively used by what's called the SLURM manager on a remote computer (in my case this was the BlueShark supercomputer). It contains information regarding how your program should be ran.
