# M4
Updated M4 by Kellen in Spring 2022.


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

