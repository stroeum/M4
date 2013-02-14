DATE  = `date +%Y%m%d-%H%M`
FLAGS =-Wall -m64 -Os # -wd981 -g -Werror -wd593 -Wwrite-strings -Wmissing-declarations -Wuninitialized -Wstrict-prototypes -Wmissing-prototypes -std=c99 -Wshadow
CFLAGS = $(FLAGS)

FOLDER = $(CURDIR)
NAME = $(notdir $(CURDIR))
MYDIR=/nv/hp5/jriousset6/data/M4/${NAME}
NODES=1
PPN=64
PROCS= $(shell echo ${NODES}*${PPN} | bc)

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

#SOURCES = $(wildcard *.c)

SOURCES = MyProfiles.c MyCtx.c MyMonitor.c MyPDEs.c

OBJECTS = $(SOURCES:.c=.o)

all: $(OBJECTS) main_exe conv_exe

main_exe: main.o chkopts ${OBJECTS}
	${CLINKER} -o main main.o ${OBJECTS} ${PETSC_LIB}
	
conv_exe: convert.o chkopts ${OBJECTS}
	${CLINKER} -o convert convert.o ${OBJECTS} ${PETSC_LIB}

run:
	#cp ~/M4_RK.5 ~/${NAME}.sh
	cp ~/M4_RK.6 ~/${NAME}.sh
	qsub -N log${NAME} -l nodes=${NODES}:ppn=${PPN} -l walltime=100:00:00 -l pmem=3gb -v MYDIR=${MYDIR},MYPROCS=${PROCS} -z ~/${NAME}.sh 
	#mpirun -n 8 ./main -options_file input/main.in # > viz_dir/log.out
	
conv:
	./convert -options_file input/main.in
	
run_help:	
	reset
	./main -options_file input/main.in -help

clr_dat:
	rm -rfv .nfs* .*.swp input/.*.swp viz_dir/*.dat viz_dir/*.general viz_fun/*.*if* viz_fun/*.png viz_fun/*.pdf viz_fun/*.eps viz_fun/*.ps 
clr_img:
	rm -rfv .nfs* .*.swp input/.*.swp viz_fun/*.*if* viz_fun/*.png viz_fun/*.pdf viz_fun/*.eps viz_fun/*.ps 
clear:
	rm -rfv *.o .nfs* .*.swp main convert input/.*.swp output/* viz_dir/*.dat viz_dir/*.general viz_fun/*.*if* viz_fun/*.png viz_dir/log.out ~/${NAME}.sh ~/log${NAME}.*
