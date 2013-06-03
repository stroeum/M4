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
	qsub -N log${NAME} -l nodes=${NODES}:ppn=${PPN} -l walltime=40:00:00 -l pmem=3gb -v MYDIR=${MYDIR},MYPROCS=${PROCS} -z ~/${NAME}.sh 
	#mpirun -n 8 ./main -options_file input/main.in # > viz_dir/log.out
	
conv:
	./convert -options_file input/main.in
	
run_help:	
	reset
	./main -options_file input/main.in -help

clr_pgm:
	find .	     \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
	find .	     \( -name "main"   -o -name "convert" -o -name "*.o" -o -name "*.out" \) -exec rm -rfv {} \;

clr_dat:
	find .	     \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
	find .       \( -name "*.?if*" -o -name "*.png"     \) -exec rm -rfv {} \;
	find viz_dir \( -name "*.dat"  -o -name "*.general" \) -exec rm -rfv {} \;
	find output  \( -name "*.dat" 			    \) -exec rm -rfv {} \;

clr_img:
	find .	     \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
	find .       \( -name "*.?if*" -o -name "*.png"     \) -exec rm -rfv {} \;
	find .       \( -name "*.pdf"  -o -name "*ps" 	    \) -exec rm -rfv {} \;

clear:
	find .	     \( -name ".*.swp" -o -name ".nfs*"     \) -exec rm -rfv {} \;
	find .	     \( -name "main"   -o -name "convert" -o -name "*.o" -o -name "*.out" \) -exec rm -rfv {} \;
	find .       \( -name "*.?if*" -o -name "*.png"     \) -exec rm -rfv {} \;
	find viz_dir \( -name "*.dat"  -o -name "*.general" \) -exec rm -rfv {} \;
	find output  \( -name "*.bin*" -o -name "*.dat"     \) -exec rm -rfv {} \;
	rm -rfv  ~/${NAME}.sh ~/log${NAME}.*

stop:
	sh scripts/stop.sh
