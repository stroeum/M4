#!/bin/bash
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

LOCALDIR := /usr/lib/petsc

SRCDIR := src
BUILDDIR := build
TARGETDIR := bin
TARGET1 := $(TARGETDIR)/M4

SRCEXT := c

SOURCES1 := $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJECTS1 := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES1:.$(SRCEXT)=.o))

CFLAGS := -g -w -m64 -Os -Wwrite-strings -Wmissing-declarations -Wuninitialized
COMPILER = mpicc

LIB := \
 	-lpetsc -lmpi -lm -ldl -lc \
	-L$(PETSC_DIR)/lib

INCLUDE := \
	-Iinclude/ \
	-I$(PETSC_DIR)/include

all: $(TARGET1)
	@echo "main srcs: $(SOURCES1)"
	@echo "main objs: $(OBJECTS1)"

run:
ifeq ($(PETSC_DIR), $(LOCALDIR))
ifdef in
	mpiexec --oversubscribe -np 16 -H localhost $(TARGET1) -options_file ./input/$(in) | tee ./output/out.txt
else
	mpiexec --oversubscribe -np 16 -H localhost $(TARGET1) -options_file ./input/main.in | tee ./output/out.txt
endif
else
ifdef in
	mpiexec $(TARGET1) -mca btl_tcp_if_include eno1 -mca btl tcp,self -options_file ./input/$(in) > ./output/out.txt
else
	mpiexec $(TARGET1) -mca btl_tcp_if_include eno1 -mca btl tcp,self -options_file ./input/main.in > ./output/out.txt
endif
endif

save:
ifdef n
	mkdir -p SaveData/$(n);
	rm -f ./output/*.info
	rm -f ./SaveData/$(n)/*
	cp start.job ./SaveData/$(n)/
	cp ./output/* ./SaveData/$(n)/
	cp ./input/Profiles.in ./SaveData/$(n)/
	cp ./input/Partition.in ./SaveData/$(n)/
	cp ./input/main.in ./SaveData/$(n)/
else
	@echo "You must provide a directory name n=<...>"
endif

$(TARGET1): $(OBJECTS1)
	@mkdir -p $(TARGETDIR)
	@echo " Linking..."
	$(COMPILER) $(OBJECTS1) -o $(TARGET1) $(PETSC_LIB) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@echo $(PETSC_DIR)
	@mkdir -p $(BUILDDIR)
	$(COMPILER) $(CFLAGS) $(INCLUDE) -c -o $@ $<

clean_data:
	@echo " Cleaning...";
	$(RM) -r $(BUILDDIR) $(TARGETDIR) output/* M4*.err M4*.out
