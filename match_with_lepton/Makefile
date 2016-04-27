
# .SUFFIX = .so

CPP             = g++
CPPFLAGS        = -g -O1 -Wall -fPIC -D_REENTRANT -Wno-deprecated -I. -L.

ROOTCFLAGS      := $(shell root-config --cflags)
ROOTLIBS        := $(shell root-config --libs) -lMinuit -lEG #-lg2c
# ROOTLIBS        += $(shell echo ${MY_LDLib})
ROOTGLIBS       := $(shell root-config --glibs)

CPPFLAGS        += $(ROOTCFLAGS)
LIBS            = $(ROOTLIBS) -lm
GLIBS           = $(ROOTGLIBS)

TARGET = main

#BASIC_OBJS =  jet.o base.o  main.o 
BASIC_OBJS = base_jet.o base_event.o base_hlt.o base_muons.o base_electrons.o main.o

all: ${TARGET} 
	echo "- - - - - - - - - - - - - - \n- - - - Successfully - - - - \n- - - - - - - - - - - - - -"

main : main.o ${BASIC_OBJS}
	$(CPP) -o $@ $(CPPFLAGS)  $(LIBS) $(BASIC_OBJS)


.cc.o:
	$(CPP) -c $(CPPFLAGS) $<
.C.o:
	$(CPP) -c $(CPPFLAGS) $<

clean: 
	-rm ${TARGET}
	-rm *~
	-rm *.o
	-rm *.pdf *.root 
