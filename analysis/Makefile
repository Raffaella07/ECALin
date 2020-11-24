CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags)

CC = g++
std = c++17
CFLAGS = -Wall -c -g -fpermissive 

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I.. -I$(ROOTSYS)/include  -I../interface/   -I../../CommonTools/ -I../../CommonTools/interface/ -I../../../XGBoost-FastForest/build/ -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src 
ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -lRooFit -lTMVA -lRooFitCore -lHtml -lMathCore -lGenVector -lMinuit -lTreePlayer  

#LDLIBSOPTIONS=-llibxgboost -l../../../xgboost/build/rabit/librabit -l$(CMSSW_BASE)/src/xgboost/build/dmlc-core/libdmlc

MassRatios: masses_ratio.cpp BParkTools.o DileptonGeneralClass.o DiLeptonFromBClass.o
	$(CC) -Wall -Wextra -g $(INCLUDES) -o MassRatios masses_ratio.cpp BParkTools.o DileptonGeneralClass.o  DiLeptonFromBClass.o $(ROOTFLAG) $(EXTRALIBS)

Skim: skimmer.cpp DileptonClass.o DiLeptonClass.o MeeClass.o 
	$(CC) -Wall -Wextra -g $(INCLUDES) -o Skim skimmer.cpp DileptonClass.o DiLeptonClass.o MeeClass.o $(ROOTFLAG) $(EXTRALIBS)

SkimZ: skimmer.cpp DileptonGeneralClass.o MeeClass.o 
	$(CC) -Wall -Wextra -g $(INCLUDES) -o SkimZ skimmer.cpp DileptonGeneralClass.o MeeClass.o $(ROOTFLAG) $(EXTRALIBS)

AngularStudies: AngularVarStudies.cpp DileptonClass.o DiLeptonClass.o MeeClass.o 
	$(CC) -Wall -Wextra -g $(INCLUDES) -o AngularStudies AngularVarStudies.cpp DileptonClass.o DiLeptonClass.o MeeClass.o $(ROOTFLAG) $(EXTRALIBS)



BParkTools.o: ../src/BParkTools.cc DileptonGeneralClass.o

	$(CC) -Wall $(INCLUDES) -c  ../src/BParkTools.cc DileptonGeneralClass.o $(ROOTFLAG) $(EXTRALIBS)

BPark_fitUtils.o: ../src/BPark_fitUtils.cc 

	$(CC) -Wall $(INCLUDES) -c  ../src/BPark_fitUtils.cc $(ROOTFLAG) $(EXTRALIBS)

DiLeptonFromBClass.o: ../src/DiLeptonFromBClass.cc
	$(CC) -Wall $(INCLUDES) -c ../src/DiLeptonFromBClass.cc $(ROOTFLAG) $(EXTRALIBS)

DiLeptonClass.o: ../src/DiLeptonClass.cc
	$(CC) -Wall $(INCLUDES) -c ../src/DiLeptonClass.cc $(ROOTFLAG) $(EXTRALIBS)

DileptonClass.o: ../src/DileptonClass.cc
	$(CC) -Wall $(INCLUDES) -c ../src/DileptonClass.cc $(ROOTFLAG) $(EXTRALIBS)
MeeClass.o: ../src/MeeClass.cc
	$(CC) -Wall $(INCLUDES) -c ../src/MeeClass.cc $(ROOTFLAG) $(EXTRALIBS)
clean:
	rm *.o BParkTools.o
