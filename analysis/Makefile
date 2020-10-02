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

MassRatios: masses_ratio.cpp BParkTools.o DiLeptonMassClass.o
	$(CC) -Wall -Wextra -g $(INCLUDES) -o MassRatios masses_ratio.cpp BParkTools.o DiLeptonMassClass.o  $(ROOTFLAG) $(EXTRALIBS)







BParkTools.o: ../src/BParkTools.cc DiLeptonMassClass.o

	$(CC) -Wall $(INCLUDES) -c  ../src/BParkTools.cc DiLeptonMassClass.o $(ROOTFLAG) $(EXTRALIBS)

BPark_fitUtils.o: ../src/BPark_fitUtils.cc 

	$(CC) -Wall $(INCLUDES) -c  ../src/BPark_fitUtils.cc $(ROOTFLAG) $(EXTRALIBS)


DiLeptonMassClass.o: ../src/DiLeptonMassClass.cc
	$(CC) -Wall $(INCLUDES) -c ../src/DiLeptonMassClass.cc $(ROOTFLAG) $(EXTRALIBS)


clean:
	rm *.o BParkTools.o
