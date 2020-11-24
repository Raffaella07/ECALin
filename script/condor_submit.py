import os
import string 
letter ='ABCD'


def makefile(part, let,dir):
	filename = "skim_part.%d_%s.condor"%(part,let)
	file = open(filename,"w")
	file.write("universe 	      = vanilla\n")
	file.write("executable            = /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/script/exec.sh\n")
	file.write("arguments             =  %s %d $(ProcId) 1 \n"%(let,part))
	file.write("output                = /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/script/output/SkimLin_%d%s$(ProcId).out\n"%(part,let))		
	file.write("error                 = /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/script/error/SkimLin_%d%s.$(ProcId).err\n"%(part,let))		
	file.write("log                   = /afs/cern.ch/work/r/ratramon/ECAL_linearity/ECALinAnalysis/ECALlinearity/script/log/SkimLin_%d%s.$(ProcId).log\n"%(part,let))	
	file.write("\n")
	file.write("request_memory        = 2000 \n")
	file.write("+MaxRuntime           = 244000\n")
	file.write("queue %d\n"%(dir))
	file.close()


for i in range(1,7):
	for j in letter:
		if j in ('A', 'B','C'):
			k = 2
			makefile(i,j,k)
			print("%d %s %d"%(i,j,k))
	#		if not os.path.isdir("/eos/home-r/ratramon/BToKEE_mini/%d_%s_%d_check/"%(i,j,k)):
	#		os.mkdir("/eos/home-r/ratramon/BToKEE_mini/%d_%s_%d_check/"%(i,j,k))
	
		elif j == 'D' :
			k =6
			makefile(i,j,k)
			print("%d %s %d"%(i,j,k))
	#		if not os.path.isdir("/eos/home-r/ratramon/BToKEE_mini/%d_%s_%d_check/"%(i,j,k)):
	# 			os.mkdir("/eos/home-r/ratramon/BToKEE_mini/%d_%s_%d_check/"%(i,j,k))
