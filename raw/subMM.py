import time,sys
import numpy as np
import subprocess
import os
from multiprocessing import Pool


os.system("rm -rf ./out/*.out ")
def subOneAnalyzeJob(filename):
  analyzeexe = "/home/hesam/softwares/Tinker-CPU/source/analyze.x"
  cmdstr="%s ./txyz/%s.txyz -key tinker.key E> ./out/%s.out"%(analyzeexe, filename, filename)
  subprocess.run(cmdstr, shell=True)
  return

def subParallelJobs(filelist):
  filenames = [line.split(".txyz")[0] for line in open(filelist).readlines()]
  p = Pool(20)
  p.map(subOneAnalyzeJob, filenames)
  p.close()
  return

def getEnergy(params, filelist, savetxt=False):
  keys = [line.split()[0] for line in open("p0.txt").readlines()]
  oFile = open("temp.prm", 'w')
  lines = open("template.prm").readlines()
  for line in lines:
    if "PRM_" in line:
      for i in range(len(keys)):	
        search = "PRM_"+keys[i]+"_"
        if search in line and (line[0:1] != "#"):
          params[i] = "%10.8f"%(float(params[i]))
          line=line.replace(search, str(params[i]))
    oFile.write(line)
  oFile.close()
  subprocess.run("mv temp.prm amoeba09_Ln2+.prm", shell=True)
  subParallelJobs(filelist)
  #Check whether all analyze jobs finished!
  files = [line.split(".txyz")[0] for line in open(filelist).readlines()]
  nFile = len(files)
  readFlag = 0
  while readFlag==0:
    cmdstr1 = "grep 'Total Potential Energy' out/*.out > result.p"
    subprocess.run(cmdstr1, shell=True)
    nLines1 = sum(1 for line in open("result.p"))
    if (nLines1==nFile):
      readFlag = 1
      break
    else:
      time.sleep(0.5)

  MM_inter = [] 
  if readFlag==1:
    for eachFile in files:
      filename = "out/%s.out"%eachFile
      for line1 in open(filename).readlines():
        if "Intermolecular Energy " in line1:
          inter = float(line1.split()[-2])
          MM_inter.append(inter)
    MM_inter = np.array(MM_inter) 
    if savetxt:
      np.savetxt("MM-energy.dat", np.transpose([MM_inter]), fmt="%15.4f")
  return MM_inter

if len(sys.argv) == 2:
  getEnergy(np.loadtxt("p0.txt", usecols=(-1,)), "filelist", savetxt=True)
