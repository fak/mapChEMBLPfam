
def master(release, user, pword, host, port):
  ## Set threshold for all calculations.
  import numpy as np
  threshold = -np.log10(50*10**(-6))
  import export
  ## Load the propDict.
  import pickle
  infile = open('data/propDict_%s.pkl'% release, 'r')
  propDict = pickle.load(infile)
  infile.close()
  selected = ['Pkinase','Pkinase_Tyr','p450','SNF','Trypsin', 'RVP']
  export.exportFPs(selected, propDict, threshold, release)




if __name__ == '__main__':
  import sys

  if len(sys.argv) < 5:  # the program name and the two arguments

    sys.exit("Must specify release, user, pword, host, port")


  release = sys.argv[1]
  user = sys.argv[2]
  pword = sys.argv[3]
  host = sys.argv[4]
  port = int(sys.argv[5])

  master(release, user, pword, host, port)

