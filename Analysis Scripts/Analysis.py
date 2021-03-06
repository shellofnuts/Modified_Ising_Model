import pandas as pd
import numpy as np
from scipy.stats import moment

import sys

opts = np.array([opt for opt in sys.argv[1:] if opt.startswith("-")])
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]


fileroute = "" # Top level file route
outputfile = ""
inputfiles = 10 # Number of files to average over
outputfile = ""
latticesizes = []

anisVal = 0

if "-f" in opts:
	fileroute = args[np.where(opts=="-f")[0][0]]
if "-n" in opts:
    inputfiles = int(args[np.where(opts=="-n")[0][0]])
if "-o" in opts:
	outputfile = args[np.where(opts=="-o")[0][0]]
if "-i" in opts:
	latticesizes = map(int, args[np.where(opts=="-i")[0][0]].strip('[]').split(','))
if "-a" in opts:
	anisVal = float(args[np.where(opts=="-a")[0][0]])
#else:
#    raise SystemExit(f"Usage: {sys.argv[0]} (-f | -n | -o) <arguments>...")

# Ideally will add a way to loop through all directories within the top level directory.

globalSystemAverage = pd.DataFrame()

systemAverages = pd.DataFrame()
systemSusc = pd.DataFrame()
systemBinder = pd.DataFrame()
systemSTD = pd.DataFrame()

for j in latticesizes:
	for i in range(inputfiles):
		inputTable = pd.read_csv(fileroute+'{}x{}_spinDist_{}_{}.csv'.format(j,j,anisVal,i+1), sep=',', header=None, index_col=0)
		inputTable = inputTable.fillna(value=np.nan)
		inputTable = inputTable.astype('float')
		systemAverages['System {}'.format(i+1)] = inputTable.mean(axis=1, skipna=True)
		systemSTD['System {}'.format(i+1)] = inputTable.std(axis=1, skipna=True)
		systemSusc['System {}'.format(i+1)] = inputTable.var(axis = 1, skipna=True) / inputTable.index
		systemBinder['System {}'.format(i+1)] = inputTable.apply(lambda x: moment(x, moment=2, nan_policy='omit'), axis=1) / inputTable.mean(axis=1, skipna=True)**2
	
	globalSystemAverage['{}x{} Average Magnetic Moment'.format(j,j)] = systemAverages.mean(axis=1, skipna=True)
	globalSystemAverage['{}x{} Average Susc'.format(j,j)] = systemSusc.mean(axis=1, skipna=True)
	globalSystemAverage['{}x{} Average U2'.format(j,j)] = systemBinder.mean(axis=1, skipna=True)
	globalSystemAverage['{}x{} Magnetic Moment std'.format(j,j)] = systemSTD.mean(axis=1, skipna=True)
	
	systemAverages = systemAverages.iloc[0:0]
	systemBinder = systemBinder.iloc[0:0]
	systemSusc = systemSusc.iloc[0:0]
	systemSTD = systemSTD.iloc[0:0]

globalSystemAverage.rename_axis(index="Temperature", inplace=True)
globalSystemAverage.to_csv('{}_{}.csv'.format(outputfile, anisVal), sep=',')