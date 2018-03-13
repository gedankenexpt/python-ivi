import os
import sys
import numpy as np
# begin: fix due to changed pyqt in conda
# import matplotlib
# matplotlib.use('qt5agg')
# end: fix due to changed pyqt in conda
# import matplotlib.pyplot as plt
import platform 
from time import strftime, gmtime, sleep

print('OS name %s' %os.name)
print('platform system %s' %platform.system())
print('platform release %s' %platform.release())

### FILE NAMING & PATH SETTING VARIABLES ###
if (platform.system() == 'Windows' and platform.release() == '7'): #Lab PC
    rootDir='C:\\Users\\njain\\data'
elif (platform.system() == 'Windows' and platform.release() == '10'): #office PC
    rootDir='C:\\Users\\nitjai\\research\\data'
elif (platform.system() == 'Linux'): #home PC
    rootDir = "/home/nitin/research/data" 

# # str42day=strftime("%Y-%m-%d")
# str4thyr=strftime("%Y")
# str4tmth=strftime("%B")
# # print ('The year is = %s ' %str4thyr)
# # print ('This month (in short) = %s ' %str4tmth[:3])
# rootDataDir=os.path.join(rootDir,str4thyr+'-'+str4tmth[:3])
# print(rootDataDir)
# if (os.path.exists(rootDataDir)):
    # print ('Path to the monthly directory exists')
# else:
    # print('something may be wrong')
    
# strSetday=gv.dirdayInfo[0:10]
#print 'first 10 chars of data dir name ', strSetday

global facGHz, facMHz, fackHz
facGHz=10**-9
facMHz=10**-6
fackHz=10**-3

global facmV, facuV
facuV=10**6
facmV=10**3

global facmW, facuW
facuW=10**6
facmW=10**3

global facns, facus, facms
facns=10**9
facus=10**6
facms=10**3

global comExtList
comExtList = ['','.csv','.txt','.CSV','.TXT']
# pathExptData = os.path.join(dirDataRoot,dirdayInfo)
# pathThryData=os.path.join(dirProgRoot, 'Theory-InterferenceCurves') #theory curves for signal interference.

def constrainVal(acv,mnv,mxv):
    if mnv > mxv:
        print ('max lim is smaller than min lim. Exiting')
        sys.exit(1)
    if acv < mnv:
        acv = mnv
    if acv > mxv:
        acv = mxv
    return acv

global fsz
fsz=24

### Various port-initializations (MAY CHANGE ON RE-CONNECTIONS, CHECK EVERYTIME!!!)###
portSFGpm = '/dev/ttyUSB1'
addrSFGpm = 4 #Check the instrument rear panel and pg 52 of 1830-R_Manual_RevA for the address and baud rate settings
portVODL='/dev/ttyUSB0'
portSIGpm = '/dev/ThorLabs_PM100D'
portPhCTR='/dev/ttyS0'
baudrate=19200 #Photon Counter

### Other global parameters ###
global stblzTime, mmtopsVODL, nPowMsmnts, timeDelPowMet

stblzTime = 0 #Time to let the system stabilize after a delay is set 
nPowMsmnts = 6 #choose higher value to get a better average powers
timeDelPowMet = 2 #time delay between succsessive calls to read the powermeter (in secs).
mmtopsVODL = float(20/3)

###?? Hardcoding of OWON oscilloscope parameters, such as ip-address, port, etc.
global dtype, skip, get, portOWON, ipaddr, numFpDel
dtype = "CH4"
skip = 0 #int(args_dict['-s'])
get = 0 #int(args_dict['-g'])
portOWON = 3000 #int(args_dict['-p'])
ipaddr = '192.168.1.72' #args_dict['-i']
numFpDel = 10 #No. of files to be captured per delay point
################################################################################
