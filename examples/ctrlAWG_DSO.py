
## Need to be resolved with Tobias : ?? 
## Need to be changed before the actual experiment : !!
## Calling scope.measurement.initiate a 2nd time gives the updated voltage offset but not impedance

## Acc to the docs, acquisition.record_length variable returns the actual number of points the DSO acquires for a given channel. Also, acquisition.record_length must be => acquisition.number_of_points_minimum, or the minimum number of points specified by the user. However, this does not seemingly work and it also screws up the calculation of the sampling rate. See below

# Min num points using inbuilt funct = 5.000e+06
# Record length using inbuilt funct = 1.000e+06 [WRONG]
# Record length acc to downloaded waveform = 5.000e+06
# fs_dso using inbuilt functn: 2.000006 GHz [WRONG]
# fs_dso based on inverse of mean time difference in downloaded waveform = 10.000 GHz (also the value seen on the DSO itself)
# Time per record using inbuilt funct = 5.000e-04 sec
# Time window acc to downloaded waveform = 5.000e-04 sec

## In the next run of the program, the record length value gets corrected. 

import time
import numpy
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
import globVarsAll as gv
import matplotlib.pyplot as plt
import ivi.tektronix.tektronixAWG5000 as awg5000
import ivi.agilent.agilentDSOS204A as dsos204a

facXt = gv.facus
facYv = gv.facmV

fs_awg = 500e6
fs_dso = 10e9
number_of_points_minimum = 1e6
time_per_record_DSO = 50e-6

awg = awg5000.tektronixAWG5000('TCPIP::10.54.4.47::INSTR')
scope = dsos204a.agilentDSOS204A('TCPIP0::WINDOWS-6TBN82B.local::hislip0::INSTR')

# setup awg
chI = 0
chQ = 1
t = numpy.arange(0,1e-3,1/fs_awg)
data1 = numpy.sin(2*numpy.pi*5e6*t)
data2 = numpy.cos(2*numpy.pi*5e6*t)

marker = numpy.ones(len(t))
marker[:100000] = 0
marker[-100000:] = 0
data1 = numpy.column_stack((data1, marker)) # using marker1 only
# data2 = numpy.column_stack((data2, marker)) # using marker1 only

awg.arbitrary.sample_rate = fs_awg
# awg.outputs[ch2].enabled = False
print('create waveform for chI')
wfm_handle1 = awg.arbitrary.waveform.create(data1)
awg.outputs[chI].arbitrary.waveform.handle = wfm_handle1
awg.outputs[chI].arbitrary.offset = 0.0
awg.outputs[chI].arbitrary.gain = 0.5
awg.outputs[chI].enabled = True

print('create waveform for chQ')
wfm_handle2 = awg.arbitrary.waveform.create(data2)
awg.outputs[chQ].arbitrary.waveform.handle = wfm_handle2
awg.outputs[chQ].arbitrary.offset = 0.0
awg.outputs[chQ].arbitrary.gain = 0.5
awg.outputs[chQ].enabled = True

# setup oscilloscope
ch1 = 0 
ch2 = 1 
ch3 = 2 #for marker 
ch4 = 3

scope.acquisition.type = 'normal'
scope.acquisition.input_frequency_max = 500e6 #Why is this needed?? And does it actually work??
scope.acquisition.number_of_points_minimum = number_of_points_minimum
scope.acquisition.time_per_record = time_per_record_DSO
print('Sample rate: {0} GHz'.format(scope.acquisition.sample_rate*gv.facGHz))

#default settings for the channels
scope.channels[ch1].input_impedance = 50 #change to 50 ohm with IQmodulator
scope.channels[ch1].coupling = 'dc'
scope.channels[ch1].offset = 0
scope.channels[ch1].range = 1
scope.channels[ch1].enabled = True

scope.channels[ch2].input_impedance = 50 #change to 50 ohm with IQmodulator
scope.channels[ch2].coupling = 'dc'
scope.channels[ch2].offset = 0
scope.channels[ch2].range = 1
scope.channels[ch2].enabled = True

scope.channels[ch4].enabled = False

scope.channels[ch3].input_impedance = 50 #change to 50 ohm with IQmodulator
scope.channels[ch3].coupling = 'dc'
scope.channels[ch3].offset = 0
scope.channels[ch3].range = 3
scope.channels[ch3].enabled = True

chenbld=[0,0,0,0]
onchs=[]
chnum=0
while (chnum<3):
    chenbld[chnum]=int(scope.channels[chnum].enabled)
    onchs.append('ch'+str(chnum+1))
    chnum = chnum + 1

# start awg
awg.initiate_generation()

time.sleep(0.5)

# start scope acquisition
scope.measurement.initiate() #May need to be able to do it continuously
# wait until scope is ready
while scope.measurement.status == 'in_progress':
    time.sleep(0.1)
    print('#')
    # should add some timeout here??

input('\nAdjust the scope vertical and horizontal settings, if reqd, and then press enter to continue\n')

print('Set ch1 impedance and coupling = %.2f and %s' %(scope.channels[ch1].input_impedance, scope.channels[ch1].coupling))
print('Set ch1 offset and range = %.2f and %.2f' %(scope.channels[ch1].offset, scope.channels[ch1].range))

print('Set ch2 impedance and coupling = %.2f and %s' %(scope.channels[ch2].input_impedance, scope.channels[ch2].coupling))
print('Set ch2 offset and range = %.2f and %.2f' %(scope.channels[ch2].offset, scope.channels[ch2].range))

# start scope acquisition
scope.measurement.initiate() 
# wait until scope is ready
while scope.measurement.status == 'in_progress':
    time.sleep(0.1)

print('Measurement status = %s' %scope.measurement.status)
print('Sample rate: {0} GHz'.format(scope.acquisition.sample_rate*gv.facGHz))
print('Continuing acquisition of %.2e points in a time window of %.2e ms \n' %(scope.acquisition.record_length, gv.facms*scope.acquisition.time_per_record))

if (chenbld[ch1]):
    measured_data1 = scope.channels[ch1].measurement.fetch_waveform()
else:
    print('Channel disabled / no data found on DSO channel ' + str(ch1+1))

if (chenbld[ch2]):
    measured_data2 = scope.channels[ch2].measurement.fetch_waveform()
else:
    print('Channel disabled / no data found on DSO channel ' + str(ch2+1))

if (chenbld[ch3]):
    measured_marker = scope.channels[ch3].measurement.fetch_waveform()
else:
    print('Channel disabled / no data found on DSO channel ' + str(ch3+1))

# stop awg and clear handles
awg.abort_generation()
if (wfm_handle1):
    awg.arbitrary.waveform.clear(wfm_handle1)
if (wfm_handle2):
    awg.arbitrary.waveform.clear(wfm_handle2)

# plot measured results and do some basic calculations
if (chenbld[ch1]):
    plt.plot(facXt*measured_data1.t, facYv*measured_data1.y, color='tab:orange', label='ch1')
    md1=numpy.array(measured_data1.y)
    numNaNmd1=numpy.count_nonzero(numpy.isnan(measured_data1.y))
    perNaNd1 = 100*numNaNmd1/len(measured_data1.y)
    print('Percentage of NaNs in ch1 acquired data = %.3f ' %(perNaNd1))
if (chenbld[ch2]):
    plt.plot(facXt*measured_data2.t, facYv*measured_data2.y, color='tab:green', label='ch2')
    md2=numpy.array(measured_data2.y)
    numNaNmd2=numpy.count_nonzero(numpy.isnan(measured_data2.y))
    perNaNd2 = 100*numNaNmd2/len(measured_data2.y)
    print('Percentage of NaNs in ch2 acquired data = %.3f ' %(perNaNd2))
if (chenbld[ch3]):
    plt.plot(facXt*measured_marker.t, facYv*measured_marker.y, color='tab:blue', label='marker')
if (True in chenbld):
    plt.ylabel('voltage (mV)',fontsize=20)
    plt.xlabel('time (us)',fontsize=20)
    plt.legend(fontsize=16)
    plt.grid()
    plt.show()

# shdbe1=md1**2+md2**2
# print(measured_data1.y)
# print('Sum of square of the two arrays; mean = '+str(numpy.nanmean(shdbe1))+' and std = '+str(numpy.nanstd(shdbe1)))

# hack to prevent exception errors in the end
del awg
del scope
print('out')