import time
import matlab.engine
import numpy
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
import globVarsAll as gv
import matplotlib.pyplot as plt
import ivi.tektronix.tektronixAWG5000 as awg5000
import ivi.agilent.agilentDSOS204A as dsos204a

print ('Matlab is starting....')
# Create a MatlabEngine object for communicating with the MATLAB process
eng = matlab.engine.start_matlab() #for options etc, see https://se.mathworks.com/help/matlab/apiref/matlab.engine.start_matlab.html
time.sleep(1)

# For redirecting standard output and error to Python
import io
out = io.StringIO()
err = io.StringIO()

# misc/decision variables
facXt = gv.facus
facYv = gv.facmV
plotYorN = True

#*#############################################################################################*#
#* Modulation related constants *#
#*#############################################################################################*#
symbol_rate = 10e6;
symbol_length = 10000.0; # number of samples (decimal reqd because of Python->Matlab's passing idiosyncracy)
M = 4.0;# Size of the signal constellation (decimal reqd because of Python->Matlab's passing idiosyncracy)

fpilot = 100e6;   # pilot tone frequency
Apilot = 0.1;    # amplitude of pilot

fcarrier = 1.1e9; #or 1.2e9     # carrier frequency in Hz set on ADF4351
fpilot_offset = 100e6; #100e6 # chosen offset in Hz of the pilot with respect to carrier

#*#############################################################################################*#
#* AWG and DSO paramters and initializations  *#
#*#############################################################################################*#
nbits_awg = 13.0# number of bits of AWG
fs_awg = 500e6  # sampling frequency of the AWG
fs_dso = 10e9   # sampling frequency of DSO

duration_chunk_awg = 1000e-6; # duration in seconds of one chunk

npts_min_DSO = 1e6
time_per_record_DSO = 50e-6

awg = awg5000.tektronixAWG5000('TCPIP::10.54.4.47::INSTR')
scope = dsos204a.agilentDSOS204A('TCPIP0::WINDOWS-6TBN82B.local::hislip0::INSTR')

# initialize awg channels, marker, sampling rate
chI = 0
chQ = 1
t = numpy.arange(0,duration_chunk_awg,1/fs_awg)

marker = numpy.ones(len(t))
marker[:100000] = 0
marker[-100000:] = 0

awg.arbitrary.sample_rate = fs_awg

# setup DSO channels
ch1 = 0 
ch2 = 1 
ch3 = 2 #for marker 
ch4 = 3

scope.acquisition.type = 'normal'
scope.acquisition.input_frequency_max = 500e6 #Why is this needed?? And does it actually work??
scope.acquisition.number_of_points_minimum = npts_min_DSO
scope.acquisition.time_per_record = time_per_record_DSO
print('Sampling rate: {0} GHz'.format(scope.acquisition.sample_rate*gv.facGHz))

#default settings for the channels
scope.channels[ch1].input_impedance = 1e6 #change to 50 ohm with IQmodulator
scope.channels[ch1].coupling = 'dc'
scope.channels[ch1].offset = 0
scope.channels[ch1].range = 1
scope.channels[ch1].enabled = True

scope.channels[ch2].input_impedance = 1e6 #change to 50 ohm with IQmodulator
scope.channels[ch2].coupling = 'dc'
scope.channels[ch2].offset = 0
scope.channels[ch2].range = 1
scope.channels[ch2].enabled = True

scope.channels[ch3].input_impedance = 1e6 #change to 50 ohm with IQmodulator
scope.channels[ch3].coupling = 'dc'
scope.channels[ch3].offset = 0
scope.channels[ch3].range = 4
scope.channels[ch3].enabled = True

scope.channels[ch4].enabled = False

chenbld=[0,0,0,0]
chnum=0
while (chnum<3):
    chenbld[chnum]=int(scope.channels[chnum].enabled)
    chnum = chnum + 1

#RRC filters for modualtion and demodulation
beta=0.35;
span=20.0;
rrc_filter_coeff_awg = eng.rcosdesign(beta, span, fs_awg/symbol_rate, 'sqrt');
rrc_filter_coeff_dso = eng.rcosdesign(beta, span, fs_dso/symbol_rate, 'sqrt');

tExecArr = []
evmArr = []
jj = 1
jjmax = 1
tf = -1
while True:
    try:
        # A = matlab.int8([1,2,3,4,5])
        # print(A[0][1:4])
        # # A = matlab.uint64([37])
        # # print('checking for primes\n')
        # # tf = eng.isprime(A)
        # # print(tf)
        # rrc_filter_coeff = eng.rcosdesign(beta, span, sps, 'sqrt'); #parr[0], parr[1], parr[2]
        # print('\nrrc_filter_coeff : ')
        # print(numpy.round(rrc_filter_coeff,3))
        # szrrcf = eng.size(rrc_filter_coeff)
        # print('\nDimensions returned by matlab : '+str(szrrcf))
        # print('\nDimensions returned by python : '+str(len(rrc_filter_coeff))+' x '+str(len(rrc_filter_coeff[0])))

        print('Creating IQ data + pilot to supply to the AWG\n')
        stm = time.clock() #deprecated, change to perf_counter() or process_time()
        #Having nargout as the last argument probably ensures the returned output is a tuple
        (awg_data, IQsamples) = eng.mod4py(M, symbol_rate, symbol_length, fpilot, Apilot, fs_awg, nbits_awg, rrc_filter_coeff_awg, duration_chunk_awg, nargout=2)
        data1 = numpy.real(awg_data)
        data2 = numpy.imag(awg_data)

        data1 = numpy.column_stack((data1, marker)) # using marker1 only

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

        # start awg
        awg.initiate_generation()

        time.sleep(0.5)

        # start scope acquisition to see waveforms on screeen
        scope.measurement.initiate() #May need to be able to do it continuously
        # wait until scope is ready
        while scope.measurement.status == 'in_progress': #isn't scope.measurement.status != 'complete' more reliable??
            time.sleep(0.1)
            print('#')
            # should add some timeout here??

        input('\nAdjust the scope vertical and horizontal settings, if reqd, and then press enter to continue\n')

        # print('Set ch1 impedance and coupling = %.2f and %s' %(scope.channels[ch1].input_impedance, scope.channels[ch1].coupling))
        # print('Set ch1 offset and range = %.2f and %.2f' %(scope.channels[ch1].offset, scope.channels[ch1].range))

        # print('Set ch2 impedance and coupling = %.2f and %s' %(scope.channels[ch2].input_impedance, scope.channels[ch2].coupling))
        # print('Set ch2 offset and range = %.2f and %.2f' %(scope.channels[ch2].offset, scope.channels[ch2].range))

        # start scope acquisition for acquiring waveforms
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
        if (plotYorN):
            plt.figure('DSO output for iteration no. '+str(jj))
            if (chenbld[ch1]):
                plt.plot(facXt*measured_data1.t, facYv*measured_data1.y, label='ch1')
                md1=numpy.array(measured_data1.y)
                numNaNmd1=numpy.count_nonzero(numpy.isnan(measured_data1.y))
                perNaNd1 = 100*numNaNmd1/len(measured_data1.y)
                print('Percentage of NaNs in ch1 acquired data = %.3f ' %(perNaNd1))
            if (chenbld[ch2]):
                plt.plot(facXt*measured_data2.t, facYv*measured_data2.y, label='ch2')
                md2=numpy.array(measured_data2.y)
                numNaNmd2=numpy.count_nonzero(numpy.isnan(measured_data2.y))
                perNaNd2 = 100*numNaNmd2/len(measured_data2.y)
                print('Percentage of NaNs in ch2 acquired data = %.3f ' %(perNaNd2))
            if (chenbld[ch3]):
                plt.plot(facXt*measured_marker.t, facYv*measured_marker.y, label='marker')
            if (True in chenbld):
                plt.ylabel('voltage (mV)',fontsize=20)
                plt.xlabel('time (us)',fontsize=20)
                plt.legend(fontsize=16)
                plt.grid()
                plt.show()

        # reconstructed = eng.demod4py(rawdata, fs_dso, fs_awg, symbol_rate, symbol_length, rrc_filter_coeff_awg, rrc_filter_coeff_dso, fcarrier, fpilot_offset, duration_chunk, IQsamples)
        # print(len(rawdata[0]))
        # print(IQsamples[0:4])

        ftm = time.clock() #deprecated, change to perf_counter() or process_time()

        # Analysis of rxdata and txdata
        # drop a few symbols in the beginning and the end
        # st = 10;
        # ed = symbol_length-10;

        # txdata = IQsamples[st:ed];
        # rxdata = reconstructed[0][st:ed];

        tExecArr.append(ftm-stm)
        # evmArr.append(tf)
        # print('\nIteration ' + str(jj) + ': Measured EVM = ' + str(numpy.round(tf,3)))
        print('\nIteration ' + str(jj) + ': Execution time = ' + str(numpy.round(ftm-stm,4)) + ' seconds\n')

        resp = input('Quit Matlab?\n')
        vb=str(resp)
        if (('y' in vb) or ('Y' in vb)):
            break
        else:
            print('Continuing...')

        # if (jj > jjmax):
            # break
        # else:
            # jj=jj+1


    except KeyboardInterrupt:
        print ('Keyboard interrupt encountered')
        resp = input('\nPausing... Would you like to: \n 1. Continue \n 2. Exit \n')
        if resp == '1':
            print('You have asked to resume')
            continue
        elif resp == '2':
            print('Exiting loop now')
            break

    except ValueError:
        print('Unknown value error encountered, quitting')
        break

    except IndexError:
        print('Unknown index error encountered, quitting')
        break

    except NameError:
        print('Unknown index error encountered, quitting')
        break

    except matlab.engine.MatlabExecutionError:
        print('Unknown engine error encountered, quitting')
        break


if (len(tExecArr) > 1):
    texm=numpy.mean(tExecArr[1:]) #drop first element for better stats 
    texs=numpy.std(tExecArr[1:]) #drop first element for better stats
    print('Mean and std dev execution time = %.3f and %.3f seconds\n' %(texm,texs))

# if (len(evmArr) > 1):
    # evmm=numpy.mean(evmArr[1:]) #drop first element for better stats 
    # evms=numpy.std(evmArr[1:]) #drop first element for better stats
    # print('Mean and std dev EVM value = %.3f and %.3f \n' %(evmm,evms))

# hack to prevent exception errors in the end
del awg
del scope

print('Out of loop and shutting Matlab down\n')
eng.quit()
