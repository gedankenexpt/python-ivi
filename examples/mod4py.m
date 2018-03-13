function [awg_data, IQsamples] = mod4py(M, symbol_rate, symbol_length, fpilot, Apilot, ...
    fs_awg, nbits_awg, rrc_filter_coeff_awg, duration_chunk)
%%paths to add
% addpath('../matlab-lib')
% addpath('../include') %for MAT files

%% decision variables
plotYorN=0; %0 if no plots are to be generated
procIndvSymbLvl=1;
enbld100MHzLPF = 'N'; %whether the 100 MHz LPF was enabled ('Y' for enabled, 'N' for disabled)
useExistingMATfiles = 'N';

%% file and folder info
% monDir='2018-Jan';
% dateDir='20180123 - upconverted IQdata';
% dateDir='2017-11-20-QPSK-electrical-upconversion';
% dateDir='20171214 - Elec IQmodulation'; %'2017-11-20-QPSK-electrical-upconversion' to be changed on daily basis
subdirname='Apil=0.1_beta=0.35_span=20';%Use '' if data is not inside any sub directory 
% subdirname='';

%% initialize the RNG
disp('Setting the RNG')
rng(1234568);

%% define root raised cosine filter
beta=0.35;
span=20;
if (rrc_filter_coeff_awg == 0.0)
    rrc_filter_coeff_awg = rcosdesign(beta, span, fs_awg/symbol_rate, 'sqrt');
end

%% load from file or make I/Q data
if (strcmpi(useExistingMATfiles,['Y'])) 
    IQsymbfname='IQsamples';% raw symbols
    IQdatafname='awg_data';%  data played by AWG (upsampled symbols + pilot)
    extn='.mat';
    if (~isempty(subdirname))
        IQsymbfname = fullfile(baseDir,monDir,dateDir,subdirname,strcat(IQsymbfname,extn));
        IQdatafname = fullfile(baseDir,monDir,dateDir,subdirname,strcat(IQdatafname,extn));
    else
        IQsymbfname = fullfile(baseDir,monDir,dateDir,strcat(IQsymbfname,extn));
        IQdatafname = fullfile(baseDir,monDir,dateDir,strcat(IQdatafname,extn));
    end
    s = load(IQdatafname,'-mat'); %save(IQsymbfname,'IQsamples')
    awg_data = s.awg_data;
    s = load(IQsymbfname,'-mat');
    IQsamples = s.IQsamples;
else
    [ad, iq] = make_data(M, symbol_rate, symbol_length, ...
        fpilot, Apilot, fs_awg, rrc_filter_coeff_awg, enbld100MHzLPF);
end

awg_data=ad;
IQsamples=iq;

end
