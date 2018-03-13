function [data, symbols] = make_data(M, symbol_rate, ...
    symbol_length, fpilot, Apilot, fs, rrc_filter_coeff, statusLPF)
% M: size of signal constellation (4 for QPSK)
% symbol_rate: rate in Hz of symbols
% symbol_length: number of symbols
% fpilot: pilot tone frequency in Hz
% Apilot: amplitude of pilot tone
% fs: sampling frequency of AWG
% rrc_filter_coeff: root raised cosine filter coefficients
% statusLPF : if the 100 MHz LPF had been enabled on the AWG, then we need to invert its effect

%% make QPSK I/Q baseband data
symbols = randi([0 M-1], symbol_length, 1);
symbols = pskmod(symbols, M, pi/M);

%% upsample
data = resample(symbols, fs/symbol_rate, 1, rrc_filter_coeff);
t = (0:1/fs:(length(data)-1)/fs)';

%% add (single sideband) pilot tone
% figure;
% t1 = 0:1/symbol_rate:(symbol_length-1)/symbol_rate;
% plot(t1,rms(real(data))*real(symbols)/rms(real(symbols)),'-*',t,real(data),'-r')
% xlabel('Time (s)')
% ylabel('Signal')
% legend('symbol','data','Location','NorthWest')

if (fpilot) 
    data = data + Apilot*exp(1j*2*pi*fpilot*t);
end

%% scale data to +-1
scale = 0.95 / max(max(real(data)), max(imag(data)));
data = data * scale;

if (strcmpi(statusLPF,['Y'])) 
    load ch1_digFiltCoeffs_170MHz_st=12_n=3_m=6.mat b a; %I channel
    yI=filter(b,a,real(data));
    clear a b
    load ch2_digFiltCoeffs_170MHz_st=8_n=7_m=2.mat b a; %Q channel
    yQ=filter(b,a,imag(data));
    data = yI + 1j*yQ;    
end

end