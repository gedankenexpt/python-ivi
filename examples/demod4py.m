function demodulated = demod4py(data, fs_dso, fs_awg, symbol_rate, symbol_length, ...
    rrc_filter_awg, rrc_filter_dso, carrier, pilot_offset, chunk_length, original_samples)
% data: 1d array of raw data
% fs_dso: sampling rate of DSO in Hz
% fs_awg: sampling rate of AWG in Hz
% symbol_rate: symbol rate in Hz
% symbol_length: number of symbols
% rrc_filter_dso: root raised cosine filter coefficients (from fs_dso to symb_rate)
% rrc_filter_awg: root raised cosine filter coefficients (from fs_awg to symb_rate)
% carrier: carrier frequency in Hz (for bandpass filter)
% pilot_offset: frequency offset in Hz of pilot with respect to carrier
% chunk_length: length of a chunk in seconds
% original_samples: I/Q complex samples

disp('Inside demodulate function now')

nancond=isnan(data);
if (sum(nancond) > 0)
    disp('replacing NaNs by zeros in the DSO data')
    data(nancond)=0;
end

if (chunk_length > symbol_length/symbol_rate) %sanity check
    warning('Limiting chunk duration to max possible value')
    chunk_length = symbol_length/symbol_rate; %choose the maximal chunk duration possible 
end

% list to vectors
if (size(data,2) ~= 1)
    disp('transposing data before demodulating')
    data = transpose(data);
end

if (size(original_samples,2) ~= 1)
    disp('transposing IQsamples before demodulating')
    original_samples = transpose(original_samples);
end

dbgplots=0;
if (dbgplots)
    fsz=20;
    nmzFac=rms(real(original_samples));
end

sgn=1;
data = data - mean(data);

%% Bin into chunks of chunk_length
N = round(chunk_length*fs_dso);
Nbins = floor(length(data)/N);
data = data(1:N*Nbins);

t = (0:1/fs_dso:(numel(data)-1)/fs_dso)';

data = reshape(data, N, Nbins);
t = reshape(t, N, Nbins);

disp('Extracting pilot and estimate frequency and phase')
%% Extract pilot and estimate frequency and phase
Hd = pilot_bandpass(fs_dso, carrier+pilot_offset, 2e6, 50e6);
Hd.arithmetic='double';
pilot = filtfilt(Hd.sosMatrix,Hd.ScaleValues,data);
phase = unwrap(angle(hilbert(pilot)));
pilot = reshape(pilot, N, Nbins);
phase = reshape(phase, N, Nbins);

% estimate pilot frequency
pilot_frequencies = zeros(1, Nbins);
for ii = (1:Nbins)
    p = polyfit(t(:,1), phase(:,ii), 1);
    pilot_frequencies(ii) = p(1)/(2*pi);
end

% transform pilot to baseband and determine phase
pilot_down = pilot .* exp(1j*2*pi*pilot_frequencies.*t(:,1));
Hd = pilot_bb_lowpass_iir(fs_dso);
Hd.arithmetic = 'double';
pilot_down = filtfilt(Hd.sosMatrix,Hd.ScaleValues,pilot_down);
pilot_phases = angle(pilot_down);
pilot_phases = unwrap(mean(pilot_phases(floor(N/10):end,:), 1)); % just drop some initial values
%pilot_phases = unwrap(mean(angle(pilot_down), 1)); % uncut version


disp('Down-converting quantum signal')
%% Down-convert quantum signal
data = data - mean(data, 1);
quantum = data .* exp(1j*2*pi*(pilot_frequencies-pilot_offset).*t(:,1));

%% compare at the AWG sampling rate to estimate the propagation delay
if (numel(rrc_filter_awg) > 1) %call demodulate() with rrc_filter_awg = 0 to skip delay compensation
    [us,ds]=rat(fs_awg/fs_dso);
    % sanity check
    if (us-floor(us)~=0) || (ds-floor(ds)~=0)
        error('Ratio of sampling rates of DSO and AWG cannot be easily rationalized')
    else
        dsfac=ds/us;
    end

    % upsample symbols to AWG rate 
    awg_data = resample(original_samples, fs_awg/symbol_rate, 1, rrc_filter_awg);
    % downsample DSO data to AWG rate 
    downsampled2i = resample(quantum, us, ds); 
    % remove pilot tone, either with bandstop or with LPF
    Hd = pilot_bandstop(fs_awg, abs(pilot_offset), 2e6, 50e6);
    Hd.arithmetic='double';
    downsampled2i = filtfilt(Hd.sosMatrix,Hd.ScaleValues,downsampled2i); %filtfilt does zero-phase filtering so no new delay should be created here!? 
    
    dataI = real(downsampled2i);
    dataQ = imag(downsampled2i);
    awgdataI = real(awg_data);
    awgdataQ = imag(awg_data);

    [acorI,lagI] = xcorr(awgdataI, rms(awgdataI)*dataI/rms(dataI));
    [~,I1] = max(abs(acorI));
    [acorQ,lagQ] = xcorr(awgdataQ, rms(awgdataQ)*dataQ/rms(dataQ));
    [~,Q1] = max(abs(acorQ));

    lagdI = lagI(I1);
    lagdQ = lagQ(Q1);
    sgnI=sign(acorI(I1));
    sgnQ=sign(acorQ(Q1));
    if (sgnQ ~= sgnI)
        warning('I and Q waveforms do not result in the same sign of the xcorrelation peak')
        sgn=-1;
    end
    if (lagdI ~= lagdQ)
        warning('I and Q waveforms do not result in the same lag for the xcorrelation peak')
    end
    
    lag = 0.5*(lagdQ+lagdI); %lagdQ

    sprintf('Delay to be compensated = %.3f', lag)
    
    if (dbgplots)
        t_awg = (0:1/fs_awg:(length(awg_data)-1)/fs_awg)';
        dataI=circshift(dataI,lagdI);
        dataQ=circshift(dataQ,lagdQ);
        figure;
        subplot(211)
        st=100;
        ed=ceil(length(awg_data)/50);
        plot(t_awg(st:ed),nmzFac*awgdataI(st:ed)/rms(awgdataI),'-sr',t_awg(st:ed),nmzFac*dataI(st:ed)/rms(dataI),'-xb','LineWidth',1.0)
        xlabel('Time (s)')
        ylabel('I waveform')
        set(gca,'FontSize', fsz);
        legend('txData to AWG','delay-compensated rxData', 'Location','NorthEast')
        grid on
        subplot(212)
        plot(t_awg(st:ed),nmzFac*awgdataQ(st:ed)/rms(awgdataQ),'-sr',t_awg(st:ed),nmzFac*dataQ(st:ed)/rms(dataQ),'-xb','LineWidth',1.0)
        xlabel('Time (s)')
        ylabel('Q waveform')
        set(gca,'FontSize', fsz);
        legend('txData to AWG','delay-compensated rxData', 'Location','NorthEast')
        grid on

        figure;
        rawdata = dataI + 1j*dataQ;
        subplot(211)
        plotPwelchEstimate(awg_data,fs_awg,'twosided');    
        legend('input-AWG')
        title('Comparison : AWG input with (resampled) DSO upconverted output')
        subplot(212)
        plotPwelchEstimate(rawdata,fs_awg,'twosided');
        legend('output-DSO-resampled')        
    end
    
    %% compensate delay of the dso data using above
    shiftedI = circshift(real(quantum),floor(dsfac*lag));
    shiftedQ = circshift(imag(quantum),floor(dsfac*lag));

    shifted = shiftedI + sgn*1j*shiftedQ;
else
    shifted = quantum;
end
%% Downsample to symbol rate with the rrc filter
downsampled2b = resample(shifted, 1, fs_dso/symbol_rate, rrc_filter_dso);

N_down = floor(chunk_length * symbol_rate);
sprintf('arguments info: length = %d,  N_down = %d, Nbins = %d ', length(original_samples), N_down, Nbins)
samples = reshape(original_samples, N_down, Nbins);

% 
% t = (0:1/fs_awg:(length(awgdataI)-1)/fs_awg)';
% figure;
% st=1000;
% ed=6000;
% subplot(211)
% plot(t(st:ed),rms(awgdataI)*dataI(st:ed)/rms(dataI),'-r',t(st:ed),awgdataI(st:ed),':k*')
% legend('rotated','original')
% xlabel('Time (s)')
% ylabel('Signal')
% subplot(212)
% plot(t(st:ed),rms(awgdataI)*shiftedI(st:ed)/rms(shiftedI),'-r',t(st:ed),awgdataI(st:ed),':k*')
% legend('rotated+shifted','original')

% t1 = 0:1/symbol_rate:(symbol_length-1)/symbol_rate;
% data = quantum;
% t = (0:1/fs_dso:(length(data)-1)/fs_dso)';
% figure;
% plot(t,real(data),'-r',t1,rms(real(data))*real(original_samples)/rms(real(original_samples)),'-*')
% xlabel('Time (s)')
% ylabel('Signal')
% legend('RawSymbols','DwnSampledData','Location','NorthWest')


%% Phase rotation
% Determine phase offset using data from 1st chunk
fun = @(theta)corr_rotated(downsampled2b(:,1), samples(:,1), theta, 'evm');
optimal_phase = fminbnd(fun, -pi, pi);

sprintf('Calculated optimal phase = %.3f', optimal_phase)

% Rotate phase
reconstructed = downsampled2b .* exp(1j*(optimal_phase-pilot_phases+pilot_phases(1)));

%% Reshape and normalize (I and Q separately)
reconstructed = reshape(reconstructed, symbol_length, 1);
reconstructed = reconstructed - mean(real(reconstructed));
reconstructed = reconstructed - mean(imag(reconstructed));
demodulated = real(reconstructed) / std(real(reconstructed)) + ...
    1j*imag(reconstructed) / std(imag(reconstructed));

%equalize in rms power
demodulated = demodulated * rms(samples)/rms(demodulated);

%% Analysis of rxdata and txdata
%drop a few symbols in the beginning and the end
st = 10;
ed = symbol_length-10;

txdata = samples(st:ed);
rxdata = demodulated(st:ed);

evm = comm.EVM; %Read here https://se.mathworks.com/help/comm/ref/comm.evm-system-object.html;jsessionid=577b23876abe0c0d97d5e8687fff#bsnan5l-2_1 for more
rmsEVM2 = step(evm,txdata,rxdata);
sprintf('Evm value = %.3f', rmsEVM2)

end
