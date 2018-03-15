function Hd = pilot_bandpass(Fs, fpilot, passband, stopband)
%PILOT_BANDPASS Returns a discrete-time filter object.
%
% Fs: sampling frequency in Hz
% fpilot: pilot tone frequency in Hz
% passband: bandwidth of passband in Hz
% stopband: frequency difference in Hz from pilot and stopband of filter
%
% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

Fstop1 = fpilot-stopband;        % First Stopband Frequency
Fpass1 = fpilot-passband/2;        % First Passband Frequency
Fpass2 = fpilot+passband/2;        % Second Passband Frequency
Fstop2 = fpilot+stopband;        % Second Stopband Frequency
Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 60;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, ...
    Astop1, Apass, Astop2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
