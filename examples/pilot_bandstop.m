function Hd = pilot_bandstop(Fs, fpilot, stopband, passband)
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
Apass1 = 1;          % First Stopband Attenuation (dB)
Astop  = 60;           % Passband Ripple (dB)
Apass2 = 1;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, ...
    Apass1, Astop, Apass2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

% [EOF]
