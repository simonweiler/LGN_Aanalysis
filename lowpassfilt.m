function filtered = lowpassfilt(signal, order, fc, fs, type)
%LOWPASSFILT  Switchboard for ML lowpass bessel and butterworth filters
%   filtered = LOWPASSFILT(signal, order, fc, fs, type) filters the input
%   signal at a cutoff freqeuency (fc, in Hz) of sampling rate (fs) using either a Bessel (type = 'Bessel') or
%   Butterworth (type = 'Butter') filter of 'order' order (aka pole).

if nargin < 5
    type = 'Bessel';
end

if type == 'Butter'
    [b,a] = butter(order,fc/(fs/2), 'low');
elseif type == 'Bessel'
    [b,a] = besself(order,fc * 2 * pi(), 'low');
    [b,a] = bilinear(b,a,fs);                       % digital to analog conversion (besself is analog only)
end

filtered = filtfilt(b,a,signal);                    % bidirectional filtering w/o phase distortion

  