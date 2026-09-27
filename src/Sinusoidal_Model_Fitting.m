function [meanValue, cosineCoefficients, sineCoefficients, fundamentalFrequency, fittedSignal, residualSumSquares] = ...
    Sinusoidal_Model_Fitting(timeData, signalData, highestHarmonic)
% SINUSOIDAL_MODEL_FITTING
% Least-squares sinusoidal model fitting for equally spaced data.
%
% This function approximates a periodic signal using a finite sinusoidal
% expansion up to a selected harmonic:
%
%   f(t) = A0 + Σ [Aj cos(j*w0*t) + Bj sin(j*w0*t)]
%
% The method computes the mean value, cosine coefficients, sine coefficients,
% fundamental angular frequency, fitted signal values, and residual error.
%
% Inputs:
%   timeData          - vector of equally spaced time points
%   signalData        - measured signal values
%   highestHarmonic   - highest harmonic included in the fit
%
% Outputs:
%   meanValue             - constant term A0
%   cosineCoefficients    - cosine coefficient array A
%   sineCoefficients      - sine coefficient array B
%   fundamentalFrequency  - fundamental angular frequency w0
%   fittedSignal          - fitted signal evaluated at input time points
%   residualSumSquares    - sum of squared residuals

%% Input validation

if ~(isscalar(highestHarmonic) && isnumeric(highestHarmonic) && ...
        isfinite(highestHarmonic) && highestHarmonic == floor(highestHarmonic) && ...
        highestHarmonic > 0)
    error('SinusoidalModelFitting:InvalidHarmonic', ...
        'highestHarmonic must be a positive integer.');
end

timeData = timeData(:);
signalData = signalData(:);

if numel(timeData) ~= numel(signalData)
    error('SinusoidalModelFitting:SizeMismatch', ...
        'timeData and signalData must have the same length.');
end

numSamples = numel(timeData);

if numSamples < 2
    error('SinusoidalModelFitting:TooFewSamples', ...
        'At least two data points are required.');
end

timeSteps = diff(timeData);
referenceStep = timeSteps(1);
if ~isfinite(referenceStep) || referenceStep <= 0 || ...
        any(~isfinite(timeData)) || any(~isfinite(signalData))
    error('Time must increase and all inputs must be finite.');
end

if any(abs(timeSteps - referenceStep) > 1e-10 * max(1, abs(referenceStep)))
    error('SinusoidalModelFitting:NonUniformSpacing', ...
        'timeData must be equally spaced.');
end

%% Fundamental frequency calculation

% The observed grid contains numSamples samples; its final point is not
% repeated as the first sample of the next period.
signalPeriod = numSamples * referenceStep;
fundamentalFrequency = 2*pi / signalPeriod;

%% Coefficient calculation

if 2*highestHarmonic + 1 > numSamples
    error('Too many harmonics for the number of samples.');
end
designMatrix = ones(numSamples, 2*highestHarmonic + 1);
for harmonicIndex = 1:highestHarmonic
    phase = harmonicIndex * fundamentalFrequency * (timeData - timeData(1));
    designMatrix(:, 2*harmonicIndex) = cos(phase);
    designMatrix(:, 2*harmonicIndex + 1) = sin(phase);
end
if rank(designMatrix) < size(designMatrix, 2)
    error('The requested harmonics are not identifiable from these samples.');
end
coefficients = designMatrix \ signalData;
meanValue = coefficients(1);
cosineCoefficients = coefficients(2:2:end).';
sineCoefficients = coefficients(3:2:end).';

%% Reconstruct fitted signal

fittedSignalColumn = designMatrix * coefficients;

fittedSignal = fittedSignalColumn.';

%% Residual analysis

residual = signalData - fittedSignalColumn;
residualSumSquares = sum(residual.^2);

end
