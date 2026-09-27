% Horner_Method_Validation.m
% Validation of Polynomial Evaluation and Derivative using Horner's Method
%
% This script verifies the correctness of a custom Horner implementation
% by comparing results against MATLAB's built-in polynomial functions.
%
% The analysis evaluates:
% - Polynomial values P(x)
% - Derivative values P'(x)
% across multiple test points.

clear; clc; close all;

%% Polynomial definition

% Coefficients (highest degree → constant term)
polynomialCoefficients = [12, -3, 13, 0, -8, 19, -9, 211];

% Evaluation points
evaluationPoints = [-13, -7, 0, 5, 23, 37];

%% Horner method evaluation

[hornerValues, hornerDerivatives] = ...
    Horner_Polynomial_Evaluation(polynomialCoefficients, evaluationPoints);

%% MATLAB reference evaluation

matlabValues = polyval(polynomialCoefficients, evaluationPoints);

derivativeCoefficients = polyder(polynomialCoefficients);
matlabDerivatives = polyval(derivativeCoefficients, evaluationPoints);

%% Display comparison table

fprintf('\n============================================================\n');
fprintf('Horner Method Validation Against MATLAB Built-in Functions\n');
fprintf('============================================================\n');

fprintf('%-10s %-20s %-20s %-20s %-20s\n', ...
    'x', 'P(x) Horner', 'P(x) MATLAB', 'P''(x) Horner', 'P''(x) MATLAB');

fprintf('%s\n', repmat('-', 1, 100));

for index = 1:length(evaluationPoints)

    fprintf('%-10.2f %-20.6e %-20.6e %-20.6e %-20.6e\n', ...
        evaluationPoints(index), ...
        hornerValues(index), matlabValues(index), ...
        hornerDerivatives(index), matlabDerivatives(index));

end

%% Accuracy verification

valueError = max(abs(hornerValues - matlabValues));
derivativeError = max(abs(hornerDerivatives - matlabDerivatives));

fprintf('\n------------------------------------------------------------\n');
fprintf('Maximum absolute error (values):      %.3e\n', valueError);
fprintf('Maximum absolute error (derivatives): %.3e\n', derivativeError);

valueScale = max(1, max(abs(matlabValues)));
derivativeScale = max(1, max(abs(matlabDerivatives)));
if valueError / valueScale < 1e-12 && derivativeError / derivativeScale < 1e-12
    fprintf('\nValidation successful: Horner implementation is correct.\n');
else
    fprintf('\nWarning: Discrepancy detected. Review implementation.\n');
end
