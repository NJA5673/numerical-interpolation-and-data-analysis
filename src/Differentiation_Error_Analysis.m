% Differentiation_Error_Analysis.m
% Finite-Difference Derivative Estimation and Error Convergence Study
%
% This script estimates derivatives of highly nonlinear functions using a
% centered second-order finite-difference approximation. The step size is
% reduced by a factor of 10 at each iteration to study convergence behavior.
%
% The numerical derivative is compared against symbolic differentiation to
% evaluate approximate relative error and true relative error.

clear; clc; close all;

%% Analysis settings

tolerance = 1e-9;
maxIterations = 20;
initialStepSize = 1;

%% Case A: Special-function derivative

functionA = @(x) exp(erf(x.^2)) .* sin(cos(besselj(1, x) ./ (x.^3)));
evaluationPointA = 1 / sqrt(2);

syms x
symbolicFunctionA = exp(erf(x^2)) * sin(cos(besselj(1, x) / x^3));

AnalyzeDerivativeCase( ...
    functionA, symbolicFunctionA, evaluationPointA, tolerance, ...
    maxIterations, initialStepSize, ...
    'Case A: Special-function expression');

%% Case B: Trigonometric nonlinear derivative

functionB = @(x) tan(cos(sqrt(5) + sin(x)) ./ (1 + sin(x.^3)));
evaluationPointB = (1 + sqrt(5)) / 3;

symbolicFunctionB = tan(cos(sqrt(5) + sin(x)) / (1 + sin(x^3)));

AnalyzeDerivativeCase( ...
    functionB, symbolicFunctionB, evaluationPointB, tolerance, ...
    maxIterations, initialStepSize, ...
    'Case B: Trigonometric nonlinear expression');

%% Case C: Nested exponential-power derivative

functionC = @(x) x.^(x.^(x.^x));
symbolicFunctionC = x^(x^(x^x));

evaluationPointsC = [1e-4, 1.0002, 1.0003];

for pointIndex = 1:numel(evaluationPointsC)

    AnalyzeDerivativeCase( ...
        functionC, symbolicFunctionC, evaluationPointsC(pointIndex), tolerance, ...
        maxIterations, initialStepSize, ...
        sprintf('Case C: Nested power expression at x0 = %.4f', evaluationPointsC(pointIndex)));

end

%% Local function

function AnalyzeDerivativeCase( ...
    targetFunction, symbolicFunction, evaluationPoint, tolerance, ...
    maxIterations, initialStepSize, caseLabel)
% ANALYZEDERIVATIVECASE
% Computes finite-difference derivative estimates and compares them against
% a symbolic derivative reference.

    fprintf('\n============================================================\n');
    fprintf('%s\n', caseLabel);
    fprintf('============================================================\n');

    % Numerical derivative convergence table
    [numIterations, convergenceTable] = Numerical_Differentiation( ...
        targetFunction, evaluationPoint, tolerance, maxIterations, initialStepSize);

    % Symbolic derivative reference
    syms x
    symbolicDerivative = diff(symbolicFunction, x);
    symbolicDerivativeFunction = matlabFunction(symbolicDerivative, 'Vars', x);
    trueDerivative = symbolicDerivativeFunction(evaluationPoint);

    % Extract convergence data
    stepSize = convergenceTable(:, 1);
    derivativeApproximation = convergenceTable(:, 2);
    approximateRelativeError = convergenceTable(:, 3);

    trueRelativeError = abs( ...
        (derivativeApproximation - trueDerivative) ./ max(abs(trueDerivative), eps));

    % Print summary
    fprintf('Evaluation point: %.10g\n', evaluationPoint);
    fprintf('Symbolic reference derivative: %.10e\n', trueDerivative);
    fprintf('Iterations used: %d\n\n', numIterations);

    fprintf('%-14s %-18s %-22s %-22s\n', ...
        'h', 'D(approx)', 'Approx. Rel. Error', 'True Rel. Error');

    fprintf('%s\n', repmat('-', 1, 82));

    for row = 1:numel(stepSize)
        fprintf('%-14.4e %-18.8e %-22.4e %-22.4e\n', ...
            stepSize(row), derivativeApproximation(row), ...
            approximateRelativeError(row), trueRelativeError(row));
    end

end