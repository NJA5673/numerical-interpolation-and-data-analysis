% Differentiation_Error_Analysis.m
% Numerical Differentiation and Error Convergence Analysis
%
% This script evaluates derivatives of highly nonlinear functions using a
% centered finite-difference method. The step size is repeatedly reduced by
% a factor of 10 to study convergence behavior, approximate relative error,
% and true relative error compared with symbolic differentiation.
%
% This analysis demonstrates the tradeoff between truncation error and
% round-off error in numerical differentiation.

clear; clc; close all;

%% Analysis settings

tolerance = 1e-9;
maxIterations = 20;
initialStepSize = 1;

%% Case A: Special-function expression

functionA = @(x) exp(erf(x.^2)) .* sin(cos(besselj(1, x) ./ (x.^3)));
evaluationPointA = 1 / sqrt(2);

fprintf('\n============================================================\n');
fprintf('Case A: Special-function derivative analysis\n');
fprintf('============================================================\n');

AnalyzeDerivativeConvergence(functionA, evaluationPointA, tolerance, maxIterations, initialStepSize, ...
    'exp(erf(x^2))*sin(cos(J_1(x)/x^3))');

%% Case B: Trigonometric nonlinear expression

functionB = @(x) tan(cos(sqrt(5) + sin(x)) ./ (1 + sin(x.^3)));
evaluationPointB = (1 + sqrt(5)) / 3;

fprintf('\n============================================================\n');
fprintf('Case B: Trigonometric nonlinear derivative analysis\n');
fprintf('============================================================\n');

AnalyzeDerivativeConvergence(functionB, evaluationPointB, tolerance, maxIterations, initialStepSize, ...
    'tan(cos(sqrt(5)+sin(x))/(1+sin(x^3)))');

%% Case C: Nested exponential-power expression

functionC = @(x) x.^(x.^(x.^x));
evaluationPointsC = [1e-4, 1.0002, 1.0003];

for pointIndex = 1:numel(evaluationPointsC)

    evaluationPointC = evaluationPointsC(pointIndex);

    fprintf('\n============================================================\n');
    fprintf('Case C: Nested power derivative analysis at x0 = %.4f\n', evaluationPointC);
    fprintf('============================================================\n');

    AnalyzeDerivativeConvergence(functionC, evaluationPointC, tolerance, maxIterations, initialStepSize, ...
        'x^(x^(x^x))');

end

%% Local analysis function

function AnalyzeDerivativeConvergence(targetFunction, evaluationPoint, tolerance, maxIterations, initialStepSize, functionLabel)
% ANALYZEDERIVATIVECONVERGENCE
% Runs centered finite-difference differentiation and compares the numerical
% derivative against MATLAB's symbolic derivative.
%
% Inputs:
%   targetFunction   - function handle to evaluate numerically
%   evaluationPoint  - point where derivative is estimated
%   tolerance        - stopping tolerance for approximate relative error
%   maxIterations    - maximum number of step-size refinements
%   initialStepSize  - starting finite-difference step size
%   functionLabel    - text label for printed output

    % Run numerical differentiation using decreasing step size
    [~, convergenceTable] = Numerical_Differentiation( ...
        targetFunction, evaluationPoint, tolerance, maxIterations, initialStepSize);

    % Symbolic derivative used as reference
    syms xSymbolic
    symbolicExpression = ConvertToSymbolicExpression(functionLabel, xSymbolic);
    symbolicDerivative = diff(symbolicExpression, xSymbolic);
    derivativeFunction = matlabFunction(symbolicDerivative, 'Vars', xSymbolic);

    trueDerivative = derivativeFunction(evaluationPoint);

    % Add true relative error column
    stepSize = convergenceTable(:, 1);
    derivativeApproximation = convergenceTable(:, 2);
    approximateRelativeError = convergenceTable(:, 3);

    trueRelativeError = abs((derivativeApproximation - trueDerivative) ./ max(abs(trueDerivative), eps));

    % Display results
    fprintf('Function: %s\n', functionLabel);
    fprintf('Evaluation point: %.8f\n', evaluationPoint);
    fprintf('Symbolic reference derivative: %.10e\n\n', trueDerivative);

    fprintf('%-14s %-18s %-20s %-20s\n', ...
        'h', 'D(approx)', 'Approx. Rel. Error', 'True Rel. Error');

    for row = 1:numel(stepSize)
        fprintf('%-14.4e %-18.8e %-20.4e %-20.4e\n', ...
            stepSize(row), derivativeApproximation(row), ...
            approximateRelativeError(row), trueRelativeError(row));
    end

end

function symbolicExpression = ConvertToSymbolicExpression(functionLabel, x)
% CONVERTTOSYMBOLICEXPRESSION
% Converts selected function labels into symbolic expressions for exact
% derivative comparison.

    switch functionLabel

        case 'exp(erf(x^2))*sin(cos(J_1(x)/x^3))'
            symbolicExpression = exp(erf(x^2)) * sin(cos(besselj(1, x) / x^3));

        case 'tan(cos(sqrt(5)+sin(x))/(1+sin(x^3)))'
            symbolicExpression = tan(cos(sqrt(5) + sin(x)) / (1 + sin(x^3)));

        case 'x^(x^(x^x))'
            symbolicExpression = x^(x^(x^x));

        otherwise
            error('Unknown symbolic expression label.');

    end

end