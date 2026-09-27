% Neville_Method_Analysis.m
% Adaptive Function Approximation Using Neville Interpolation
%
% This script applies Neville's interpolation method to estimate function
% values at selected points using progressively refined interpolation nodes.
% The analysis evaluates convergence behavior, approximate relative error,
% and true relative error for nonlinear functions.

clear; clc; close all;

%% Analysis settings

tolerance = 1e-9;
maxIterations = 30;

%% Verification case: ln(x)

testX = [2, 2.2, 2.3];
testY = log(testX);
testEvaluationPoint = 2.1;

[verificationTable, verificationValue] = Neville_Interpolation_Method( ...
    testX, testY, testEvaluationPoint);

fprintf('\n============================================================\n');
fprintf('Neville Interpolation Verification\n');
fprintf('============================================================\n');
fprintf('Function: ln(x)\n');
fprintf('Evaluation point: %.4f\n', testEvaluationPoint);
fprintf('Interpolated value: %.12g\n\n', verificationValue);
fprintf('Neville interpolation table:\n');
disp(verificationTable);

%% Case A: Oscillatory exponential-polynomial function

functionA = @(x) x.^3 .* exp(x) .* cos(10*x) + 1;
evaluationPointA = -1;
domainA = [-3, 6];

[estimateA, approximateErrorA, iterationsA] = Neville_Function_Evaluation( ...
    functionA, evaluationPointA, domainA, tolerance, maxIterations);

trueValueA = functionA(evaluationPointA);
trueRelativeErrorA = abs((estimateA - trueValueA) / max(abs(trueValueA), eps));

fprintf('\n============================================================\n');
fprintf('Case A: Oscillatory Nonlinear Function\n');
fprintf('============================================================\n');
fprintf('Function: f(x) = x^3 e^x cos(10x) + 1\n');
fprintf('Domain: [%.2f, %.2f]\n', domainA(1), domainA(2));
fprintf('Evaluation point: %.4f\n', evaluationPointA);
fprintf('Interpolated value: %.12g\n', estimateA);
fprintf('Approximate relative error: %.3e\n', approximateErrorA);
fprintf('True relative error: %.3e\n', trueRelativeErrorA);
fprintf('Iterations required: %d\n', iterationsA);

%% Case B: Trigonometric rational function

functionB = @(x) tan(cos(sqrt(5) + sin(x)) ./ (1 + x.^3));
evaluationPointB = 5;
% Keep the interpolation interval away from the pole at x = -1.
domainB = [4, 6];

[estimateB, approximateErrorB, iterationsB] = Neville_Function_Evaluation( ...
    functionB, evaluationPointB, domainB, tolerance, maxIterations);

trueValueB = functionB(evaluationPointB);
trueRelativeErrorB = abs((estimateB - trueValueB) / max(abs(trueValueB), eps));

fprintf('\n============================================================\n');
fprintf('Case B: Trigonometric Rational Function\n');
fprintf('============================================================\n');
fprintf('Function: f(x) = tan(cos(sqrt(5)+sin(x))/(1+x^3))\n');
fprintf('Domain: [%.2f, %.2f]\n', domainB(1), domainB(2));
fprintf('Evaluation point: %.4f\n', evaluationPointB);
fprintf('Interpolated value: %.12g\n', estimateB);
fprintf('Approximate relative error: %.3e\n', approximateErrorB);
fprintf('True relative error: %.3e\n', trueRelativeErrorB);
fprintf('Iterations required: %d\n', iterationsB);
