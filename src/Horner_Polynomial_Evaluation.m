function [polynomialValue, derivativeValue] = Horner_Polynomial_Evaluation(coefficients, evaluationPoints)
% HORNER_POLYNOMIAL_EVALUATION
% Evaluates a polynomial and its derivative using Horner's method.
%
% This function efficiently evaluates a polynomial of the form:
%
%   P(x) = a_n*x^n + a_(n-1)*x^(n-1) + ... + a_1*x + a_0
%
% It also evaluates the derivative P'(x) at the same points.
%
% Horner's method reduces computational cost by rewriting the polynomial
% in nested form, avoiding repeated exponentiation.
%
% Inputs:
%   coefficients      - row vector of polynomial coefficients ordered from
%                       highest degree to constant term
%   evaluationPoints  - row vector of x-values where P(x) and P'(x)
%                       are evaluated
%
% Outputs:
%   polynomialValue   - values of P(x) at each evaluation point
%   derivativeValue   - values of P'(x) at each evaluation point

% --- Determine polynomial degree and number of evaluation points ---
numPoints = length(evaluationPoints);

% --- Preallocate output arrays for efficiency ---
polynomialValue = zeros(1, numPoints);
derivativeValue = zeros(1, numPoints);

% --- Evaluate polynomial and derivative at each point ---
for pointIndex = 1:numPoints

    currentX = evaluationPoints(pointIndex);

    % Initialize Horner recursion with leading coefficient
    polynomialAccumulator = coefficients(1);
    derivativeAccumulator = 0;

    % Nested polynomial evaluation
    for coefficientIndex = 2:length(coefficients)

        derivativeAccumulator = derivativeAccumulator * currentX + polynomialAccumulator;
        polynomialAccumulator = polynomialAccumulator * currentX + coefficients(coefficientIndex);

    end

    polynomialValue(pointIndex) = polynomialAccumulator;
    derivativeValue(pointIndex) = derivativeAccumulator;

end

end