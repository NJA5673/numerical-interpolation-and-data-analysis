function [polynomialCoefficients, basisPolynomials] = Lagrange_Interpolation_Method(xData, yData)
% LAGRANGE_INTERPOLATION_METHOD
% Constructs the interpolating polynomial using Lagrange basis functions.
%
% This function computes:
%   1) The full interpolating polynomial coefficients
%   2) Each individual Lagrange basis polynomial contribution
%
% Given data points (x_i, y_i), the interpolating polynomial is:
%
%   P(x) = Σ y_k * L_k(x)
%
% where each L_k(x) is the Lagrange basis polynomial defined as:
%
%   L_k(x) = Π (x - x_i) / (x_k - x_i),   for i ≠ k
%
% Inputs:
%   xData  - vector of x-coordinates of data points
%   yData  - vector of corresponding y-values
%
% Outputs:
%   polynomialCoefficients - coefficients of the full interpolating polynomial
%   basisPolynomials       - matrix where each row represents one Lagrange
%                            basis polynomial contribution

% --- Ensure row vectors ---
xData = xData(:).';
yData = yData(:).';

% --- Validate input ---
if numel(xData) ~= numel(yData)
    error('Input vectors xData and yData must have the same length.');
end

numPoints = numel(xData);
polynomialDegree = numPoints - 1;

% --- Preallocate matrix for basis polynomials ---
basisPolynomials = zeros(numPoints, polynomialDegree + 1);

% --- Construct Lagrange basis polynomials ---
for k = 1:numPoints

    numeratorPolynomial = 1;
    denominatorScalar = 1;

    for i = 1:numPoints
        if i == k
            continue;
        end

        % Multiply polynomial by (x - x_i)
        numeratorPolynomial = conv(numeratorPolynomial, [1, -xData(i)]);

        % Compute denominator product
        denominatorScalar = denominatorScalar * (xData(k) - xData(i));
    end

    % Scale basis polynomial by corresponding y-value
    basisPolynomials(k, :) = (yData(k) / denominatorScalar) * numeratorPolynomial;

end

% --- Combine all basis polynomials ---
polynomialCoefficients = sum(basisPolynomials, 1);

end