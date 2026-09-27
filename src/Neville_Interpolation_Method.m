function [interpolationTable, interpolatedValue] = Neville_Interpolation_Method(dataPointsX, dataPointsY, evaluationPoint)
% NEVILLE_INTERPOLATION_METHOD
% Performs polynomial interpolation using Neville's algorithm.
%
% Neville's method constructs an interpolating polynomial recursively using
% a tabular approach. It avoids explicitly computing polynomial coefficients
% and instead directly evaluates the interpolated value at a desired point.
%
% Given data points (x_i, y_i), the method builds a triangular table where:
%
%   Q(i,1) = y_i
%   Q(i,j) = [(x - x_(i-j+1)) * Q(i,j-1) - (x - x_i) * Q(i-1,j-1)]
%            / (x_i - x_(i-j+1))
%
% The final interpolated value is Q(n,n).
%
% Inputs:
%   dataPointsX     - vector of known x-values
%   dataPointsY     - vector of corresponding y-values
%   evaluationPoint - x-value where interpolation is performed
%
% Outputs:
%   interpolationTable - Neville table (triangular matrix)
%   interpolatedValue  - interpolated value at evaluationPoint

% --- Ensure inputs are row vectors ---
dataPointsX = dataPointsX(:).';
dataPointsY = dataPointsY(:).';

% --- Validate input sizes ---
if numel(dataPointsX) ~= numel(dataPointsY)
    error('Input vectors must have the same length.');
end

numPoints = numel(dataPointsX);

% --- Initialize Neville table ---
interpolationTable = zeros(numPoints, numPoints);

% First column = given function values
interpolationTable(:,1) = dataPointsY(:);

% --- Construct table recursively ---
for rowIndex = 2:numPoints
    for columnIndex = 2:rowIndex

        numerator = ...
            (evaluationPoint - dataPointsX(rowIndex - columnIndex + 1)) * interpolationTable(rowIndex, columnIndex - 1) ...
          - (evaluationPoint - dataPointsX(rowIndex)) * interpolationTable(rowIndex - 1, columnIndex - 1);

        denominator = dataPointsX(rowIndex) - dataPointsX(rowIndex - columnIndex + 1);

        interpolationTable(rowIndex, columnIndex) = numerator / denominator;

    end
end

% --- Final interpolated value ---
interpolatedValue = interpolationTable(numPoints, numPoints);

end
