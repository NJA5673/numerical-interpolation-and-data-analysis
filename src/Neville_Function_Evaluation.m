function [interpolatedValue, approximateError, iterationsUsed] = Neville_Function_Evaluation( ...
    targetFunction, evaluationPoint, domainBounds, tolerance, maxDivisions)
% NEVILLE_FUNCTION_EVALUATION
% Adaptive function approximation using Neville interpolation.
%
% This function repeatedly refines the interpolation grid over a specified
% domain until the approximate relative error falls below a desired
% tolerance. For each refinement level, Neville's method is used to estimate
% the function value at the selected evaluation point.
%
% Inputs:
%   targetFunction     - function handle for f(x)
%   evaluationPoint    - point where f(x) is approximated
%   domainBounds       - interval [a, b] used to generate interpolation nodes
%   tolerance          - stopping criterion for approximate relative error
%   maxDivisions       - maximum number of interpolation refinements
%
% Outputs:
%   interpolatedValue  - final Neville interpolation estimate
%   approximateError   - approximate relative error between refinements
%   iterationsUsed     - number of refinements required

% --- Domain setup ---
leftBound = domainBounds(1);
rightBound = domainBounds(2);

% Initialize output values
approximateError = NaN;
interpolatedValue = NaN;
iterationsUsed = maxDivisions;
previousEstimate = NaN;

% --- Adaptive refinement loop ---
for numDivisions = 1:maxDivisions

    % Generate equally spaced interpolation nodes
    dataPointsX = linspace(leftBound, rightBound, numDivisions + 1);
    dataPointsY = targetFunction(dataPointsX);

    % Evaluate interpolation estimate using Neville method
    [~, currentEstimate] = Neville_Interpolation_Method( ...
        dataPointsX, dataPointsY, evaluationPoint);

    interpolatedValue = currentEstimate;

    % Compare complete interpolants on successive grids.
    if numDivisions >= 2
        denominator = max(abs(currentEstimate), eps);
        approximateError = abs((currentEstimate - previousEstimate) / denominator);

        if approximateError <= tolerance
            iterationsUsed = numDivisions;
            return;
        end

    end

    previousEstimate = currentEstimate;

end

end
