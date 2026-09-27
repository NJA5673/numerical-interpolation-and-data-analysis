function [numIterations, convergenceTable] = Numerical_Differentiation( ...
    targetFunction, evaluationPoint, tolerance, maxIterations, initialStepSize)
% NUMERICAL_DIFFERENTIATION
% Estimates the derivative of a function using a centered finite-difference
% approximation with second-order accuracy.
%
% The step size is reduced by a factor of 10 at each iteration to study
% convergence behavior:
%
%   f'(x0) ≈ [f(x0 + h) - f(x0 - h)] / (2h)
%
% This method demonstrates the tradeoff between truncation error and
% round-off error as h becomes increasingly small.
%
% Inputs:
%   targetFunction   - function handle for f(x)
%   evaluationPoint  - point x0 where derivative is estimated
%   tolerance        - stopping tolerance for approximate relative error
%   maxIterations    - maximum number of step-size refinements
%   initialStepSize  - starting value of h
%
% Outputs:
%   numIterations    - number of iterations performed
%   convergenceTable - table with columns:
%                      [step size h, derivative estimate, approximate error]

if nargin < 5
    initialStepSize = 1;
end

convergenceTable = zeros(maxIterations, 3);
previousDerivative = NaN;

for iteration = 1:maxIterations

    stepSize = initialStepSize / 10^(iteration - 1);

    derivativeEstimate = ...
        (targetFunction(evaluationPoint + stepSize) - ...
         targetFunction(evaluationPoint - stepSize)) / (2 * stepSize);

    if iteration == 1
        approximateRelativeError = NaN;
    else
        denominator = max(abs(derivativeEstimate), eps);
        approximateRelativeError = ...
            abs(derivativeEstimate - previousDerivative) / denominator;

        if approximateRelativeError <= tolerance
            convergenceTable(iteration, :) = ...
                [stepSize, derivativeEstimate, approximateRelativeError];

            numIterations = iteration;
            convergenceTable = convergenceTable(1:iteration, :);
            return;
        end
    end

    convergenceTable(iteration, :) = ...
        [stepSize, derivativeEstimate, approximateRelativeError];

    previousDerivative = derivativeEstimate;

end

numIterations = maxIterations;
convergenceTable = convergenceTable(1:numIterations, :);

end