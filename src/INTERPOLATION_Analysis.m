% Interpolation_Analysis.m
% Polynomial Interpolation and Function Approximation Study
%
% This script analyzes polynomial interpolation using the Lagrange method.
% Several functions are reconstructed from discrete nodes, and the resulting
% interpolating polynomials are compared against the original functions.
%
% The analysis demonstrates:
% - Construction of interpolating polynomials from sampled data
% - Accuracy of interpolation within a bounded interval
% - Effect of polynomial degree on stability
% - High-order oscillation behavior near interval boundaries

clear; clc; close all;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end

%% Figure formatting

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Verification Case: cosine interpolation

sampleX = [0, 0.6, 1.2];
sampleY = cos(sampleX);

[cosineCoefficients, cosineBasisTerms] = Lagrange_Interpolation_Method(sampleX, sampleY);

evaluationGrid = linspace(min(sampleX) - 0.2, max(sampleX) + 0.2, 800);
cosineInterpolant = Horner_Polynomial_Evaluation(cosineCoefficients, evaluationGrid);

figure('Name', 'Cosine Interpolation Verification', 'Position', [100 100 900 600]);

plot(evaluationGrid, cos(evaluationGrid), 'LineWidth', 1.3); hold on;
plot(evaluationGrid, cosineInterpolant, '--', 'LineWidth', 1.3);
plot(sampleX, sampleY, 'o', 'MarkerFaceColor', [0.1 0.1 0.1]);

grid on;
xlabel('$x$');
ylabel('$y$');
title('Verification: Lagrange Interpolation of $\cos(x)$');
legend({'$\cos(x)$', 'Interpolating Polynomial $P_2$', 'Interpolation Nodes'}, ...
    'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'cosine_interpolation_verification.png'), 'Resolution', 300);

%% Basis contribution visualization

figure('Name', 'Cosine Basis Polynomial Contributions', 'Position', [100 100 900 600]);
hold on;

for termIndex = 1:size(cosineBasisTerms, 1)
    basisContribution = Horner_Polynomial_Evaluation(cosineBasisTerms(termIndex, :), evaluationGrid);
    plot(evaluationGrid, basisContribution, 'LineWidth', 1.0);
end

plot(evaluationGrid, cosineInterpolant, 'k', 'LineWidth', 1.6);

grid on;
xlabel('$x$');
ylabel('$y$');
title('Lagrange Basis Contributions for $\cos(x)$');
legend({'Term 1', 'Term 2', 'Term 3', '$P_2$'}, 'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'cosine_basis_contributions.png'), 'Resolution', 300);

%% Linear data consistency test

linearX = 1:10;
linearY = 3*linearX + 2;

[linearCoefficients, ~] = Lagrange_Interpolation_Method(linearX, linearY);

linearGrid = linspace(1, 10, 800);
linearInterpolant = Horner_Polynomial_Evaluation(linearCoefficients, linearGrid);

figure('Name', 'Linear Data Interpolation', 'Position', [100 100 900 600]);

plot(linearGrid, 3*linearGrid + 2, 'LineWidth', 1.3); hold on;
plot(linearGrid, linearInterpolant, '--', 'LineWidth', 1.3);
plot(linearX, linearY, 'o', 'MarkerFaceColor', [0.1 0.1 0.1]);

grid on;
xlabel('$x$');
ylabel('$y$');
title('Interpolation Consistency Test: Linear Data');
legend({'Original Line $3x+2$', 'Interpolating Polynomial', 'Nodes'}, ...
    'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'linear_interpolation_consistency.png'), 'Resolution', 300);

%% Eighth-order interpolation on smooth functions

nodes = linspace(0, 2, 9);
evaluationGrid = linspace(0, 2, 1200);

smoothCases = { ...
    struct('label', '$e^x$', 'fileTag', 'exp_interpolation', ...
           'functionHandle', @(x) exp(x)); ...
    struct('label', '$\sin(x)$', 'fileTag', 'sin_interpolation', ...
           'functionHandle', @(x) sin(x)); ...
    struct('label', '$(x+1)^{x+1}$', 'fileTag', 'power_interpolation', ...
           'functionHandle', @(x) (x+1).^(x+1)) ...
};

for caseIndex = 1:numel(smoothCases)

    currentCase = smoothCases{caseIndex};
    targetFunction = currentCase.functionHandle;

    [coefficients, basisTerms] = Lagrange_Interpolation_Method(nodes, targetFunction(nodes));
    interpolant = Horner_Polynomial_Evaluation(coefficients, evaluationGrid);

    figure('Name', ['Smooth Function Interpolation: ', currentCase.fileTag], ...
        'Position', [100 100 900 600]);

    plot(evaluationGrid, targetFunction(evaluationGrid), 'LineWidth', 1.3); hold on;
    plot(evaluationGrid, interpolant, '--', 'LineWidth', 1.3);
    plot(nodes, targetFunction(nodes), 'o', 'MarkerFaceColor', [0.1 0.1 0.1]);

    grid on;
    xlabel('$x$');
    ylabel('$y$');
    title(['Eighth-Order Lagrange Interpolation of ', currentCase.label]);
    legend({'Original Function', '$P_8$', 'Nodes'}, 'Location', 'best');

    exportgraphics(gcf, fullfile(resultsDir, sprintf('%s_overlay.png', currentCase.fileTag)), ...
        'Resolution', 300);

    figure('Name', ['Basis Contributions: ', currentCase.fileTag], ...
        'Position', [100 100 1000 650]);
    hold on;

    for termIndex = 1:size(basisTerms, 1)
        basisContribution = Horner_Polynomial_Evaluation(basisTerms(termIndex, :), evaluationGrid);
        plot(evaluationGrid, basisContribution, 'LineWidth', 1.0);
    end

    plot(evaluationGrid, interpolant, 'k', 'LineWidth', 1.6);

    grid on;
    xlabel('$x$');
    ylabel('$y$');
    title(['Lagrange Basis Contributions for ', currentCase.label]);

    legendEntries = arrayfun(@(k) sprintf('Term %d', k), ...
        1:size(basisTerms, 1), 'UniformOutput', false);
    legendEntries{end+1} = '$P_8$';

    legend(legendEntries, 'Location', 'bestoutside');

    exportgraphics(gcf, fullfile(resultsDir, sprintf('%s_basis_contributions.png', currentCase.fileTag)), ...
        'Resolution', 300);

end

%% High-order interpolation stability study

rungeFunction = @(x) exp(-x) .* sin(5*pi*x) .* sin(2*pi*x);
polynomialOrders = [2 5 8 11 14 17 20];

evaluationGrid = linspace(0, 1, 1600);

figure('Name', 'High-Order Interpolation Stability Study', ...
    'Position', [100 100 1200 800]);

tileLayout = tiledlayout(2, 4);
title(tileLayout, 'High-Order Polynomial Interpolation Stability Study');

for orderIndex = 1:numel(polynomialOrders)

    polynomialOrder = polynomialOrders(orderIndex);
    interpolationNodes = linspace(0, 1, polynomialOrder + 1);
    nodeValues = rungeFunction(interpolationNodes);

    coefficients = Lagrange_Interpolation_Method(interpolationNodes, nodeValues);
    interpolant = Horner_Polynomial_Evaluation(coefficients, evaluationGrid);

    nexttile;
    plot(evaluationGrid, interpolant, 'LineWidth', 1.0); hold on;
    plot(interpolationNodes, nodeValues, 'o', 'MarkerFaceColor', [0.1 0.1 0.1]);

    grid on;
    xlabel('$x$');
    ylabel('$y$');
    title(sprintf('$n = %d$', polynomialOrder));

end

nexttile;
plot(evaluationGrid, rungeFunction(evaluationGrid), 'LineWidth', 1.3);
grid on;
xlabel('$x$');
ylabel('$y$');
title('Original Function');

exportgraphics(gcf, fullfile(resultsDir, 'high_order_interpolation_stability.png'), 'Resolution', 300);

%% Final overlay for n = 20

finalOrder = 20;
finalNodes = linspace(0, 1, finalOrder + 1);
finalNodeValues = rungeFunction(finalNodes);

finalCoefficients = Lagrange_Interpolation_Method(finalNodes, finalNodeValues);
finalInterpolant = Horner_Polynomial_Evaluation(finalCoefficients, evaluationGrid);

figure('Name', 'High-Order Interpolation Boundary Oscillation', ...
    'Position', [100 100 900 600]);

plot(evaluationGrid, rungeFunction(evaluationGrid), 'LineWidth', 1.4); hold on;
plot(evaluationGrid, finalInterpolant, '--', 'LineWidth', 1.4);
plot(finalNodes, finalNodeValues, 'o', 'MarkerFaceColor', [0.1 0.1 0.1]);

grid on;
xlabel('$x$');
ylabel('$y$');
title('High-Order Lagrange Interpolation: $n = 20$');
legend({'Original Function', '$P_{20}$', 'Nodes'}, 'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'high_order_n20_overlay.png'), 'Resolution', 300);