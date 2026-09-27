% Projectile_Trajectory_Modeling.m
% Finite-Difference Analysis of Projectile Motion Data
%
% This script analyzes projectile motion using clean and noisy position data.
% Finite-difference methods are used to estimate velocity and acceleration
% components from discrete measurements. The analysis compares how noise
% affects derivative estimates and physical quantities such as gravity,
% launch speed, and launch angle.

clear; clc; close all;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(projectRoot, 'data');
resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end

%% Load data

cleanData = readmatrix(fullfile(dataDir, 'projectile_data.csv'));
noisyData = readmatrix(fullfile(dataDir, 'projectile_data_noise.csv'));
if size(cleanData, 2) < 3 || size(noisyData, 2) < 3 || ...
        size(cleanData, 1) < 5 || size(noisyData, 1) < 5
    error('Each data file needs at least five rows and three columns: time, x, y.');
end

timeClean = cleanData(:, 1);
xClean = cleanData(:, 2);
yClean = cleanData(:, 3);

timeNoisy = noisyData(:, 1);
xNoisy = noisyData(:, 2);
yNoisy = noisyData(:, 3);

%% Estimate time step

timeStepClean = mean(diff(timeClean));
timeStepNoisy = mean(diff(timeNoisy));
if any(diff(timeClean) <= 0) || any(diff(timeNoisy) <= 0) || ...
        any(abs(diff(timeClean) - timeStepClean) > 1e-8 * timeStepClean) || ...
        any(abs(diff(timeNoisy) - timeStepNoisy) > 1e-8 * timeStepNoisy)
    error('Time values must be increasing and equally spaced.');
end

%% Compute velocity and acceleration using finite differences

vxClean = FirstDerivative_Oh2(xClean, timeStepClean);
vyClean = FirstDerivative_Oh2(yClean, timeStepClean);

axClean = SecondDerivative_Oh2(xClean, timeStepClean);
ayClean = SecondDerivative_Oh2(yClean, timeStepClean);

vxNoisy = FirstDerivative_Oh2(xNoisy, timeStepNoisy);
vyNoisy = FirstDerivative_Oh2(yNoisy, timeStepNoisy);

axNoisy = SecondDerivative_Oh2(xNoisy, timeStepNoisy);
ayNoisy = SecondDerivative_Oh2(yNoisy, timeStepNoisy);

%% Figure 1: Projectile trajectory

figure('Name', 'Projectile Trajectory Comparison', 'Position', [100 100 900 600]);

plot(xClean, yClean, 'LineWidth', 1.6); hold on;
plot(xNoisy, yNoisy, '--', 'LineWidth', 1.3);

grid on;
xlabel('Horizontal Position, x (m)', 'FontSize', 12);
ylabel('Vertical Position, y (m)', 'FontSize', 12);
title('Projectile Trajectory: Clean Data vs Noisy Data', 'FontSize', 14);
legend({'Clean Data', 'Noisy Data'}, 'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'projectile_trajectory_comparison.png'), 'Resolution', 300);

%% Figure 2: Velocity components

figure('Name', 'Projectile Velocity Components', 'Position', [100 100 900 700]);

subplot(2, 1, 1);
plot(timeClean, vxClean, 'LineWidth', 1.3); hold on;
plot(timeNoisy, vxNoisy, '--', 'LineWidth', 1.1);
grid on;
xlabel('Time, t (s)');
ylabel('v_x (m/s)');
title('Horizontal Velocity Component');
legend({'Clean Data', 'Noisy Data'}, 'Location', 'best');

subplot(2, 1, 2);
plot(timeClean, vyClean, 'LineWidth', 1.3); hold on;
plot(timeNoisy, vyNoisy, '--', 'LineWidth', 1.1);
grid on;
xlabel('Time, t (s)');
ylabel('v_y (m/s)');
title('Vertical Velocity Component');
legend({'Clean Data', 'Noisy Data'}, 'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'projectile_velocity_components.png'), 'Resolution', 300);

%% Figure 3: Acceleration components

figure('Name', 'Projectile Acceleration Components', 'Position', [100 100 900 700]);

subplot(2, 1, 1);
plot(timeClean, axClean, 'LineWidth', 1.3); hold on;
plot(timeNoisy, axNoisy, '--', 'LineWidth', 1.1);
grid on;
xlabel('Time, t (s)');
ylabel('a_x (m/s^2)');
title('Horizontal Acceleration Component');
legend({'Clean Data', 'Noisy Data'}, 'Location', 'best');

subplot(2, 1, 2);
plot(timeClean, ayClean, 'LineWidth', 1.3); hold on;
plot(timeNoisy, ayNoisy, '--', 'LineWidth', 1.1);
grid on;
xlabel('Time, t (s)');
ylabel('a_y (m/s^2)');
title('Vertical Acceleration Component');
legend({'Clean Data', 'Noisy Data'}, 'Location', 'best');

exportgraphics(gcf, fullfile(resultsDir, 'projectile_acceleration_components.png'), 'Resolution', 300);

%% Summary calculations

interiorClean = 3:numel(timeClean)-2;
interiorNoisy = 3:numel(timeNoisy)-2;

meanAxClean = mean(axClean(interiorClean));
meanAyClean = mean(ayClean(interiorClean));

meanAxNoisy = mean(axNoisy(interiorNoisy));
meanAyNoisy = mean(ayNoisy(interiorNoisy));

gravityEstimateClean = -meanAyClean;
gravityEstimateNoisy = -meanAyNoisy;

initialVxClean = vxClean(1);
initialVyClean = vyClean(1);

initialVxNoisy = vxNoisy(1);
initialVyNoisy = vyNoisy(1);

launchSpeedClean = hypot(initialVxClean, initialVyClean);
launchAngleClean = atan2d(initialVyClean, initialVxClean);

launchSpeedNoisy = hypot(initialVxNoisy, initialVyNoisy);
launchAngleNoisy = atan2d(initialVyNoisy, initialVxNoisy);

%% Print summary table

fprintf('\n================ PROJECTILE MOTION SUMMARY ================\n');
fprintf('%-15s %-18s %-18s %-15s %-15s %-15s\n', ...
    'Data Source', 'Mean ax (m/s^2)', 'Mean ay (m/s^2)', ...
    'g Estimate', 'v0 (m/s)', 'theta0 (deg)');
fprintf('--------------------------------------------------------------------------\n');

fprintf('%-15s %+18.6f %+18.6f %15.6f %15.4f %15.3f\n', ...
    'Clean Data', meanAxClean, meanAyClean, ...
    gravityEstimateClean, launchSpeedClean, launchAngleClean);

fprintf('%-15s %+18.6f %+18.6f %15.6f %15.4f %15.3f\n', ...
    'Noisy Data', meanAxNoisy, meanAyNoisy, ...
    gravityEstimateNoisy, launchSpeedNoisy, launchAngleNoisy);

fprintf('==========================================================================\n');

%% Local functions

function derivative = FirstDerivative_Oh2(signal, h)
% FIRSTDERIVATIVE_OH2
% Computes a second-order accurate finite-difference estimate of the first
% derivative for uniformly spaced data.
%
% Centered differences are used for interior points. One-sided second-order
% formulas are used near the boundaries.

numPoints = numel(signal);
derivative = zeros(numPoints, 1);

% Forward one-sided estimates near beginning
derivative(1) = (-3*signal(1) + 4*signal(2) - signal(3)) / (2*h);
derivative(2) = (signal(3) - signal(1)) / (2*h);

% Centered estimates for interior points
for index = 3:numPoints-2
    derivative(index) = (signal(index+1) - signal(index-1)) / (2*h);
end

% Backward one-sided estimates near end
derivative(numPoints-1) = (signal(numPoints) - signal(numPoints-2)) / (2*h);
derivative(numPoints) = (3*signal(numPoints) - 4*signal(numPoints-1) + signal(numPoints-2)) / (2*h);

end

function secondDerivative = SecondDerivative_Oh2(signal, h)
% SECONDDERIVATIVE_OH2
% Computes a second-order accurate finite-difference estimate of the second
% derivative for uniformly spaced data.
%
% Centered differences are used for interior points. One-sided second-order
% formulas are used near the boundaries.

numPoints = numel(signal);
secondDerivative = zeros(numPoints, 1);

% Forward one-sided estimates near beginning
secondDerivative(1) = (2*signal(1) - 5*signal(2) + 4*signal(3) - signal(4)) / h^2;
secondDerivative(2) = (signal(3) - 2*signal(2) + signal(1)) / h^2;

% Centered estimates for interior points
for index = 3:numPoints-2
    secondDerivative(index) = ...
        (signal(index+1) - 2*signal(index) + signal(index-1)) / h^2;
end

% Backward one-sided estimates near end
secondDerivative(numPoints-1) = ...
    (signal(numPoints) - 2*signal(numPoints-1) + signal(numPoints-2)) / h^2;

secondDerivative(numPoints) = ...
    (2*signal(numPoints) - 5*signal(numPoints-1) + 4*signal(numPoints-2) - signal(numPoints-3)) / h^2;

end
