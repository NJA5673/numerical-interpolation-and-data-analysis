% Reproducible synthetic projectile data for the trajectory analysis.
% These are simulated positions, not experimental observations.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(projectRoot, 'data');
if ~isfolder(dataDir), mkdir(dataDir); end

rng(132);
time = (0:0.05:4).';
launchSpeed = 25;
launchAngle = 40;
gravity = 9.81;
x = launchSpeed*cosd(launchAngle)*time;
y = launchSpeed*sind(launchAngle)*time - 0.5*gravity*time.^2;
cleanData = [time, x, y];
noisyData = [time, x + 0.02*randn(size(time)), ...
    y + 0.02*randn(size(time))];

writematrix(cleanData, fullfile(dataDir, 'projectile_data.csv'));
writematrix(noisyData, fullfile(dataDir, 'projectile_data_noise.csv'));
fprintf('Wrote synthetic demo data to %s\n', dataDir);
