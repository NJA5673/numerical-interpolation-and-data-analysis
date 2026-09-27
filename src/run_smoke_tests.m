% Small checks for the reusable numerical functions.
addpath(fileparts(mfilename('fullpath')));

[p, dp] = Horner_Polynomial_Evaluation([2, -3, 1], [-1, 0, 2]);
assert(max(abs(p - [6, 1, 3])) < 1e-12);
assert(max(abs(dp - [-7, -3, 5])) < 1e-12);

[coefficients, ~] = Lagrange_Interpolation_Method([0, 1, 2], [1, 0, 3]);
assert(max(abs(polyval(coefficients, [0, 1, 2]) - [1, 0, 3])) < 1e-12);
[~, estimate] = Neville_Interpolation_Method([0, 1, 2], [1, 0, 3], 1.5);
assert(abs(estimate - polyval(coefficients, 1.5)) < 1e-12);
[estimate, ~, ~] = Neville_Function_Evaluation(@(x) x.^2, 1.5, [0, 2], 1e-10, 10);
assert(abs(estimate - 2.25) < 1e-10);

[~, table] = Numerical_Differentiation(@(x) x.^2, 2, 1e-8, 8, 0.1);
assert(abs(table(end, 2) - 4) < 1e-8);

t = (0:31)'/32;
signal = 2 + 3*cos(2*pi*t) - sin(2*pi*t);
[a0, a, b, omega, fitted, rss] = Sinusoidal_Model_Fitting(t, signal, 1);
assert(abs(a0 - 2) < 1e-10 && abs(a - 3) < 1e-10);
assert(abs(b + 1) < 1e-10 && abs(omega - 2*pi) < 1e-10);
assert(max(abs(fitted(:) - signal)) < 1e-10 && rss < 1e-18);

disp('All smoke tests passed.');
