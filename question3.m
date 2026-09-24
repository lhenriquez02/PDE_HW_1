clear;
clc;
close all;

u = @(x) sin(x);
f = @(x) sin(x).*cos(x);

Ns = [10 20 40 80 160 320];

errors = zeros(size(Ns));
hs = zeros(size(Ns));

for k = 1:length(Ns)

    N = Ns(k);
    h = 5/N;
    hs(k) = h;

    x = linspace(0,5,N+1)';

    xi = x(2:end-1);

    n = N - 1;

    a = 1/h^2 - sin(xi)/(2*h);    % U_{i-1}
    b = 1 - 2/h^2;                % U_i
    c = 1/h^2 + sin(xi)/(2*h);    % U_{i+1}

    A = diag(b*ones(n,1)) ...
      + diag(a(2:end),-1) ...
      + diag(c(1:end-1),1);

    rhs = f(xi);

    u0 = 0;
    uN = sin(5);

    rhs(1)   = rhs(1)   - a(1)*u0;
    rhs(end) = rhs(end) - c(end)*uN;

    Uinterior = A\rhs;

    U = [u0; Uinterior; uN];

    Uexact = u(x);

    errors(k) = norm(U - Uexact, inf);

end


%% Convergence plot

% Okabe-Ito colorblind-safe colors
cols = [ ...
      0 114 178;    % blue
    230 159   0 ... % orange
    ] / 255;

fs = 24;

% Make sure everything is a column vector
hs     = hs(:);
errors = errors(:);

fig = figure( ...
    'Theme', 'light', ...
    'Position', [100 100 1000 700]);

hold on

% ---------------------------------------------------------
% Numerical error
% ---------------------------------------------------------
loglog(hs, errors, '-o', ...
    'Color', cols(1,:), ...
    'LineWidth', 2.5, ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', cols(1,:), ...
    'MarkerFaceColor', cols(2,:), ...
    'DisplayName', 'Finite difference error');

% ---------------------------------------------------------
% Second-order reference line
% Scale it so it appears near the numerical error
% ---------------------------------------------------------
C = 1.5*errors(1) / hs(1)^2;
reference = C * hs.^2;

loglog(hs, reference, '--d', ...
    'Color', 'k', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', '$\mathcal{O}(\Delta x^2)$');


% ---------------------------------------------------------
% Formatting
% ---------------------------------------------------------
ax = gca;

set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'LineWidth', 1.5, ...
    'FontSize', fs, ...
    'TickLabelInterpreter', 'latex');

grid on
box on

xlabel('$\Delta x$', ...
    'Interpreter', 'latex', ...
    'FontSize', fs);

ylabel('$\|U-u\|_\infty$', ...
    'Interpreter', 'latex', ...
    'FontSize', fs);

title('Convergence Plot', ...
    'Interpreter', 'latex', ...
    'FontSize', fs);

legend('Location', 'southeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 0.8*fs);

axis tight

hold off

exportgraphics(fig, 'convergencePlot.pdf', ...
    'ContentType', 'vector');