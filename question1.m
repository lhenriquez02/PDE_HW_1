clear;
clc;
close all;

%% Part a: Approximate derivative

N = 100;

a = 0;
b = 2*pi;
h = (b-a)/N;
x = linspace(a,b,N+1);

f = exp(sin(x));
df = exp(sin(x)).*cos(x);

df_approx = zeros(size(x));

df_approx(1) = (f(2)-f(1))/h;
df_approx(N+1) = (f(N+1)-f(N))/h;
df_approx(2:N) = (f(3:N+1)-f(1:N-1))/(2*h);


%% Plot exact and approximate derivative

fs = 24;

cols = [ ...
      0 114 178;
    230 159   0
    ] / 255;

fig = figure( ...
    'Theme', 'light', ...
    'Position', [100, 100, 1000, 1000]);

hold on

plot(x, df, '-', ...
    'Color', cols(1,:), ...
    'LineWidth', 3, ...
    'DisplayName', 'Exact derivative');

plot(x(1:3:end), df_approx(1:3:end), '--o', ...
    'Color', cols(2,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', cols(2,:), ...
    'MarkerFaceColor', 'white', ...
    'DisplayName', 'Finite difference');

ax = gca;
set(ax, ...
    'LineWidth', 2, ...
    'FontSize', fs, ...
    'TickLabelInterpreter', 'latex');

grid on
box on

xlabel('$x$', 'Interpreter','latex','FontSize',fs);
ylabel('$f''(x)$', 'Interpreter','latex','FontSize',fs);

title('Derivative of $f(x)=e^{\sin(x)}$', ...
    'Interpreter','latex', ...
    'FontSize',fs);

legend( ...
    'Location','best', ...
    'Interpreter','latex', ...
    'FontSize',0.75*fs);

xlim([0 2*pi]);

xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'});

hold off

exportgraphics(fig, 'finiteDifferenceQ1a.pdf', ...
    'ContentType','vector');


%% Part b: Relative error convergence

Ns = 2.^(3:10);

errInf = zeros(size(Ns));
errL2  = zeros(size(Ns));

for k = 1:length(Ns)

    N = Ns(k);

    a = 0;
    b = 2*pi;
    dx = (b-a)/N;
    x = linspace(a,b,N+1);

    f  = exp(sin(x));
    df = exp(sin(x)).*cos(x);

    df_approx = zeros(size(x));

    df_approx(1)   = (f(2)-f(1))/dx;
    df_approx(N+1) = (f(N+1)-f(N))/dx;
    df_approx(2:N) = ...
        (f(3:N+1)-f(1:N-1))/(2*dx);

    e = df-df_approx;

    errInf(k) = norm(e,inf)/norm(df,inf);
    errL2(k)  = norm(e,2)/norm(df,2);

end


%% Convergence plot

fig = figure( ...
    'Theme','light', ...
    'Position',[100 100 1200 800]);

hold on

loglog(Ns,errInf,'-o', ...
    'Color',cols(1,:), ...
    'LineWidth',2, ...
    'MarkerSize',9, ...
    'MarkerFaceColor','white', ...
    'DisplayName','$L^\infty$ error');

loglog(Ns,errL2,'--s', ...
    'Color',cols(2,:), ...
    'LineWidth',2, ...
    'MarkerSize',9, ...
    'MarkerFaceColor','white', ...
    'DisplayName','$L^2$ error');

ax = gca;
set(ax, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'LineWidth', 2, ...
    'FontSize', fs, ...
    'TickLabelInterpreter', 'latex');

grid on

xlabel('$N$','Interpreter','latex','FontSize',fs);
ylabel('Relative error','Interpreter','latex','FontSize',fs);

title('Finite Difference Convergence', ...
    'Interpreter','latex', ...
    'FontSize',fs);

legend( ...
    'Location','northeast', ...
    'Interpreter','latex', ...
    'FontSize',0.75*fs);

hold off

exportgraphics(fig,'convergenceq1.pdf', ...
    'ContentType','vector');