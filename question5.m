clear;
clc;
close all;

u_exact = @(x,y) sin(x).*cos(y);
N = 50;
a = 0;
b = 5;
h = (b-a)/N;

x = linspace(a,b,N+1);
y = linspace(a,b,N+1);

xi = x(2:N);
yi = y(2:N);

[X,Y] = meshgrid(xi,yi);
X = X.';
Y = Y.';

n = N-1;
e = ones(n,1);

Dx = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
Dy = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;

Ix = speye(n);
Iy = speye(n);

L = kron(Iy,Dx) + kron(Dy,Ix);

F = -2*sin(X).*cos(Y);

F(1,:) = F(1,:) - u_exact(0,yi)/h^2;
F(end,:) = F(end,:) - u_exact(5,yi)/h^2;

F(:,1) = F(:,1) - u_exact(xi,0)'/h^2;
F(:,end) = F(:,end) - u_exact(xi,5)'/h^2;

u_vec = L \ F(:);
Uinterior = reshape(u_vec,n,n);
Uinterior = reshape(u_vec,n,n);
U = zeros(N+1,N+1);
U(2:N,2:N) = Uinterior;
U(:,1) = u_exact(x,0)';
U(:,end) = u_exact(x,5)';
U(1,:) = u_exact(0,y);
U(end,:) = u_exact(5,y);

[Xfull,Yfull] = meshgrid(x,y);

Uexact = u_exact(Xfull,Yfull);

[Xfull,Yfull] = meshgrid(x,y);
Xfull = Xfull.';
Yfull = Yfull.';

Uexact = sin(Xfull).*cos(Yfull);

%% Okabe-Ito colorblind-safe colors
blue   = [0 114 178]/255;
orange = [230 159 0]/255;

figure('Theme','light');
hold on;

% Exact solution: smooth surface
s1 = surf(Xfull,Yfull,Uexact, ...
    'FaceColor',orange, ...
    'FaceAlpha',0.75, ...
    'EdgeColor','none');

% Numerical solution: wireframe
s2 = mesh(Xfull,Yfull,U, ...
    'EdgeColor',blue, ...
    'FaceColor','none', ...
    'LineWidth',0.8);

xlabel('$x$', ...
    'Interpreter','latex', ...
    'FontSize',22);

ylabel('$y$', ...
    'Interpreter','latex', ...
    'FontSize',22);

zlabel('$u(x,y)$', ...
    'Interpreter','latex', ...
    'FontSize',22);

title('Exact and Finite Difference Solutions', ...
    'Interpreter','latex', ...
    'FontSize',22);

legend([s1 s2], ...
    {'Exact solution', ...
    'Finite difference solution'}, ...
    'Interpreter','latex', ...
    'Location','northeast', ...
    'FontSize',16);

set(gca,'FontSize',18);

grid on;
box on;

view(45,30);

hold off;

exportgraphics(gcf, 'poisson_solution.pdf', 'ContentType', 'vector');




%% Convergence study

Ns = [10 20 40 80 160];

errors = zeros(size(Ns));
hs = zeros(size(Ns));

u_exact = @(x,y) sin(x).*cos(y);
f_fun = @(x,y) -2*sin(x).*cos(y);

for k = 1:length(Ns)

    N = Ns(k);

    a = 0;
    b = 5;

    h = (b-a)/N;
    hs(k) = h;

    x = linspace(a,b,N+1);
    y = linspace(a,b,N+1);
    xi = x(2:N);
    yi = y(2:N);
    [X,Y] = meshgrid(xi,yi);
    X = X.';
    Y = Y.';

    n = N-1;

    e = ones(n,1);

    Dx = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
    Dy = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;

    Ix = speye(n);
    Iy = speye(n);

    L = kron(Iy,Dx) + kron(Dy,Ix);

    F = f_fun(X,Y);
    F(1,:) = F(1,:) - u_exact(0,yi)/h^2;
    F(end,:) = F(end,:) - u_exact(5,yi)/h^2;
    F(:,1) = F(:,1) - u_exact(xi,0)'/h^2;
    F(:,end) = F(:,end) - u_exact(xi,5)'/h^2;
    u_vec = L \ F(:);
    Uinterior = reshape(u_vec,n,n);

    Uexact = u_exact(X,Y);

    errors(k) = max(abs(Uinterior(:) - Uexact(:)));

end

%% Convergence plot

blue   = [0 114 178]/255;
orange = [230 159 0]/255;

figure('Theme','light');
hold on;

% Numerical error
plot(hs, errors, '-d', ...
    'Color', blue, ...
    'MarkerFaceColor', blue, ...
    'LineWidth', 2, ...
    'MarkerSize', 8);

% O(h^2) reference line
C = errors(1)/hs(1)^2;

plot(hs, C*hs.^2, '--', ...
    'Color', orange, ...
    'LineWidth', 2);

% FORCE logarithmic axes
set(gca, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', 18);

xlabel('$h$', ...
    'Interpreter','latex', ...
    'FontSize',22);

ylabel('$\|u_h-u\|_\infty$', ...
    'Interpreter','latex', ...
    'FontSize',22);

title('Convergence of the Finite Difference Solution', ...
    'Interpreter','latex', ...
    'FontSize',22);

legend( ...
    '$\|u_h-u\|_\infty$', ...
    '$\mathcal{O}(h^2)$', ...
    'Interpreter','latex', ...
    'Location','southeast', ...
    'FontSize',17);

grid on;
box on;

hold off;

exportgraphics(gcf, ...
    'poisson_convergence.pdf', ...
    'ContentType','vector');