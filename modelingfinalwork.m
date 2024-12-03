%% Final Modeling Project Draft Notebook
% Marcus & Angela
 
% Evolution of a Sandy Coastline
%~~~~~~~~~~ Introduction ~~~~~~~~~~%

% We are modeling the evolution of a sandy coastline, which is important to
% study due to the current global rising sea level that leaves coastlines
% vulnerable to further erosion and flooding. Currently, according to Pham 
% Nam, there are 6 numerical modeling methods to study coastline evolution:
% (1) conceptual models, (2) shoreline evolution models, (3) profile 
% evolution models, (4) 2D horizontal morphological evolution models, 
% (5) quasi-3D morphological models, and (6) fully 3D morphological models.
% Here we use a combination of 2D and profile shoreline evolution models. 

%~~~~~~~~~~ Equations ~~~~~~~~~~%

% 1-D diffusion PDE: dy/dt - D*d^2y/dx^2 = 0 
% where "d" = partial derivative and units are m/s

% setting parameters
L = 1000; % length of the domain (m)
T = 10000; % simulation time (s)
D = 1; % diffusivity (m^2/s)
Nx = 100; % number of spatial points
Nt = 500; % number of time steps
dx = L/(Nx -1); % Spatial step size
dt = T/Nt; % time step size

% initial conditions
x = linspace(0, L, Nx);
y = zeros(Nx, 1);
y(1:Nx/2) = linspace(0, 100, Nx/2); % linear initial condition
%% 
% Figure 1 & 2, Intial vs. Final Coastline Position and Evolution over time
% 
% Finite Difference Method (do we need more than one method? how many?)

% Initial vs Final Coastline
figure(1);
plot(x, y, 'k--', 'LineWidth', 1.6, 'DisplayName', 'Initial');
hold on;

for n = 1:Nt
    y_new = y;
    for i = 2:Nx-1
        y_new(i) = y(i) + D * dt / dx^2 * (y(i+1) - 2*y(i) + y(i-1));
    end
    y = y_new;
    % save the final state for comparison
    if n == Nt
        plot(x, y, 'r', 'LineWidth', 1.6, 'DisplayName', 'Final');
    end
end
xlabel('Distance alongshore (m)'); ylabel('Coastline position (m)');
title('Initial vs Final Coastline Position'); legend; grid on;

% store shoreline positions at selected time steps for plotting
time_steps_to_plot = [1, 100, 200, 300, Nt];
y_history = zeros(Nx, length(time_steps_to_plot));
plot_index = 1;
% Reset y for evolution plotting
y = zeros(Nx, 1);
y(1:Nx/2) = linspace(0, 100, Nx/2);

for n = 1:Nt
    y_new = y;
    for i = 2:Nx-1
        y_new(i) = y(i)+D*dt/dx^2*(y(i+1)-2*y(i)+y(i-1));
    end
    y = y_new;
    if ismember(n, time_steps_to_plot)
        y_history(:, plot_index) = y;
        plot_index = plot_index + 1;
    end
end

% Plot selected time steps
figure(2);
hold on;
for i = 1:length(time_steps_to_plot)
    plot(x, y_history(:, i), 'DisplayName', sprintf('t = %.1f s', time_steps_to_plot(i) * dt));
end
xlabel('Distance alongshore (m)'); ylabel('Coastline position (m)');
title('Coastline Evolution Over Time'); legend show; grid on;
%%
% Figure 3, Error Analysis of Intial vs. Final Coastline Position and Evolution over time
%

figure(3);
% Analytical solution function
function y_analytical = analytical_solution(Nx, t, D, y_initial)
    % Analytical solution for the 1D diffusion equation
    % Inputs:
    % Nx - number of spatial points
    % t - time
    % D - diffusion coefficient

    % Ensure t > 0 to avoid division by zero
    if t == 0
        y_analytical = y_initial(Nx); % At t = 0, solution is the initial condition
        return;
    end

    % Analytical solution: Gaussian initial condition
    y_analytical = zeros(size(Nx));
    for i = 1:length(Nx)
        % integrate using numerical approximation (trapz = trapezoid)
        Nx_prime = linspace(min(Nx), max(Nx), 1000); % discretize x' for integration
        gauss = exp(-(Nx(i) - Nx_prime).^2 / (4 * D * t)) / sqrt(4 * pi * D * t);
        y_analytical(i) = trapz(Nx_prime, gauss .* y_initial(Nx_prime));
    end
end

% Reset numerical solution
y = zeros(Nx, 1);
y(1:Nx/2) = linspace(0, 100, Nx/2); % linear initial condition

% Numerical solution for final time
for n = 1:Nt
    y_new = y;
    for i = 2:Nx-1
        y_new(i) = y(i) + D * dt / dx^2 * (y(i+1) - 2*y(i) + y(i-1));
    end
    y = y_new;
end

% Analytical solution
t = T; % final time
y_initial = @(x) (x <= L/2) .* (200 * x / L) + (x > L/2) .* (200 * (L - x) / L);
analytical_y = analytical_solution(x, t, D, y_initial);

% Error calculation
absolute_error = abs(y - analytical_y);
relative_error = absolute_error/analytical_y;

% Plot numerical vs analytical solution and error
plot(x, analytical_y, 'b-', 'LineWidth', 1.6, 'DisplayName', 'Analytical Solution');
hold on;
plot(x, y, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Numerical Solution');
plot(x, relative_error, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Error');
xlabel('Distance alongshore (m)');
ylabel('Coastline position (m)');
title('Numerical vs Analytical Solution and Error');
legend({'Analytical Solution', 'Numerical Solution', 'Error'}, 'Location', 'best');
%% 
% Figure 4, Effect of Diffusivity on Final Coastline
% 
% same methods

figure(4);
D_values = [0.5, 1, 2]; % different diffusivities
hold on;
for D_test = D_values
    y = zeros(Nx, 1);
    y(1:Nx/2) = linspace(0, 100, Nx/2); % reset initial condition
    for n = 1:Nt
        y_new = y;
        for i = 2:Nx-1
            y_new(i) = y(i)+D_test*dt/dx^2*(y(i+1)-2*y(i)+y(i-1));
        end
        y = y_new;
    end
    plot(x, y, 'DisplayName', sprintf('D = %.1f', D_test));
end
xlabel('Distance alongshore (m)'); ylabel('Coastline position (m)');
title('Effect of Diffusivity on Final Coastline'); legend show; grid on;

%~~~~~~~~~~ References ~~~~~~~~~~%

% Nam, P. T., Larson, M., Hanson, H., & Hoan, L. X. (2011). 
% A numerical model of beach morphological evolution due to waves and 
% currents in the vicinity of coastal structures. Coastal Engineering, 
% 58(9), 863–876. https://doi.org/10.1016/j.coastaleng.2011.05.006 
