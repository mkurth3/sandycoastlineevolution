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
% dy/dt = D*d^2y/dx^2
% where "d" = partial derivative and units are m/s
% y is the height of the coastline, x is the distance along the shore, 
% and t is time while D is diffusivity

%~~~~~~~~~~ Method ~~~~~~~~~~%

% Finite difference methods are used with partial differential equations (PDEs) 
% to discretize derivatives, replacing them with finite differences that can 
% be evaluated at grid points in space and time

%~~~~~~~~~~ Code ~~~~~~~~~~%

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

% Figure 1 & 2, Intial vs. Final Coastline Position and Coastline Evolution over time
% Finite Difference Method
% We chose the Finite Different Method because we wanted to discritize the derivatives in the diffusion equation with finite difference approximations

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

% Figure 3, Effect of Diffusivity on Final Coastline
% Using the same method (Finite Difference) we run the previous model but with varying constant diffusivities.
% Here we compare the final coastline profiles against each other.

figure(3);
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

% As expected, these results show that the higher the diffusivity, the lower/more flat the final coastline profile will be.
% This is due to a higher diffusivity allowing further spreading of materials (sand) most likely due to erosion.
% Higher diffusivity in the example of a sandy coastline is most likely due to stronger waves, smaller sand grain size, or coastal geometry.

% Figures 4-6, Varying Diffusivites Along Distance
% Using the same method just incorporating different patterns of variation in diffusivity along the profile.
% 3 different patterns: sinusoidal (following a sin wave) , half-zero (half of profile is 0 diffusivity and other half is 1),
% and linear gradient (diffusivity is 0 at starting x and increases linearly along profile)
% Our goal here is to see how differing diffusivity patterns affects the coastline evolution.

% define diffusivity variations
diffusivity_patterns = {
    1 + 0.5 * sin(2*pi*x/L), ... % sinusoidal variation
    [zeros(1, Nx/2), ones(1, Nx/2)], ... % half-zero diffusivity
    linspace(1, 2, Nx) ... % linear gradient
};
titles = {'Sinusoidal Diffusivity', 'Half-Zero Diffusivity', 'Linear Gradient Diffusivity'};

% loop through each diffusivity pattern
for pattern_idx = 1:length(diffusivity_patterns)
    D = diffusivity_patterns{pattern_idx};
    
    % initial condition
    y = zeros(Nx, 1);
    y(1:Nx/2) = linspace(0, 100, Nx/2); % linear initial condition
    
    % time-stepping loop
    y_history = zeros(Nx, Nt/100 + 1); % store results for plotting
    history_idx = 1; % time history storage index
    for n = 1:Nt
        y_new = y;
        for i = 2:Nx-1
            % incorporate varying diffusivity
            D_i_plus_half = (D(i) + D(i+1)) / 2; % average diffusivity between i and i+1
            D_i_minus_half = (D(i) + D(i-1)) / 2; % average diffusivity between i and i-1
            % handling zero diffusivity 
            if D(i) > 0
                y_new(i) = y(i) + dt/dx^2 * ...
                    (D_i_plus_half*(y(i+1)-y(i))-D_i_minus_half*(y(i)-y(i-1)));
            else
                y_new(i) = y(i); % no evolution for zero diffusivity
            end
        end
        y = y_new;
        % save results every 100 time steps
        if mod(n, 100) == 0
            y_history(:, history_idx) = y;
            history_idx = history_idx+1;
        end
    end
    % plot diffusivity graph and shoreline evolution
    figure;
    % subplot 1: Diffusivity
    subplot(1, 2, 1);
    plot(x, D, 'b-', 'LineWidth', 1.6);
    xlabel('Distance alongshore (m)'); ylabel('Diffusivity (m^2/s)');
    title(['Diffusivity: ', titles{pattern_idx}]); grid on;
    % subplot 2: Shoreline evolution
    subplot(1, 2, 2);
    hold on;
    for t_idx = 1:size(y_history, 2)
        plot(x, y_history(:, t_idx), 'DisplayName', sprintf('t = %.1f s', (t_idx-1)*100*dt));
    end
    xlabel('Distance alongshore (m)');ylabel('Shoreline position (m)');
    title('Shoreline Evolution'); legend show; grid on;
end

% These results show that the varying diffusivity patterns do severely effect the coastline evolution.
% This is most likely more representative of real life since many factors not included in the model can alter diffusivity along a profile.
% This could be the presence of vegetation, size and angularity of material that is being diffused, water patterns, etc.

%~~~~~Final Results/Model Analysis and Model Limitations/Improvements ~~~~~~%



%~~~~~~~~~~ References ~~~~~~~~~~%

% Nam, P. T., Larson, M., Hanson, H., & Hoan, L. X. (2011). 
% A numerical model of beach morphological evolution due to waves and 
% currents in the vicinity of coastal structures. Coastal Engineering, 
% 58(9), 863–876. https://doi.org/10.1016/j.coastaleng.2011.05.006 

% Partial differential equation and finite difference. 
% Scientific Machine Learning (SciML). (n.d.). 
% https://kks32-courses.github.io/sciml/lectures/04-pde-fdm/04-pde-fdm.html 

% Slingerland, R., & Kump, L. R. (2011).
% Mathematical modeling of Earth’s dynamical systems: A Primer. 
% Princeton University Press. 
