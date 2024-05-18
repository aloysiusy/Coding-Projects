clear; clc; close all;



% Add the wake
x(end+1,1) = x(end,1) + 5000*cosd(alpha);
z(end+1,1) = z(end,1) + 5000*sind(alpha);

% Formulate panel angles beta
for k = 1:length(x)-1
    beta(k,1) = rad2deg( atan2(z(k+1)-z(k), x(k+1)-x(k)) ); % in degrees
end

% Formulate B matrix
B = zeros(N+1,1);
for k = 1:length(B)-1
    B(k,1) = -U_inf * sind(alpha - beta(k,1));
end
% Kutta condition automatically added (last zero)

% Formulate A matrix
A = zeros(N+1,N+1);
for i = 1:N+1
    for j = 1:N+1
        midpt = [(x(i+1)+x(i))/2, (z(i+1)+z(i))/2];
        plot(midpt(1), midpt(2), 'b.');
        [uij, vij] = cdoublet(midpt, [x(j) z(j)], [x(j+1) z(j+1)]);
        A(i,j) = vij*cosd(beta(i,1)) - uij*sind(beta(i,1));
    end
end
A(end,:) = [1 zeros(1, N-2) -1 1];

% Calculate the circulation intensity
mu = A \ B;

% Calculate lift coefficient
CL = -2*mu(end) / U_inf;


%
% Plot velocity vectors and streamlines

% Set graph limits
xlim([-0.2 1.2])
ylim([-0.7 0.7])

% Create a grid to compute the u and v velocities at those grid points
size = 200;
gx = linspace(-0.5, 1.5, size);
gz = linspace(-1, 1, size);
[Gx, Gz] = meshgrid(gx, gz);



% Calculate the induced velocity by the doublet on the grid points
u_grid = zeros(length(gx), length(gz)); % Velocity grid for all the grid points
v_grid = zeros(length(gx), length(gz)); % Velocity grid for all the grid points

for i = 1:length(gx)
    for j = 1:length(gz)
        for k = 1:N+1
            xc = gx(i); zc = gz(j); % set coordinates of interest
            
            [u,v] = cdoublet([xc, zc], [x(k) z(k)], [x(k+1) z(k+1)]);
            u_grid(i,j) = u_grid(i,j) + mu(k)*u;
            v_grid(i,j) = v_grid(i,j) + mu(k)*v;
        end
    end
end

% Calculate the total velocity (add free stream velocity)
u_grid = u_grid + U_inf*cosd(alpha);
v_grid = v_grid + U_inf*sind(alpha);

u_grid = u_grid';
v_grid = v_grid';

%Check whether points are in polygon
% Make velocity = 0 if points are in polygon
[in, on] = inpolygon(Gx, Gz, x, z);
for i = 1:length(in)
    for j = 1:length(in)
        if in(i,j) == 1 || on(i,j) == 1
            u_grid(i,j) = 0;
            v_grid(i,j) = 0;
        end
    end
end

%
hold on;
% quiver(Gx, Gz, u_grid, v_grid)

% [startX, startY] = meshgrid(-0.2, -0.7:0.01:0.7);
streamslice(Gx, Gz, u_grid, v_grid);
% streamline(Gx, Gz, u_grid, v_grid, startX, startY);

% Xfoil = table2array(readtable('xf-naca2412-il-1000000'));
% figure
% plot(xfoil.alpha, xfoil.CL)


