function [x, z] = panelgen(airfoil_code, N)

% This function calculates the x and z coordinates of the panels
% Inputs: Airfoil 4-digit NACA code, number of panels N
% Outputs: x and z endpoint coordinates of panels + plot visualising the
% panels


% User input of airfoil code
% airfoil_code = input('Input airfoil code: ');

% Convert to char to split
if class(airfoil_code) == "double"
    airfoil_code = num2str(airfoil_code);
elseif class(airfoil_code) == "string"
    airfoil_code = convertStringsToChars(airfoil_code);
end

% Split the numbers into the different components
m = str2double(airfoil_code(1)) / 100; % Maximum camber m
p = str2double(airfoil_code(2)) / 10 ; % Location of max camber p
t = str2double(airfoil_code(3:end)) / 100; % Max thickness t

% User input - Number of panels to be used
% N = input('Input number of panels to be used: ');

% Calculate panel end points, x_end
% N panels would have N+1 end points
for i = 1:N+1
    x_end(i,1) = 1 - 0.5*(1 - cos(2*pi*(i-1)/N));
end

% Calculate yc and yt, and dyc/dx
for j = 1:length(x_end)
    if x_end(j,1) >=0 && x_end(j,1) <= p
        yc(j,1) = m/p^2 * (2*p*x_end(j,1) - x_end(j,1)^2);
        dyc(j,1) = 2*m/p^2 * (p - x_end(j,1));
    elseif x_end(j,1) >= p && x_end(j,1) <= 1
        yc(j,1) = m/((1-p)^2) * ( (1-2*p) + 2*p*x_end(j,1) - x_end(j,1)^2 );
        dyc(j,1) = 2*m / (1-p)^2 * (p - x_end(j,1));
    end
    
    yt(j,1) = 5*t*(0.2969*sqrt(x_end(j,1)) - 0.126*x_end(j,1) - 0.3516*x_end(j,1)^2 + 0.2843*x_end(j,1)^3 - 0.1015*x_end(j,1)^4);
    
end

% Calculate theta
for j = 1:length(yc)
    theta(j,1) = atan(dyc(j,1));
end

% Calculate x and z coordinates
for j = 1:length(x_end)
    if j >= length(x_end)/2 % For lower points
        x(j,1) = x_end(j,1) - yt(j,1)*sin(theta(j,1));
        z(j,1) = yc(j,1) + yt(j,1)*cos(theta(j,1));
    else % For upper points
        x(j,1) = x_end(j,1) + yt(j,1)*sin(theta(j,1));
        z(j,1) = yc(j,1) - yt(j,1)*cos(theta(j,1));
    end    
end


% Plot the panels for visualisation
figure();
plot(x, z, 'r-x'); grid on;
xlim([-0.01 1.01])
ylim([-0.5 0.5])
end

