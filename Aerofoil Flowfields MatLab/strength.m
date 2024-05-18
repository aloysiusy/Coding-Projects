function [CL] = strength(code, N, U_inf, alpha)


[x, z] = panelgen(code, N);


    




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


end

