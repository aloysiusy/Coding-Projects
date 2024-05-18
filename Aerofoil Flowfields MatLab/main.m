clear; clc; close all;

% inputs
code = input ("Input the NACA 4 digit code: ");
N = input ("Input the number of panels: ");
U_inf = input ("Input the free-stream velocity: ");
alpha = input ("Input the angle of attack: ");


if code == 2412
    N = [50, 100, 200];
    alpha = [0:10];
    for i=1:length(N)
        parfor j = 1:length(alpha)
            CL = strength(code, N(i), U_inf, alpha(j));
            CLplot(i,j) = CL;
        end
    end
    
    xfnaca = importdata("xfnaca.txt");
    xfCL = xfnaca.data(:,2);   
    xfalpha = xfnaca.data(:,1);   
    
    figure
    plot(xfalpha, xfCL)
    hold on
    for k = 1:length(N)
        plot(alpha, CLplot(k,:))
    end
    
    

else
    CL = strength(code, N, U_inf, alpha);
end 