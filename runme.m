%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation code for the paper "Further results on stability 
% of linear systems with slow and fast time variation and switching,"
% presented at HSCC 2024: 27th ACM International Conference on 
% Hybrid Systems: Computation and Control, Hong Kong, China, May 13-16, 2024
% 
% Authors: Hyungbo Shim and Daniel Liberzon
% 
% This code is used for generating the simulation results of Section I.A,
% and Figure 2 of the paper.
%
% Tested in MATLAB 9.13 (R2020b)
% Version 0.1
% Date: Jan 24, 2024

clear all
close all

% Simulation time: [0,2] seconds
tw = [0, 2];
% Initial condition [x; x_dot; phi; phi_dot]
ic = [-0.2; 0; 0; 0];
% Dynamics of the Inverted Pendulum without perturbation is called by ip0. 
ip0 = @(t,x) ip(t,x,0);
[t0,x0] = ode45(ip0,tw,ic);
% Dynamics of the Inverted Pendulum with perturbation is called by ip1. 
ip1 = @(t,x) ip(t,x,1);
[t1,x1] = ode45(ip1,tw,ic);
% Simulation outcomes are stored in t0, x0, t1, and x1 (time and states).

% To observe the values of L(t), I(t), K(t), and K'(t) [perturbed K], 
% re-run the function ip with the mode 2 (returning internal parameters).
% The parameters are stored in the variable rec.
getvalue = @(t) ip(t,0,2);
N = size(t1,1);
rec = [];
for i=1:N,
    rec = [rec; getvalue(t1(i))];
end

% Plot the figures of the paper
h = figure;
screenSize = get(0, 'ScreenSize');
figWidth = screenSize(3)/2;
figHeight = screenSize(4)/2;
set(h, 'Position', [0, 0, figWidth, figHeight]);

subplot(4,1,1)
h = plot(t1,rec(:,2),'-',t1,rec(:,1),'--');
set(h,'linewidth', 2)
legend('Moment of Inertia','Ball position','fontsize',12)
title('Time varying parameters','fontsize',15,'fontweight','bold')

subplot(4,1,2)
h = plot(t0,x0(:,3),'-',t0,x0(:,1),'--');
set(h,'linewidth', 2)
axis([0,2,-0.25,0.1])
grid on
legend('Angle of Tube','Position of Cart','fontsize',12,'Location','SouthEast')
title('Response of Average System','fontsize',15,'fontweight','bold')

subplot(4,1,3)
h = plot(t1,-rec(:,7),'-',t1,-rec(:,3),'--');
set(h,'linewidth', 2)
legend('Used Gain K_1(t)','Designed Gain K_1(t)','fontsize',12)
title('Comparison between Designed Gain and Used Gain','fontsize',15,'fontweight','bold')

subplot(4,1,4)
h = plot(t1,x1(:,3),'-',t1,x1(:,1),'--');
set(h,'linewidth', 2)
axis([0,2,-0.25,0.1])
grid on
legend('Angle of Tube','Position of Cart','fontsize',12,'Location','SouthEast')
title('Response of Closed-loop System','fontsize',15,'fontweight','bold')
xlabel('time','fontsize',15)


function dxdt = ip(t,x,mode)
%
% states = {'x' 'x_dot' 'phi' 'phi_dot'};
% inputs = {'u'};
% outputs = {'x'; 'phi'};
% mode == 0 [unperturbed] / == 1 [perturbed] / == 2 [return parameters]
%

M = 0.5;
m = 0.2;
b = 0.1;
g = 9.8;
L = 0.1 + mod(t,0.5);
I = m*L^2;
p = (M+m)*I + M*m*L^2; 

A = [0      1              0           0;
     0 -(I+m*L^2)*b/p  (m^2*g*L^2)/p   0;
     0      0              0           1;
     0 -(m*L*b)/p       m*g*L*(M+m)/p  0];

B = [     0;
     (I+m*L^2)/p;
          0;
        m*L/p];

K = place( A, B, [-8.5 + 7.9i;   -8.5 - 7.9i; ...
    -4.8 + 0.8i;  -4.8 - 0.8i] );

Kprime = K + 20*sin(2*pi*20*t)*ones(1,4);

switch mode
    case 0
        dxdt = A*x - B*K*x;
    case 1
        dxdt = A*x - B*Kprime*x;
    case 2
        dxdt = [L, I, K, Kprime];
end

end
