% solving the 1d acoustics equation for simple wave 
% problem is 1d acoustic with a0 = 1; rho = 1/2;
% ic is p'(x,0) = 1/4 + sin(2*pi*x)/80 ; u'(x,0) = 1/4 - sin(2*pi*x)/10

% implementation of 1D linear Advection problem by Active Flux Scheme
% order of Accuracy is three
% formulas for scheme have been taken from eymann's thesis
% initial condition and problem has been taken from maeng's thesis
% wave speed is negative

function one_d_acoustic_2
clear all;                                             % clearing workspace
clc;                                                   % clearing workspace
format long;                                          % for better accuracy
hold off;                                    % to not repeat previous curve
xa = 0.0;                              % left limit of computational domain
xb = 2.0;                             % right limit of computational domain
t = 0.0;                                                    % starting time
tend = 76.8;                                  % time at which soln is reqrd
No_Cell = 50;                                             % number of cells 
Delta_x = (xb-xa)/No_Cell;                                             % dx
Cell_Centers = xa-Delta_x/2:Delta_x:xb+Delta_x/2;            % cell centers
Cell_Interfaces = xa-Delta_x:Delta_x:xb+Delta_x;          % cell interfaces
n = size(Cell_Centers,2);                          % number of cell centers
m = size(Cell_Interfaces,2);                    % number of cell interfeces
u0 = (-(1/4 + sin(2*pi*Cell_Centers)/80) + 1/4 - sin(pi*Cell_Centers)/10)/sqrt(2);                 % initial condition
v0 = (-(1/4 + sin(2*pi*Cell_Interfaces)/80) + 1/4 - sin(pi*Cell_Interfaces)/10)/sqrt(2);             % initial condition
w0 = ( v0(1:m-1) + 4*u0(1:n)  + v0(2:m) )/6;    % constructing cell average

U0 = ((1/4 + sin(2*pi*Cell_Centers)/80) + 1/4 - sin(pi*Cell_Centers)/10)/sqrt(2);                 % initial condition
V0 = ((1/4 + sin(2*pi*Cell_Interfaces)/80) + 1/4 - sin(pi*Cell_Interfaces)/10)/sqrt(2);            % initial condition
W0 = ( V0(1:m-1) + 4*U0(1:n)  + V0(2:m) )/6;    % constructing cell average

v1 = v0; w1=w0; v=0.7;                  % pre allocating v1 and w1 cfl =0.7
V1 = V0; W1=W0;

while(t<tend)
    dt = v * Delta_x;                                      % time step size
    dt = min(tend-dt,dt);                      % so that final time is tend
    % cell inrerface values update
    v1(1:m-2) = (1-v)*(1-3*v)*v0(1:m-2) + 6*v*(1-v)*w0(1:n-1)  + v*(3*v-2)*v0(2:m-1);
    v1(m) = v1(3); v1(m-1) = v1(2);        % periodic boundary conditions
    % cell average update
    w1(2:n-1) = -v*(1-v)^2*v0(2:m-2) + (1-v)^2*(1+2*v)*w0(2:n-1) + v*(1-v)*v0(3:m-1) + v^2*(3-2*v)*w0(3:n) + v^2*(v-1)*v0(4:m);
    w1(1) = w1(n-1); w1(n) = w1(2);          % periodic boundary conditions
    v0 = v1; w0 = w1;                % for next time step, to continue loop
    
    
    % cell inrerface values update
    V1(3:m) = v*(3*v-2)*V0(2:m-1) + 6*v*(1-v)*W0(2:n)  + (1-v)*(1-3*v)*V0(3:m);
    V1(1) = V1(m-2); V1(2) = V1(m-1);        % periodic boundary conditions
    % cell average update
    W1(2:n-1) = v^2*(v-1)*V0(1:m-3) + v^2*(3-2*v)*W0(1:n-2) + v*(1-v)*V0(2:m-2) + (1-v)^2*(1+2*v)*W0(2:n-1) - v*(1-v)^2*V0(3:m-1);
    W1(1) = W1(n-1); W1(n) = W1(2);          % periodic boundary conditions
    V0 = V1; W0 = W1;                % for next time step, to continue loop
    
    exact_neg = (-(1/4 + sin(2*pi*(Cell_Centers+t+dt))/80) + 1/4 - sin(pi*(Cell_Centers+t+dt))/10)/sqrt(2);
    exact_pos = ((1/4 + sin(2*pi*(Cell_Centers-t-dt))/80) + 1/4 - sin(pi*(Cell_Centers-t-dt))/10)/sqrt(2);;
    
    pressure = (-w1 + W1)/(2*sqrt(2));
    velocity = (w1 + W1)/sqrt(2);
    
    exact_pressure = (-exact_neg + exact_pos)/(2*sqrt(2));
    exact_velocity = (exact_neg + exact_pos)/sqrt(2);
    
    %ploting of solution and exact solution;
    figure(1)
    hold off;
    plot(Cell_Centers,pressure,'o')
    grid on;
    hold on;
    plot(Cell_Centers,exact_pressure);
    ylabel('p','fontsize', 16)
    xlabel('x','fontsize', 16)
    title(sprintf('time = %f',t),'fontsize',16)
    drawnow;
    
    figure(2)
    hold off;
    plot(Cell_Centers,velocity,'o')
    grid on;
    hold on;
    plot(Cell_Centers,exact_velocity);
    ylabel('u','fontsize', 16)
    xlabel('x','fontsize', 16)
    title(sprintf('time = %f',t),'fontsize',16)
    drawnow;
    t=t+dt;
end
L2_pressure=sqrt((1/((xb-xa+2*Delta_x)*(dt))).*sum(abs(pressure*Delta_x*dt-exact_pressure*Delta_x*dt).^2))
L2_velocity=sqrt((1/((xb-xa+2*Delta_x)*(dt))).*sum(abs(velocity*Delta_x*dt-exact_velocity*Delta_x*dt).^2))