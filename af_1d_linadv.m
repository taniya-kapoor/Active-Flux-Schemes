% implementation of 1D linear Advection problem by Active Flux Scheme
% order of Accuracy is three
% formulas for scheme have been taken from eymann's thesis
% initial condition and problem has been taken from maeng's thesis

function arijit
clear all;                                             % clearing workspace
clc;                                                   % clearing workspace
format long;                                          % for better accuracy
hold off;                                    % to not repeat previous curve
xa = 0.0;                              % left limit of computational domain
xb = 1.0;                             % right limit of computational domain
t = 0.0;                                                    % starting time
tend = 1.0;                                   % time at which soln is reqrd
No_Cell = 80;                                             % number of cells 
Delta_x = (xb-xa)/No_Cell;                                             % dx
Cell_Centers = xa-Delta_x/2:Delta_x:xb+Delta_x/2;            % cell centers
Cell_Interfaces = xa-Delta_x:Delta_x:xb+Delta_x;          % cell interfaces
n = size(Cell_Centers,2);                          % number of cell centers
m = size(Cell_Interfaces,2);                    % number of cell interfeces
u0 = (1/(2*pi))*sin(2*pi*Cell_Centers);                 % initial condition
v0 = (1/(2*pi))*sin(2*pi*Cell_Interfaces);              % initial condition
w0 = ( v0(1:m-1) + 4*u0(1:n)  + v0(2:m) )/6;    % constructing cell average
v1 = v0; w1=w0; v=0.7;                  % pre allocating v1 and w1 cfl =0.7

while(t<tend)
    dt = v * Delta_x;                                      % time step size
    dt = min(tend-dt,dt);                      % so that final time is tend
    % cell inrerface values update
    v1(3:m) = v*(3*v-2)*v0(2:m-1) + 6*v*(1-v)*w0(2:n)  + (1-v)*(1-3*v)*v0(3:m);
    v1(1) = v1(m-2); v1(2) = v1(m-1);        % periodic boundary conditions
    % cell average update
    w1(2:n-1) = v^2*(v-1)*v0(1:m-3) + v^2*(3-2*v)*w0(1:n-2) + v*(1-v)*v0(2:m-2) + (1-v)^2*(1+2*v)*w0(2:n-1) - v*(1-v)^2*v0(3:m-1);
    w1(1) = w1(n-1); w1(n) = w1(2);          % periodic boundary conditions
    v0 = v1; w0 = w1;                % for next time step, to continue loop
    
    exact = (1/(2*pi))*sin(2*pi*(Cell_Centers-t-dt));
    
    %ploting of solution and exact solution;
    hold off;
    plot(Cell_Centers,w1,'o')
    grid on;
    hold on;
    plot(Cell_Centers,exact);
    ylabel('u','fontsize', 16)
    xlabel('x','fontsize', 16)
    title(sprintf('time = %f',t),'fontsize',16)
    drawnow;
    t=t+dt;
end
L2=sqrt((1/((xb-xa+2*Delta_x)*(dt))).*sum(abs(w1*Delta_x*dt-exact*Delta_x*dt).^2))
