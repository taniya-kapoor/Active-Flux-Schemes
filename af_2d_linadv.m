% two dimensional advection code solved by active flux method
% desired order of accuracy is between 2.5 and 3
% solving u_t + u_x + u_y = 0
% ic is gaussian type u(x,y,0) = exp(-0.5*(x^2+y^2))
% bc depending on the advection is applied
% cartesian mesh is used

clear all;                                             % clearing workspace
clc;                                                   % clearing workspace
format long                                           % for better accuracy
xa = -5;                                             % left limit of domain
xb = 5;                                             % right limit of domain
ya = -5;                                            % below limit of domain
yb = 5;                                               % top limit of domain
No_Cell_x = 100;                                      % no of cells along x
No_Cell_y = 100;                                      % no of cells along y
Delta_x = (xb-xa)/No_Cell_x;                                    % step size
Delta_y = (yb-ya)/No_Cell_y;                                    % step size
v = 0.2;                                                       % cfl number
t = 0;                                                         % time start
tend = 1;                                                        % time end
dt = (v*Delta_x*Delta_y)/(Delta_x+Delta_y);                 % cfl condition
lambda = [1;1];                                                % wave speed

dfbx = xa - Delta_x/2 : Delta_x : xb + Delta_x/2;   % descritization for black dot
dfby = ya - Delta_y/2 : Delta_y : yb + Delta_y/2;   % descritization for black dot

dfcx = xa - Delta_x : Delta_x : xb + Delta_x;    % descritization for cross
dfcy = ya - Delta_y : Delta_y : yb + Delta_y;    % descritization for cross

dfsx = xa - 3*Delta_x/2 : Delta_x : xb + Delta_x/2;% descritization for star
dfsy = ya - Delta_y : Delta_y : yb + Delta_y;      % descritization for star

dfwx = xa - Delta_x : Delta_x : xb + Delta_x;       % descritization for white dot
dfwy = ya - 3*Delta_y/2 : Delta_y : yb + Delta_y/2; % descritization for white dot

[X1,Y1] = meshgrid(dfbx,dfby);                              % creating grid
[X2,Y2] = meshgrid(dfcx,dfcy);                              % creating grid
[X3,Y3] = meshgrid(dfsx,dfsy);                              % creating grid
[X4,Y4] = meshgrid(dfwx,dfwy);                              % creating grid

lx1 = length(X1);
lx2 = length(X2);
lx3 = length(X3);
lx4 = length(X4);
%========================= work for IC (verified)==========================
u0 = zeros(lx1,lx1);                                       % pre allocating
v0 = zeros(lx2,lx2);                                       % pre allocating
w0 = zeros(lx3,lx3);                                       % pre allocating
z0 = zeros(lx4,lx4);                                       % pre allocating
exact = zeros(lx1,lx1);                                    % pre allocating

sizeu = size(u0,1)*size(u0,2);
sizev = size(v0,1)*size(v0,2);
sizew = size(w0,1)*size(w0,2);
sizez = size(z0,1)*size(z0,2);

for i=1:sizeu
    u0(i) = exp(-0.5*(X1(i)^2 + Y1(i)^2)); % initial condition on black dot(u)
end

for i=1:sizev
    v0(i) = exp(-0.5*(X2(i)^2 + Y2(i)^2)); % initial condition on cross(v)
end

for i=1:sizew
    w0(i) = exp(-0.5*(X3(i)^2 + Y3(i)^2)); % initial condition on star(w)
end

for i=1:sizez
    z0(i) = exp(-0.5*(X4(i)^2 + Y4(i)^2)); % initial condition on white dot(z)
end

w0 = w0(:,2:end);
z0 = z0(2:end,:);

%initial condition plot
figure(1)
contourf(X1,Y1,u0)
contourcbar

%==========================================================================

% finally stored only the grid points with IC not anything outside

%=================== length of each face (verified)========================
No_of_edges = 4*sizeu;                    % number of edges with repetition
l = Delta_x;                                           % length of one edge
length_of_face = l*ones(4*sizeu,1);           % storing all length together
%==========================================================================

%========================= work for normals (verified)=====================
n1 = zeros(4*sizeu,1);         % calculating first column of normal vectors
n1(1:4:end) = 0;
n1(2:4:end) = 1;
n1(3:4:end) = 0;
n1(4:4:end) = -1;

n2 = zeros(4*sizeu,1);        % calculating second column of normal vectors
n2(1:4:end) = -1;
n2(2:4:end) = 0;
n2(3:4:end) = 1;
n2(4:4:end) = 0;
normals = [n1 n2];                            % normal vector through faces
%==========================================================================


%========================= work for area (verified)========================
ar = Delta_x*Delta_y;                                   % area of each grid
area = ar*ones(sizeu,1);                  % storing all grid areas together
%==========================================================================

jac_inv = [2/Delta_x 0;0 2/Delta_y];         % inverse jacobian calculation
minus_part = dt*jac_inv*lambda;                       % minus part for xi_f

% pre allocating xi_i follows

xi_i_6_x_half = 1;
xi_i_6_y_half = 0;

xi_i_3_x_half = 1;
xi_i_3_y_half = 1;

xi_i_7_x_half = 0;
xi_i_7_y_half = 1;

xi_i_6_x_full = 1;
xi_i_6_y_full = 0;

xi_i_3_x_full = 1;
xi_i_3_y_full = 1;

xi_i_7_x_full = 0;
xi_i_7_y_full = 1;

% calculation of xi_f follows

xi_f_6_x_half = xi_i_6_x_half - minus_part(1)/2;
xi_f_6_y_half = xi_i_6_y_half - minus_part(2)/2;

xi_f_3_x_half = xi_i_3_x_half - minus_part(1)/2;
xi_f_3_y_half = xi_i_3_y_half - minus_part(2)/2;

xi_f_7_x_half = xi_i_7_x_half - minus_part(1)/2;
xi_f_7_y_half = xi_i_7_y_half - minus_part(2)/2;

xi_f_6_x_full = xi_i_6_x_full - minus_part(1);
xi_f_6_y_full = xi_i_6_y_full - minus_part(2);

xi_f_3_x_full = xi_i_3_x_full - minus_part(1);
xi_f_3_y_full = xi_i_3_y_full - minus_part(2);

xi_f_7_x_full = xi_i_7_x_full - minus_part(1);
xi_f_7_y_full = xi_i_7_y_full - minus_part(2);

% allocating values on all four type of points
u0 = u0;
u1 = v0(1:end-1,1:end-1);
u2 = v0(1:end-1,2:end);
u3 = v0(2:end,2:end);
u4 = v0(2:end,1:end-1);
u5 = w0(1:end-1,:);
u6 = z0(:,2:end);
u7 = w0(2:end,:);
u8 = z0(:,1:end-1);

% calculating cell average
u9 = ( u1+u2+u3+u4 + 4*(u5+u6+u7+u8) +16*u0 )/36;

% pre allocating the values at half time step
u1_1 = u1;
u2_1 = u2;
u3_1 = u3;
u4_1 = u4;
u5_1 = u5;
u6_1 = u6;
u7_1 = u7;
u8_1 = u8;

% pre allocating vales at full time step
u1_2 = u1;
u2_2 = u2;
u3_2 = u3;
u4_2 = u4;
u5_2 = u5;
u6_2 = u6;
u7_2 = u7;
u8_2 = u8;
u9_2 = u9;

% loop starts
for ii = t:dt:tend
    % calculating the half time interface update of pt. 6
    u6_1 = 0.25*(xi_f_6_x_half-1)*(xi_f_6_y_half-1)*(xi_f_6_x_half)*(xi_f_6_y_half)*u1 + 0.25*(xi_f_6_x_half+1)*(xi_f_6_y_half-1)*(xi_f_6_x_half)*(xi_f_6_y_half)*u2 + 0.25*(xi_f_6_x_half+1)*(xi_f_6_y_half+1)*(xi_f_6_x_half)*(xi_f_6_y_half)*u3 + 0.25*(xi_f_6_x_half-1)*(xi_f_6_y_half+1)*(xi_f_6_x_half)*(xi_f_6_y_half)*u4 + 0.5*(1-xi_f_6_x_half^2)*(xi_f_6_y_half-1)*(xi_f_6_y_half)*u5 + 0.5*(1-xi_f_6_y_half^2)*(xi_f_6_x_half+1)*(xi_f_6_x_half)*u6 + 0.5*(1-xi_f_6_x_half^2)*(xi_f_6_y_half+1)*(xi_f_6_y_half)*u7 + 0.5*(1-xi_f_6_y_half^2)*(xi_f_6_x_half-1)*(xi_f_6_x_half)*u8 + (1-xi_f_6_x_half^2)*(1-xi_f_6_y_half^2)*u9;
    % calculating the full time interface update of pt. 6
    u6_2 = 0.25*(xi_f_6_x_full-1)*(xi_f_6_y_full-1)*(xi_f_6_x_full)*(xi_f_6_y_full)*u1 + 0.25*(xi_f_6_x_full+1)*(xi_f_6_y_full-1)*(xi_f_6_x_full)*(xi_f_6_y_full)*u2 + 0.25*(xi_f_6_x_full+1)*(xi_f_6_y_full+1)*(xi_f_6_x_full)*(xi_f_6_y_full)*u3 + 0.25*(xi_f_6_x_full-1)*(xi_f_6_y_full+1)*(xi_f_6_x_full)*(xi_f_6_y_full)*u4 + 0.5*(1-xi_f_6_x_full^2)*(xi_f_6_y_full-1)*(xi_f_6_y_full)*u5 + 0.5*(1-xi_f_6_y_full^2)*(xi_f_6_x_full+1)*(xi_f_6_x_full)*u6 + 0.5*(1-xi_f_6_x_full^2)*(xi_f_6_y_full+1)*(xi_f_6_y_full)*u7 + 0.5*(1-xi_f_6_y_full^2)*(xi_f_6_x_full-1)*(xi_f_6_x_full)*u8 + (1-xi_f_6_x_full^2)*(1-xi_f_6_y_full^2)*u9;
    
    % calculating the half time interface update of pt. 3
    u3_1 = 0.25*(xi_f_3_x_half-1)*(xi_f_3_y_half-1)*(xi_f_3_x_half)*(xi_f_3_y_half)*u1 + 0.25*(xi_f_3_x_half+1)*(xi_f_3_y_half-1)*(xi_f_3_x_half)*(xi_f_3_y_half)*u2 + 0.25*(xi_f_3_x_half+1)*(xi_f_3_y_half+1)*(xi_f_3_x_half)*(xi_f_3_y_half)*u3 + 0.25*(xi_f_3_x_half-1)*(xi_f_3_y_half+1)*(xi_f_3_x_half)*(xi_f_3_y_half)*u4 + 0.5*(1-xi_f_3_x_half^2)*(xi_f_3_y_half-1)*(xi_f_3_y_half)*u5 + 0.5*(1-xi_f_3_y_half^2)*(xi_f_3_x_half+1)*(xi_f_3_x_half)*u6 + 0.5*(1-xi_f_3_x_half^2)*(xi_f_3_y_half+1)*(xi_f_3_y_half)*u7 + 0.5*(1-xi_f_3_y_half^2)*(xi_f_3_x_half-1)*(xi_f_3_x_half)*u8 + (1-xi_f_3_x_half^2)*(1-xi_f_3_y_half^2)*u9;
    % calculating the full time interface update of pt. 3
    u3_2 = 0.25*(xi_f_3_x_full-1)*(xi_f_3_y_full-1)*(xi_f_3_x_full)*(xi_f_3_y_full)*u1 + 0.25*(xi_f_3_x_full+1)*(xi_f_3_y_full-1)*(xi_f_3_x_full)*(xi_f_3_y_full)*u2 + 0.25*(xi_f_3_x_full+1)*(xi_f_3_y_full+1)*(xi_f_3_x_full)*(xi_f_3_y_full)*u3 + 0.25*(xi_f_3_x_full-1)*(xi_f_3_y_full+1)*(xi_f_3_x_full)*(xi_f_3_y_full)*u4 + 0.5*(1-xi_f_3_x_full^2)*(xi_f_3_y_full-1)*(xi_f_3_y_full)*u5 + 0.5*(1-xi_f_3_y_full^2)*(xi_f_3_x_full+1)*(xi_f_3_x_full)*u6 + 0.5*(1-xi_f_3_x_full^2)*(xi_f_3_y_full+1)*(xi_f_3_y_full)*u7 + 0.5*(1-xi_f_3_y_full^2)*(xi_f_3_x_full-1)*(xi_f_3_x_full)*u8 + (1-xi_f_3_x_full^2)*(1-xi_f_3_y_full^2)*u9;
    
    % calculating the half time interface update of pt. 7
    u7_1 = 0.25*(xi_f_7_x_half-1)*(xi_f_7_y_half-1)*(xi_f_7_x_half)*(xi_f_7_y_half)*u1 + 0.25*(xi_f_7_x_half+1)*(xi_f_7_y_half-1)*(xi_f_7_x_half)*(xi_f_7_y_half)*u2 + 0.25*(xi_f_7_x_half+1)*(xi_f_7_y_half+1)*(xi_f_7_x_half)*(xi_f_7_y_half)*u3 + 0.25*(xi_f_7_x_half-1)*(xi_f_7_y_half+1)*(xi_f_7_x_half)*(xi_f_7_y_half)*u4 + 0.5*(1-xi_f_7_x_half^2)*(xi_f_7_y_half-1)*(xi_f_7_y_half)*u5 + 0.5*(1-xi_f_7_y_half^2)*(xi_f_7_x_half+1)*(xi_f_7_x_half)*u6 + 0.5*(1-xi_f_7_x_half^2)*(xi_f_7_y_half+1)*(xi_f_7_y_half)*u7 + 0.5*(1-xi_f_7_y_half^2)*(xi_f_7_x_half-1)*(xi_f_7_x_half)*u8 + (1-xi_f_7_x_half^2)*(1-xi_f_7_y_half^2)*u9;
    % calculating the full time interface update of pt. 7
    u7_2 = 0.25*(xi_f_7_x_full-1)*(xi_f_7_y_full-1)*(xi_f_7_x_full)*(xi_f_7_y_full)*u1 + 0.25*(xi_f_7_x_full+1)*(xi_f_7_y_full-1)*(xi_f_7_x_full)*(xi_f_7_y_full)*u2 + 0.25*(xi_f_7_x_full+1)*(xi_f_7_y_full+1)*(xi_f_7_x_full)*(xi_f_7_y_full)*u3 + 0.25*(xi_f_7_x_full-1)*(xi_f_7_y_full+1)*(xi_f_7_x_full)*(xi_f_7_y_full)*u4 + 0.5*(1-xi_f_7_x_full^2)*(xi_f_7_y_full-1)*(xi_f_7_y_full)*u5 + 0.5*(1-xi_f_7_y_full^2)*(xi_f_7_x_full+1)*(xi_f_7_x_full)*u6 + 0.5*(1-xi_f_7_x_full^2)*(xi_f_7_y_full+1)*(xi_f_7_y_full)*u7 + 0.5*(1-xi_f_7_y_full^2)*(xi_f_7_x_full-1)*(xi_f_7_x_full)*u8 + (1-xi_f_7_x_full^2)*(1-xi_f_7_y_full^2)*u9;
    
    % updating pt. 1 half time step
    u1_1(2:end,2:end) = u3_1(1:end-1,1:end-1);
    u1_1(1,:) = u3_1(:,end)'; % bc
    u1_1(:,1) = u3_1(end,:)'; % bc
    
    % updating pt. 1 full time step
    u1_2(2:end,2:end) = u3_2(1:end-1,1:end-1);
    u1_2(1,:) = u3_2(:,end)'; % bc
    u1_2(:,1) = u3_2(end,:)'; % bc
    
    % updating pt. 5 half time step
    u5_1(2:end,:) = u7_1(1:end-1,:);
    u5_1(1,:) = u7_1(:,end)'; % bc
    
    % updating pt. 5 full time step
    u5_2(2:end,:) = u7_2(1:end-1,:);
    u5_2(1,:) = u7_2(:,end)'; % bc
    
    % updating pt. 2 half time step
    u2_1(2:end,:) = u3_1(1:end-1,:);
    u2_1(1,:) = u3_1(:,end)'; % bc
    
    % updating pt. 2 full time step
    u2_2(2:end,:) = u3_2(1:end-1,:);
    u2_2(1,:) = u3_2(:,end)'; % bc
    
    % updating pt. 4 half time step
    u4_1(:,2:end) = u3_1(:,1:end-1);
    u4_1(:,1) = u3_1(end,:)'; % bc
    
    % updating pt. 4 full time step
    u4_2(:,2:end) = u3_2(:,1:end-1);
    u4_2(:,1) = u3_2(end,:)'; % bc
    
    % updating pt. 8 half time step
    u8_1(:,2:end) = u6_1(:,1:end-1);
    u8_1(:,1) = u6_1(end,:)'; % bc
    
    % updating pt. 8 full time step
    u8_2(:,2:end) = u6_2(:,1:end-1);
    u8_2(:,1) = u6_2(end,:)'; % bc
    
    % parameters for flux through face one

    uln1 = u1;
    umn1 = u5;
    urn1 = u2;
    ulnplushalf1 = u1_1; 
    umnplushalf1 = u5_1;
    urnplushalf1 = u2_1;
    ulnplusone1 = u1_2;
    umnplusone1 = u5_2;
    urnplusone1 = u2_2;

    % parameters for flux through face three

    uln3 = u2;
    umn3 = u6;
    urn3 = u3;
    ulnplushalf3 = u2_1;
    umnplushalf3 = u6_1;
    urnplushalf3 = u3_1;
    ulnplusone3 = u2_2;
    umnplusone3 = u6_2;
    urnplusone3 = u3_2;

    % parameters for flux through face two

    uln2 = u3;
    umn2 = u7;
    urn2 = u4;
    ulnplushalf2 = u3_1;
    umnplushalf2 = u7_1;
    urnplushalf2 = u4_1;
    ulnplusone2 = u3_2;
    umnplusone2 = u7_2;
    urnplusone2 = u4_2;

    % parameters for flux through face four

    uln4 = u4;
    umn4 = u8;
    urn4 = u1;
    ulnplushalf4 = u4_1;
    umnplushalf4 = u8_1;
    urnplushalf4 = u1_1;
    ulnplusone4 = u4_2;
    umnplusone4 = u8_2;
    urnplusone4 = u1_2;

    % flux through face 1
    flux1 = (( uln1 + urn1 + ulnplusone1 + urnplusone1 ) + 4*( umn1 + umnplusone1 + ulnplushalf1 + urnplushalf1 ) + 16*umnplushalf1 )/36;

    % flux through face 2
    flux2 = (( uln2 + urn2 + ulnplusone2 + urnplusone2 ) + 4*( umn2 + umnplusone2 + ulnplushalf2 + urnplushalf2 ) + 16*umnplushalf2 )/36;

    % flux through face 3
    flux3 = (( uln3 + urn3 + ulnplusone3 + urnplusone3 ) + 4*( umn3 + umnplusone3 + ulnplushalf3 + urnplushalf3 ) + 16*umnplushalf3 )/36;

    % flux through face 4
    flux4 = (( uln4 + urn4 + ulnplusone4 + urnplusone4 ) + 4*( umn4 + umnplusone4 + ulnplushalf4 + urnplushalf4 ) + 16*umnplushalf4 )/36;
    
    %flux1(1,:) = flux2(:,end);
    %flux4(:,1) = flux3(end,:);
    
    sz = size(flux1,1)*size(flux1,2);
    flux = zeros(4*sz,1);

    flux(1:4:end) = flux1;
    flux(2:4:end) = flux2;
    flux(3:4:end) = flux3;
    flux(4:4:end) = flux4;
    final_flux = [flux flux];

    flux_dot_normals = final_flux.*normals;
    final_dot = flux_dot_normals(:,1) + flux_dot_normals(:,2); 
    under_summation = final_dot.*length_of_face;

    for i=1:sizeu
        total(i) = under_summation(4*i-3) + under_summation(4*i-2) + under_summation(4*i-1) + under_summation(4*i);
    end

    for i=1:sizeu
        u9_2(i) = u9(i) -(dt/area(i))*total(i);
    end

    %u9_2(1,:) = u9_2(:,end)';
    %u9_2(1,:) = u9_2(end,:)';

    u9 = u9_2;
    u1 = u1_2;
    u2 = u2_2;
    u3 = u3_2;
    u4 = u4_2;
    u5 = u5_2;
    u6 = u6_2;
    u7 = u7_2;
    u8 = u8_2;

    hold off;
    figure(2)
    contourf(X1,Y1,u9_2)
    contourcbar
    drawnow;
end

for i=1:sizeu
    exact(i) = exp(-0.5*((X1(i)-1)^2 + (Y1(i)-1)^2)); % initial condition on black dot(u)
end
u9_2;
exact;
forl2 = Delta_x*Delta_y*dt*ones(lx1,lx1);
L2=sqrt((1/((xb-xa+2*Delta_x)*(yb-ya+2*Delta_y)*(dt))).*sum(sum(abs(u9_2.*forl2-exact.*forl2).^2)))
ncell = sizeu;
nface = No_of_edges;
nvertex = sizev;
DOF = ncell + nface/2 + nvertex/4
h = 1/sqrt(DOF)
    




















