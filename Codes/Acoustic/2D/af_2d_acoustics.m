% two dimensional acoustics code solved by active flux method
% desired order of accuracy is 3
% solving p_t + div(v) = 0; v_t + grad(p) = 0;
% IC is p(x,y,0) = [sin(2*pi*x) + sin(2*pi*y)]/c ; u(x,y,0) = 0 = v(x,y,0)
% BC periodic (direct what goes out from right enters through left and,
% what goes out from top enters through bottom )
% cartesian mesh is used

clear all;                                             % clearing workspace
clc;                                                   % clearing workspace
format short                                          % for better accuracy
xa = -1;                                             % left limit of domain
xb = 1;                                             % right limit of domain
ya = -1;                                            % below limit of domain
yb = 1;                                               % top limit of domain
No_Cell_x = 50;                                       % no of cells along x
No_Cell_y = 50;                                       % no of cells along y
Delta_x = (xb-xa)/No_Cell_x;                                    % step size
Delta_y = (yb-ya)/No_Cell_y;                                    % step size
v = 0.95;                                                      % cfl number
t = 0;                                                         % time start
tend = 1/8;                                                       % time end
dt = Delta_x/2;                                             % cfl condition

dfbx = xa - Delta_x/2 : Delta_x : xb + Delta_x/2;   % descritization for black dot
dfby = ya - Delta_y/2 : Delta_y : yb + Delta_y/2;   % descritization for black dot

dfcx = xa - Delta_x : Delta_x : xb + Delta_x;    % descritization for cross
dfcy = ya - Delta_y : Delta_y : yb + Delta_y;    % descritization for cross

dfsx = xa - 3*Delta_x/2 : Delta_x : xb + Delta_x/2; % descritization for star
dfsy = ya - Delta_y : Delta_y : yb + Delta_y;       % descritization for star

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
u10 = zeros(lx1,lx1);                                       % pre allocating pressure
v10 = zeros(lx2,lx2);                                       % pre allocating pressure
w10 = zeros(lx3,lx3);                                       % pre allocating pressure
z10 = zeros(lx4,lx4);                                       % pre allocating pressure

u20 = zeros(lx1,lx1);                                       % pre allocating vel1
v20 = zeros(lx2,lx2);                                       % pre allocating vel1
w20 = zeros(lx3,lx3);                                       % pre allocating vel1
z20 = zeros(lx4,lx4);                                       % pre allocating vel1

u30 = zeros(lx1,lx1);                                       % pre allocating vel2
v30 = zeros(lx2,lx2);                                       % pre allocating vel2
w30 = zeros(lx3,lx3);                                       % pre allocating vel2
z30 = zeros(lx4,lx4);                                       % pre allocating vel2

exactp = zeros(lx1,lx1);                                    % pre allocating exact pressure
exactvel1 = zeros(lx1,lx1);                                 % pre allocating exact vel1
exactvel2 = zeros(lx1,lx1);                                 % pre allocating exact vel2

sizeu = size(u10,1)*size(u10,2);
sizev = size(v10,1)*size(v10,2);
sizew = size(w10,1)*size(w10,2);
sizez = size(z10,1)*size(z10,2);

% initial condition for pressure

for i=1:sizeu
    u10(i) = sin(2*pi*X1(i)) + sin(2*pi*Y1(i)); % initial condition on black dot(u) pressure
end

for i=1:sizev
    v10(i) = sin(2*pi*X2(i)) + sin(2*pi*Y2(i)); % initial condition on cross(v) pressure
end

for i=1:sizew
    w10(i) = sin(2*pi*X3(i)) + sin(2*pi*Y3(i)); % initial condition on star(w) pressure
end

for i=1:sizez
    z10(i) = sin(2*pi*X4(i)) + sin(2*pi*Y4(i)); % initial condition on white dot(z) pressure
end

w10 = w10(:,2:end);
z10 = z10(2:end,:);

% initial condition for velocity 1

for i=1:sizeu
    u20(i) = 0; % initial condition on black dot(u) vel1
end

for i=1:sizev
    v20(i) = 0; % initial condition on cross(v) vel1
end

for i=1:sizew
    w20(i) = 0; % initial condition on star(w) vel1
end

for i=1:sizez
    z20(i) = 0; % initial condition on white dot(z) vel1
end

w20 = w20(:,2:end);
z20 = z20(2:end,:);

% initial condition for velocity 2

for i=1:sizeu
    u30(i) = 0; % initial condition on black dot(u) vel2
end

for i=1:sizev
    v30(i) = 0; % initial condition on cross(v) vel2
end

for i=1:sizew
    w30(i) = 0; % initial condition on star(w) vel2
end

for i=1:sizez
    z30(i) = 0; % initial condition on white dot(z) vel2
end

w30 = w30(:,2:end);
z30 = z30(2:end,:);

%initial condition plot
%figure(1)
%contourf(X1,Y1,u0)
%contourcbar

figure(1)
surf(X1,Y1,u10)
%view(2)
colorbar

figure(2)
surf(X1,Y1,u20)
%view(2)
colorbar

figure(3)
surf(X1,Y1,u30)
%view(2)
colorbar
%==========================================================================

% finally stored only the grid points with IC not anything outside

%=================== length of each face (verified)========================
No_of_edges = 4*sizeu;                    % number of edges with repetition
length_of_face = Delta_x;                          % storing length of face
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

% allocating values on all four type of points for pressure
press0 = u10; 
press1 = v10(1:end-1,1:end-1); 
press2 = w10(1:end-1,:);
press3 = v10(1:end-1,2:end);
press4 = z10(:,2:end);
press5 = v10(2:end,2:end);
press6 = w10(2:end,:);
press7 = v10(2:end,1:end-1);
press8 = z10(:,1:end-1); 

% calculating cell average for pressure
press9 = 9*( 4*press0 + ( press1 + press3 + press5 + press7 -4*press2 -4*press4 -4*press6 -4*press8)/3 )/16 ;

% allocating values on all four type of points for velocity 1
vel10 = u20; 
vel11 = v20(1:end-1,1:end-1); 
vel12 = w20(1:end-1,:);
vel13 = v20(1:end-1,2:end);
vel14 = z20(:,2:end);
vel15 = v20(2:end,2:end);
vel16 = w20(2:end,:);
vel17 = v20(2:end,1:end-1);
vel18 = z20(:,1:end-1); 

% calculating cell average for velocity 1
vel19 = 9*( 4*vel10 + ( vel11 + vel13 + vel15 + vel17 -4*vel12 -4*vel14 -4*vel16 -4*vel18)/3 )/16 ;

% allocating values on all four type of points for velocity 2
vel20 = u30; 
vel21 = v30(1:end-1,1:end-1); 
vel22 = w30(1:end-1,:);
vel23 = v30(1:end-1,2:end);
vel24 = z30(:,2:end);
vel25 = v30(2:end,2:end);
vel26 = w30(2:end,:);
vel27 = v30(2:end,1:end-1);
vel28 = z30(:,1:end-1); 

% calculating cell average for velocty 2
vel29 = 9*( 4*vel20 + ( vel21 + vel23 + vel25 + vel27 -4*vel22 -4*vel24 -4*vel26 -4*vel28)/3 )/16 ;

% pre allocating the values at half time step for pressure
press1_1 = press1;
press2_1 = press2;
press3_1 = press3;
press4_1 = press4;
press5_1 = press5;
press6_1 = press6;
press7_1 = press7;
press8_1 = press8;

% pre allocating vales at full time step for pressure
press1_2 = press1;
press2_2 = press2;
press3_2 = press3;
press4_2 = press4;
press5_2 = press5;
press6_2 = press6;
press7_2 = press7;
press8_2 = press8;
press9_2 = press9;

% pre allocating the values at half time step for velocity 1
vel11_1 = vel11;
vel12_1 = vel12;
vel13_1 = vel13;
vel14_1 = vel14;
vel15_1 = vel15;
vel16_1 = vel16;
vel17_1 = vel17;
vel18_1 = vel18;

% pre allocating vales at full time step for velocity 1
vel11_2 = vel11;
vel12_2 = vel12;
vel13_2 = vel13;
vel14_2 = vel14;
vel15_2 = vel15;
vel16_2 = vel16;
vel17_2 = vel17;
vel18_2 = vel18;
vel19_2 = vel19;

% pre allocating the values at half time step for velocity 2
vel21_1 = vel21;
vel22_1 = vel22;
vel23_1 = vel23;
vel24_1 = vel24;
vel25_1 = vel25;
vel26_1 = vel26;
vel27_1 = vel27;
vel28_1 = vel28;

% pre allocating vales at full time step for velocity 2
vel21_2 = vel21;
vel22_2 = vel22;
vel23_2 = vel23;
vel24_2 = vel24;
vel25_2 = vel25;
vel26_2 = vel26;
vel27_2 = vel27;
vel28_2 = vel28;
vel29_2 = vel29;

% writing values obtained from mathematica === starts 

x11vf = 0.235179;
x12vf = 0.494801;
x13vf = 0.235179;
x14vf = 0.494801;
x15vf = 0.235179;
x16vf = 0.494801;
x17vf = 0.235179;
x18vf = 0.494801;

x11vh = 0.242545;
x12vh = 0.49745;
x13vh = 0.242545;
x14vh = 0.49745;
x15vh = 0.242545;
x16vh = 0.49745;
x17vh = 0.242545;
x18vh = 0.49745;


x21vf = 0.00608009;
x22vf = 0.0126123;
x23vf = 0.00608009;
x24vf = 0;
x25vf = -0.00608009;
x26vf = -0.0126123;
x27vf = -0.0608009;
x28vf = 0;

x21vh = 0.00311118;
x22vh = 0.00633646;
x23vh = 0.00311118;
x24vh = 0;
x25vh = -0.00311118;
x26vh = -0.00633646;
x27vh = -0.00311118;
x28vh = 0;


x31vf = 0.00608009;
x32vf = 0;
x33vf = -0.00608009;
x34vf = -0.0126123;
x35vf = -0.00608009;
x36vf = 0;
x37vf = 0.00608009;
x38vf = 0.0126123;

x31vh = 0.00311118;
x32vh = 0;
x33vh = -0.00311118;
x34vh = -0.00633646;
x35vh = -0.00311118;
x36vh = 0;
x37vh = 0.00311118;
x38vh = 0.00633646;


y21vf = -1.49678;
y22vf = -2.94747;
y23vf = -1.49678;
y24vf = -2.94106;
y25vf = -1.49678;
y26vf = -2.94747;
y27vf = -1.49678;
y28vf = -2.94106;

y21vh = -1.74188;
y22vh = -3.46058;
y23vh = -1.74188;
y24vh = -3.4573;
y25vh = -1.74188;
y26vh = -3.46058;
y27vh = -1.74188;
y28vh = -3.4573;

y31vf = -0.953769;
y32vf = 0;
y33vf = 0.953769;
y34vf = 0;
y35vf = -0.953769;
y36vf = 0;
y37vf = 0.953769;
y38vf = 0;

y31vh = -1.10936;
y32vh = 0;
y33vh = 1.10936;
y34vh = 0;
y35vh = -1.10936;
y36vh = 0;
y37vh = 1.10936;
y38vh = 0;


y41vf = -0.992916;
y42vf = -1.96111;
y43vf = -0.992916;
y44vf = -1.96111;
y45vf = -0.992916;
y46vf = -1.96111;
y47vf = -0.992916;
y48vf = -1.96111;

y41vh =-1.15877; 
y42vh = -2.30511;
y43vh = -1.15877;
y44vh = -2.30511;
y45vh = -1.15877;
y46vh = -2.30511;
y47vh = -1.15877;
y48vh = -2.30511;


z31vf = -1.49678;
z32vf = -2.94106;
z33vf = -1.49678;
z34vf = -2.94106;
z35vf =-1.49678 ;
z36vf =-2.94106 ;
z37vf = -1.49678;
z38vf = -2.94747;

z31vh = -1.74188;
z32vh = -3.4573;
z33vh = -1.74188;
z34vh = -3.46058;
z35vh = -1.74188;
z36vh = -3.4573;
z37vh = -1.74188;
z38vh = -3.46058;


y11vf = x21vf;
y12vf = x22vf;
y13vf = x23vf;
y14vf = x24vf;
y15vf = x25vf;
y16vf = x26vf;
y17vf = x27vf;
y18vf = x28vf;

y11vh = x21vh;
y12vh = x22vh;
y13vh = x23vh;
y14vh = x24vh;
y15vh = x25vh;
y16vh = x26vh;
y17vh = x27vh;
y18vh = x28vh;


z11vf = x31vf;
z12vf = x32vf;
z13vf = x33vf;
z14vf = x34vf;
z15vf = x35vf;
z16vf = x36vf;
z17vf = x37vf;
z18vf = x38vf;

z11vh = x31vh;
z12vh = x32vh;
z13vh = x33vh;
z14vh = x34vh;
z15vh = x35vh;
z16vh = x36vh;
z17vh = x37vh;
z18vh = x38vh;


z21vf = y31vf;
z22vf = y32vf;
z23vf = y33vf;
z24vf = y34vf;
z25vf = y35vf;
z26vf = y36vf;
z27vf = y37vf;
z28vf = y38vf;

z21vh = y31vh;
z22vh = y32vh;
z23vh = y33vh;
z24vh = y34vh;
z25vh = y35vh;
z26vh = y36vh;
z27vh = y37vh;
z28vh = y38vh;


z41vf = y41vf;
z42vf = y42vf;
z43vf = y43vf;
z44vf = y44vf;
z45vf = y45vf;
z46vf = y46vf;
z47vf = y47vf;
z48vf = y48vf;

z41vh = y41vh;
z42vh = y42vh;
z43vh = y43vh;
z44vh = y44vh;
z45vh = y45vh;
z46vh = y46vh;
z47vh = y47vh;
z48vh = y48vh;

% writing values obtained from mathematica === ends

% preparing spherical means for all 15 type of points for all three
% pressure, vel1 and vel2 at half and full time steps

t1phx1 = x11vh + x13vh + x15vh + x17vh;
t1phx2 = x21vh + x23vh + x25vh + x27vh;
t1phx3 = x31vh + x33vh + x35vh + x37vh;

t1pfx1 = x11vf + x13vf + x15vf + x17vf;
t1pfx2 = x21vf + x23vf + x25vf + x27vf;
t1pfx3 = x31vf + x33vf + x35vf + x37vf;

t1uhy1 = y11vh + y13vh + y15vh + y17vh;
t1uhy2 = y21vh + y23vh + y25vh + y27vh;
t1uhy3 = y31vh + y33vh + y35vh + y37vh;
t1uhy4 = y41vh + y43vh + y45vh + y47vh;

t1ufy1 = y11vf + y13vf + y15vf + y17vf;
t1ufy2 = y21vf + y23vf + y25vf + y27vf;
t1ufy3 = y31vf + y33vf + y35vf + y37vf;
t1ufy4 = y41vf + y43vf + y45vf + y47vf;

t1vhz1 = z11vh + z13vh + z15vh + z17vh;
t1vhz2 = z21vh + z23vh + z25vh + z27vh;
t1vhz3 = z31vh + z33vh + z35vh + z37vh;
t1vhz4 = z41vh + z43vh + z45vh + z47vh;

t1vfz1 = z11vf + z13vf + z15vf + z17vf;
t1vfz2 = z21vf + z23vf + z25vf + z27vf;
t1vfz3 = z31vf + z33vf + z35vf + z37vf;
t1vfz4 = z41vf + z43vf + z45vf + z47vf;


t2phx1 = x11vh + x13vh;
t2phx2 = x21vh + x23vh;
t2phx3 = x31vh + x33vh;

t2pfx1 = x11vf + x13vf;
t2pfx2 = x21vf + x23vf;
t2pfx3 = x31vf + x33vf;

t2uhy1 = y11vh + y13vh;
t2uhy2 = y21vh + y23vh;
t2uhy3 = y31vh + y33vh;
t2uhy4 = y41vh + y43vh;

t2ufy1 = y11vf + y13vf;
t2ufy2 = y21vf + y23vf;
t2ufy3 = y31vf + y33vf;
t2ufy4 = y41vf + y43vf;

t2vhz1 = z11vh + z13vh;
t2vhz2 = z21vh + z23vh;
t2vhz3 = z31vh + z33vh;
t2vhz4 = z41vh + z43vh;

t2vfz1 = z11vf + z13vf;
t2vfz2 = z21vf + z23vf;
t2vfz3 = z31vf + z33vf;
t2vfz4 = z41vf + z43vf;



t3phx1 = x15vh + x17vh;
t3phx2 = x25vh + x27vh;
t3phx3 = x35vh + x37vh;

t3pfx1 = x15vf + x17vf;
t3pfx2 = x25vf + x27vf;
t3pfx3 = x35vf + x37vf;

t3uhy1 = y15vh + y17vh;
t3uhy2 = y25vh + y27vh;
t3uhy3 = y35vh + y37vh;
t3uhy4 = y45vh + y47vh;

t3ufy1 = y15vf + y17vf;
t3ufy2 = y25vf + y27vf;
t3ufy3 = y35vf + y37vf;
t3ufy4 = y45vf + y47vf;

t3vhz1 = z15vh + z17vh;
t3vhz2 = z25vh + z27vh;
t3vhz3 = z35vh + z37vh;
t3vhz4 = z45vh + z47vh;

t3vfz1 = z15vf + z17vf;
t3vfz2 = z25vf + z27vf;
t3vfz3 = z35vf + z37vf;
t3vfz4 = z45vf + z47vf;


t4phx1 = x11vh + x17vh;
t4phx2 = x21vh + x27vh;
t4phx3 = x31vh + x37vh;

t4pfx1 = x11vf + x17vf;
t4pfx2 = x21vf + x27vf;
t4pfx3 = x31vf + x37vf;

t4uhy1 = y11vh + y17vh;
t4uhy2 = y21vh + y27vh;
t4uhy3 = y31vh + y37vh;
t4uhy4 = y41vh + y47vh;

t4ufy1 = y11vf + y17vf;
t4ufy2 = y21vf + y27vf;
t4ufy3 = y31vf + y37vf;
t4ufy4 = y41vf + y47vf;

t4vhz1 = z11vh + z17vh;
t4vhz2 = z21vh + z27vh;
t4vhz3 = z31vh + z37vh;
t4vhz4 = z41vh + z47vh;

t4vfz1 = z11vf + z17vf;
t4vfz2 = z21vf + z27vf;
t4vfz3 = z31vf + z37vf;
t4vfz4 = z41vf + z47vf;



t5phx1 = x15vh + x13vh;
t5phx2 = x25vh + x23vh;
t5phx3 = x35vh + x33vh;

t5pfx1 = x15vf + x13vf;
t5pfx2 = x25vf + x23vf;
t5pfx3 = x35vf + x33vf;

t5uhy1 = y15vh + y13vh;
t5uhy2 = y25vh + y23vh;
t5uhy3 = y35vh + y33vh;
t5uhy4 = y45vh + y43vh;

t5ufy1 = y15vf + y13vf;
t5ufy2 = y25vf + y23vf;
t5ufy3 = y35vf + y33vf;
t5ufy4 = y45vf + y43vf;

t5vhz1 = z15vh + z13vh;
t5vhz2 = z25vh + z23vh;
t5vhz3 = z35vh + z33vh;
t5vhz4 = z45vh + z43vh;

t5vfz1 = z15vf + z13vf;
t5vfz2 = z25vf + z23vf;
t5vfz3 = z35vf + z33vf;
t5vfz4 = z45vf + z43vf;



t6phx1 = x11vh;
t6phx2 = x21vh;
t6phx3 = x31vh;

t6pfx1 = x11vf;
t6pfx2 = x21vf;
t6pfx3 = x31vf;

t6uhy1 = y11vh;
t6uhy2 = y21vh;
t6uhy3 = y31vh;
t6uhy4 = y41vh;

t6ufy1 = y11vf;
t6ufy2 = y21vf;
t6ufy3 = y31vf;
t6ufy4 = y41vf;

t6vhz1 = z11vh;
t6vhz2 = z21vh;
t6vhz3 = z31vh;
t6vhz4 = z41vh;

t6vfz1 = z11vf;
t6vfz2 = z21vf;
t6vfz3 = z31vf;
t6vfz4 = z41vf;


t7phx1 = x13vh;
t7phx2 = x23vh;
t7phx3 = x33vh;

t7pfx1 = x13vf;
t7pfx2 = x23vf;
t7pfx3 = x33vf;

t7uhy1 = y13vh;
t7uhy2 = y23vh;
t7uhy3 = y33vh;
t7uhy4 = y43vh;

t7ufy1 = y13vf;
t7ufy2 = y23vf;
t7ufy3 = y33vf;
t7ufy4 = y43vf;

t7vhz1 = z13vh;
t7vhz2 = z23vh;
t7vhz3 = z33vh;
t7vhz4 = z43vh;

t7vfz1 = z13vf;
t7vfz2 = z23vf;
t7vfz3 = z33vf;
t7vfz4 = z43vf;



t8phx1 = x17vh;
t8phx2 = x27vh;
t8phx3 = x37vh;

t8pfx1 = x17vf;
t8pfx2 = x27vf;
t8pfx3 = x37vf;

t8uhy1 = y17vh;
t8uhy2 = y27vh;
t8uhy3 = y37vh;
t8uhy4 = y47vh;

t8ufy1 = y17vf;
t8ufy2 = y27vf;
t8ufy3 = y37vf;
t8ufy4 = y47vf;

t8vhz1 = z17vh;
t8vhz2 = z27vh;
t8vhz3 = z37vh;
t8vhz4 = z47vh;

t8vfz1 = z17vf;
t8vfz2 = z27vf;
t8vfz3 = z37vf;
t8vfz4 = z47vf;


t9phx1 = x15vh;
t9phx2 = x25vh;
t9phx3 = x35vh;

t9pfx1 = x15vf;
t9pfx2 = x25vf;
t9pfx3 = x35vf;

t9uhy1 = y15vh;
t9uhy2 = y25vh;
t9uhy3 = y35vh;
t9uhy4 = y45vh;

t9ufy1 = y15vf;
t9ufy2 = y25vf;
t9ufy3 = y35vf;
t9ufy4 = y45vf;

t9vhz1 = z15vh;
t9vhz2 = z25vh;
t9vhz3 = z35vh;
t9vhz4 = z45vh;


t9vfz1 = z15vf;
t9vfz2 = z25vf;
t9vfz3 = z35vf;
t9vfz4 = z45vf;




t10phx1 = x12vh + x16vh;
t10phx2 = x22vh + x26vh;
t10phx3 = x32vh + x36vh;

t10pfx1 = x12vf + x16vf;
t10pfx2 = x22vf + x26vf;
t10pfx3 = x32vf + x36vf;

t10uhy1 = y12vh + y16vh;
t10uhy2 = y22vh + y26vh;
t10uhy3 = y32vh + y36vh;
t10uhy4 = y42vh + y46vh;

t10ufy1 = y12vf + y16vf;
t10ufy2 = y22vf + y26vf;
t10ufy3 = y32vf + y36vf;
t10ufy4 = y42vf + y46vf;


t10vhz1 = z12vh + z16vh;
t10vhz2 = z22vh + z26vh;
t10vhz3 = z32vh + z36vh;
t10vhz4 = z42vh + z46vh;


t10vfz1 = z12vf + z16vf;
t10vfz2 = z22vf + z26vf;
t10vfz3 = z32vf + z36vf;
t10vfz4 = z42vf + z46vf;



t11phx1 = x12vh;
t11phx2 = x22vh;
t11phx3 = x32vh;

t11pfx1 = x12vf;
t11pfx2 = x22vf;
t11pfx3 = x32vf;

t11uhy1 = y12vh;
t11uhy2 = y22vh;
t11uhy3 = y32vh;
t11uhy4 = y42vh;


t11ufy1 = y12vf;
t11ufy2 = y22vf;
t11ufy3 = y32vf;
t11ufy4 = y42vf;



t11vhz1 = z12vh;
t11vhz2 = z22vh;
t11vhz3 = z32vh;
t11vhz4 = z42vh;


t11vfz1 = z12vf;
t11vfz2 = z22vf;
t11vfz3 = z32vf;
t11vfz4 = z42vf;


t12phx1 = x16vh;
t12phx2 = x26vh;
t12phx3 = x36vh;

t12pfx1 = x16vf;
t12pfx2 = x26vf;
t12pfx3 = x36vf;

t12uhy1 = y16vh;
t12uhy2 = y26vh;
t12uhy3 = y36vh;
t12uhy4 = y46vh;


t12ufy1 = y16vf;
t12ufy2 = y26vf;
t12ufy3 = y36vf;
t12ufy4 = y46vf;


t12vhz1 = z16vh;
t12vhz2 = z26vh;
t12vhz3 = z36vh;
t12vhz4 = z46vh;


t12vfz1 = z16vf;
t12vfz2 = z26vf;
t12vfz3 = z36vf;
t12vfz4 = z46vf;



t13phx1 = x14vh + x18vh;
t13phx2 = x24vh + x28vh;
t13phx3 = x34vh + x38vh;

t13pfx1 = x14vf + x18vf;
t13pfx2 = x24vf + x28vf;
t13pfx3 = x34vf + x38vf;

t13uhy1 = y14vh + y18vh;
t13uhy2 = y24vh + y28vh;
t13uhy3 = y34vh + y38vh;
t13uhy4 = y44vh + y48vh;



t13ufy1 = y14vf + y18vf;
t13ufy2 = y24vf + y28vf;
t13ufy3 = y34vf + y38vf;
t13ufy4 = y44vf + y48vf;

t13vhz1 = z14vh + z18vh;
t13vhz2 = z24vh + z28vh;
t13vhz3 = z34vh + z38vh;
t13vhz4 = z44vh + z48vh;

t13vfz1 = z14vf + z18vf;
t13vfz2 = z24vf + z28vf;
t13vfz3 = z34vf + z38vf;
t13vfz4 = z44vf + z48vf;


t14phx1 = x18vh;
t14phx2 = x28vh;
t14phx3 = x38vh;

t14pfx1 = x18vf;
t14pfx2 = x28vf;
t14pfx3 = x38vf;

t14uhy1 = y18vh;
t14uhy2 = y28vh;
t14uhy3 = y38vh;
t14uhy4 = y48vh;


t14ufy1 = y18vf;
t14ufy2 = y28vf;
t14ufy3 = y38vf;
t14ufy4 = y48vf;


t14vhz1 = z18vh;
t14vhz2 = z28vh;
t14vhz3 = z38vh;
t14vhz4 = z48vh;


t14vfz1 = z18vf;
t14vfz2 = z28vf;
t14vfz3 = z38vf;
t14vfz4 = z48vf;



t15phx1 = x14vh;
t15phx2 = x24vh;
t15phx3 = x34vh;

t15pfx1 = x14vf;
t15pfx2 = x24vf;
t15pfx3 = x34vf;

t15uhy1 = y14vh;
t15uhy2 = y24vh;
t15uhy3 = y34vh;
t15uhy4 = y44vh;


t15ufy1 = y14vf;
t15ufy2 = y24vf;
t15ufy3 = y34vf;
t15ufy4 = y44vf;

t15vhz1 = z14vh;
t15vhz2 = z24vh;
t15vhz3 = z34vh;
t15vhz4 = z44vh;

t15vfz1 = z14vf;
t15vfz2 = z24vf;
t15vfz3 = z34vf;
t15vfz4 = z44vf;

% spherical means ends 

% putting initial values in all 15 type of points %%%%%%%%%%%%%%%%%%%%%%

t1_press_iv = press1(2:end,2:end);
t1_vel1_iv = vel11(2:end,2:end);
t1_vel2_iv = vel21(2:end,2:end);

t2_press_iv = press1(1,2:end);
t2_vel1_iv = vel11(1,2:end);
t2_vel2_iv = vel21(1,2:end);

t3_press_iv = press7(end,2:end);
t3_vel1_iv = vel17(end,2:end);
t3_vel2_iv = vel27(end,2:end);

t4_press_iv = press1(2:end,1);
t4_vel1_iv = vel11(2:end,1);
t4_vel2_iv = vel21(2:end,1);

t5_press_iv = press3(2:end,end);
t5_vel1_iv = vel13(2:end,end);
t5_vel2_iv = vel23(2:end,end);

t6_press_iv = press1(1,1);
t6_vel1_iv = vel11(1,1);
t6_vel2_iv = vel21(1,1);

t7_press_iv = press3(1,end);
t7_vel1_iv = vel13(1,end);
t7_vel2_iv = vel23(1,end);

t8_press_iv = press7(end,1);
t8_vel1_iv = vel17(end,1);
t8_vel2_iv = vel27(end,1);

t9_press_iv = press5(end,end);
t9_vel1_iv = vel15(end,end);
t9_vel2_iv = vel25(end,end);

t10_press_iv = press2(2:end,:);
t10_vel1_iv = vel12(2:end,:);
t10_vel2_iv = vel22(2:end,:);

t11_press_iv = press2(1,:);
t11_vel1_iv = vel12(1,:);
t11_vel2_iv = vel22(1,:);

t12_press_iv = press6(end,:);
t12_vel1_iv = vel16(end,:);
t12_vel2_iv = vel26(end,:);

t13_press_iv = press4(:,1:end-1);
t13_vel1_iv = vel14(:,1:end-1);
t13_vel2_iv = vel24(:,1:end-1);

t14_press_iv = press8(:,1);
t14_vel1_iv = vel18(:,1);
t14_vel2_iv = vel28(:,1);

t15_press_iv = press4(:,end);
t15_vel1_iv = vel14(:,end);
t15_vel2_iv = vel24(:,end);

% putting initial values in all 15 type of points ends %%%%%%%%%%%%%%%

% loop starts
for ii = t:dt:tend

    % for type 1 point
    t1_press_h = (t1_press_iv * t1phx1) + 2*(t1_vel1_iv * t1phx2)/dt - 2*(t1_vel2_iv * t1phx3)/dt;
    t1_vel1_h = ( t1_vel1_iv ) + 2*(t1_press_iv * t1uhy1)/dt + (t1_vel1_iv * t1uhy2) - (t1_vel2_iv * t1uhy3) - (t1_vel1_iv * t1uhy4);
    t1_vel2_h = ( t1_vel2_iv ) - 2*(t1_press_iv * t1vhz1)/dt - (t1_vel1_iv * t1vhz2) + (t1_vel2_iv * t1vhz3) - (t1_vel2_iv * t1vhz4);
    
    t1_press_f = (t1_press_iv * t1pfx1) + (t1_vel1_iv * t1pfx2)/dt - (t1_vel2_iv * t1pfx3)/dt;
    t1_vel1_f = ( t1_vel1_iv ) + (t1_press_iv * t1ufy1)/dt + (t1_vel1_iv * t1ufy2) - (t1_vel2_iv * t1ufy3) - (t1_vel1_iv * t1ufy4);
    t1_vel2_f = ( t1_vel2_iv ) - (t1_press_iv * t1vfz1)/dt - (t1_vel1_iv * t1vfz2) + (t1_vel2_iv * t1vfz3) - (t1_vel2_iv * t1vfz4);
    
    
    % for type 2 point
    t2_press_h = (t2_press_iv * t2phx1) + 2*(t2_vel1_iv * t2phx2)/dt - 2*(t2_vel2_iv * t2phx3)/dt;
    t2_vel1_h = ( t2_vel1_iv ) + 2*(t2_press_iv * t2uhy1)/dt + (t2_vel1_iv * t2uhy2) - (t2_vel2_iv * t2uhy3) - (t2_vel1_iv * t2uhy4);
    t2_vel2_h = ( t2_vel2_iv ) - 2*(t2_press_iv * t2vhz1)/dt - (t2_vel1_iv * t2vhz2) + (t2_vel2_iv * t2vhz3) - (t2_vel2_iv * t2vhz4);
    
    t2_press_f = (t2_press_iv * t2pfx1) + (t2_vel1_iv * t2pfx2)/dt - (t2_vel2_iv * t2pfx3)/dt;
    t2_vel1_f = ( t2_vel1_iv ) + (t2_press_iv * t2ufy1)/dt + (t2_vel1_iv * t2ufy2) - (t2_vel2_iv * t2ufy3) - (t2_vel1_iv * t2ufy4);
    t2_vel2_f = ( t2_vel2_iv ) - (t2_press_iv * t2vfz1)/dt - (t2_vel1_iv * t2vfz2) + (t2_vel2_iv * t2vfz3) - (t2_vel2_iv * t2vfz4);
    
    
    % for type 3 point
    t3_press_h = (t3_press_iv * t3phx1) + 2*(t3_vel1_iv * t3phx2)/dt - 2*(t3_vel2_iv * t3phx3)/dt;
    t3_vel1_h = ( t3_vel1_iv ) + 2*(t3_press_iv * t3uhy1)/dt + (t3_vel1_iv * t3uhy2) - (t3_vel2_iv * t3uhy3) - (t3_vel1_iv * t3uhy4);
    t3_vel2_h = ( t3_vel2_iv ) - 2*(t3_press_iv * t3vhz1)/dt - (t3_vel1_iv * t3vhz2) + (t3_vel2_iv * t3vhz3) - (t3_vel2_iv * t3vhz4);
    
    t3_press_f = (t3_press_iv * t3pfx1) + (t3_vel1_iv * t3pfx2)/dt - (t3_vel2_iv * t3pfx3)/dt;
    t3_vel1_f = ( t3_vel1_iv ) + (t3_press_iv * t3ufy1)/dt + (t3_vel1_iv * t3ufy2) - (t3_vel2_iv * t3ufy3) - (t3_vel1_iv * t3ufy4);
    t3_vel2_f = ( t3_vel2_iv ) - (t3_press_iv * t3vfz1)/dt - (t3_vel1_iv * t3vfz2) + (t3_vel2_iv * t3vfz3) - (t3_vel2_iv * t3vfz4);
    
    
    % for type 4 point
    t4_press_h = (t4_press_iv * t4phx1) + 2*(t4_vel1_iv * t4phx2)/dt - 2*(t4_vel2_iv * t4phx3)/dt;
    t4_vel1_h = ( t4_vel1_iv ) + 2*(t4_press_iv * t4uhy1)/dt + (t4_vel1_iv * t4uhy2) - (t4_vel2_iv * t4uhy3) - (t4_vel1_iv * t4uhy4);
    t4_vel2_h = ( t4_vel2_iv ) - 2*(t4_press_iv * t4vhz1)/dt - (t4_vel1_iv * t4vhz2) + (t4_vel2_iv * t4vhz3) - (t4_vel2_iv * t4vhz4);
    
    t4_press_f = (t4_press_iv * t4pfx1) + (t4_vel1_iv * t4pfx2)/dt - (t4_vel2_iv * t4pfx3)/dt;
    t4_vel1_f = ( t4_vel1_iv ) + (t4_press_iv * t4ufy1)/dt + (t4_vel1_iv * t4ufy2) - (t4_vel2_iv * t4ufy3) - (t4_vel1_iv * t4ufy4);
    t4_vel2_f = ( t4_vel2_iv ) - (t4_press_iv * t4vfz1)/dt - (t4_vel1_iv * t4vfz2) + (t4_vel2_iv * t4vfz3) - (t4_vel2_iv * t4vfz4);
    
    
    % for type 5 point
    t5_press_h = (t5_press_iv * t5phx1) + 2*(t5_vel1_iv * t5phx2)/dt - 2*(t5_vel2_iv * t5phx3)/dt;
    t5_vel1_h = ( t5_vel1_iv ) + 2*(t5_press_iv * t5uhy1)/dt + (t5_vel1_iv * t5uhy2) - (t5_vel2_iv * t5uhy3) - (t5_vel1_iv * t5uhy4);
    t5_vel2_h = ( t5_vel2_iv ) - 2*(t5_press_iv * t5vhz1)/dt - (t5_vel1_iv * t5vhz2) + (t5_vel2_iv * t5vhz3) - (t5_vel2_iv * t5vhz4);
    
    t5_press_f = (t5_press_iv * t5pfx1) + (t5_vel1_iv * t5pfx2)/dt - (t5_vel2_iv * t5pfx3)/dt;
    t5_vel1_f = ( t5_vel1_iv ) + (t5_press_iv * t5ufy1)/dt + (t5_vel1_iv * t5ufy2) - (t5_vel2_iv * t5ufy3) - (t5_vel1_iv * t5ufy4);
    t5_vel2_f = ( t5_vel2_iv ) - (t5_press_iv * t5vfz1)/dt - (t5_vel1_iv * t5vfz2) + (t5_vel2_iv * t5vfz3) - (t5_vel2_iv * t5vfz4);
    
    
    % for type 6 point
    t6_press_h = (t6_press_iv * t6phx1) + 2*(t6_vel1_iv * t6phx2)/dt - 2*(t6_vel2_iv * t6phx3)/dt;
    t6_vel1_h = ( t6_vel1_iv ) + 2*(t6_press_iv * t6uhy1)/dt + (t6_vel1_iv * t6uhy2) - (t6_vel2_iv * t6uhy3) - (t6_vel1_iv * t6uhy4);
    t6_vel2_h = ( t6_vel2_iv ) - 2*(t6_press_iv * t6vhz1)/dt - (t6_vel1_iv * t6vhz2) + (t6_vel2_iv * t6vhz3) - (t6_vel2_iv * t6vhz4);
    
    t6_press_f = (t6_press_iv * t6pfx1) + (t6_vel1_iv * t6pfx2)/dt - (t6_vel2_iv * t6pfx3)/dt;
    t6_vel1_f = ( t6_vel1_iv ) + (t6_press_iv * t6ufy1)/dt + (t6_vel1_iv * t6ufy2) - (t6_vel2_iv * t6ufy3) - (t6_vel1_iv * t6ufy4);
    t6_vel2_f = ( t6_vel2_iv ) - (t6_press_iv * t6vfz1)/dt - (t6_vel1_iv * t6vfz2) + (t6_vel2_iv * t6vfz3) - (t6_vel2_iv * t6vfz4);
    
    
    % for type 7 point
    t7_press_h = (t7_press_iv * t7phx1) + 2*(t7_vel1_iv * t7phx2)/dt - 2*(t7_vel2_iv * t7phx3)/dt;
    t7_vel1_h = ( t7_vel1_iv ) + 2*(t7_press_iv * t7uhy1)/dt + (t7_vel1_iv * t7uhy2) - (t7_vel2_iv * t7uhy3) - (t7_vel1_iv * t7uhy4);
    t7_vel2_h = ( t7_vel2_iv ) - 2*(t7_press_iv * t7vhz1)/dt - (t7_vel1_iv * t7vhz2) + (t7_vel2_iv * t7vhz3) - (t7_vel2_iv * t7vhz4);
    
    t7_press_f = (t7_press_iv * t7pfx1) + (t7_vel1_iv * t7pfx2)/dt - (t7_vel2_iv * t7pfx3)/dt;
    t7_vel1_f = ( t7_vel1_iv ) + (t7_press_iv * t7ufy1)/dt + (t7_vel1_iv * t7ufy2) - (t7_vel2_iv * t7ufy3) - (t7_vel1_iv * t7ufy4);
    t7_vel2_f = ( t7_vel2_iv ) - (t7_press_iv * t7vfz1)/dt - (t7_vel1_iv * t7vfz2) + (t7_vel2_iv * t7vfz3) - (t7_vel2_iv * t7vfz4);
    
    
    % for type 8 point
    t8_press_h = (t8_press_iv * t8phx1) + 2*(t8_vel1_iv * t8phx2)/dt - 2*(t8_vel2_iv * t8phx3)/dt;
    t8_vel1_h = ( t8_vel1_iv ) + 2*(t8_press_iv * t8uhy1)/dt + (t8_vel1_iv * t8uhy2) - (t8_vel2_iv * t8uhy3) - (t8_vel1_iv * t8uhy4);
    t8_vel2_h = ( t8_vel2_iv ) - 2*(t8_press_iv * t8vhz1)/dt - (t8_vel1_iv * t8vhz2) + (t8_vel2_iv * t8vhz3) - (t8_vel2_iv * t8vhz4);
    
    t8_press_f = (t8_press_iv * t8pfx1) + (t8_vel1_iv * t8pfx2)/dt - (t8_vel2_iv * t8pfx3)/dt;
    t8_vel1_f = ( t8_vel1_iv ) + (t8_press_iv * t8ufy1)/dt + (t8_vel1_iv * t8ufy2) - (t8_vel2_iv * t8ufy3) - (t8_vel1_iv * t8ufy4);
    t8_vel2_f = ( t8_vel2_iv ) - (t8_press_iv * t8vfz1)/dt - (t8_vel1_iv * t1vfz2) + (t8_vel2_iv * t8vfz3) - (t8_vel2_iv * t8vfz4);
    
    
    % for type 9 point
    t9_press_h = (t9_press_iv * t9phx1) + 2*(t9_vel1_iv * t9phx2)/dt - 2*(t9_vel2_iv * t9phx3)/dt;
    t9_vel1_h = ( t9_vel1_iv ) + 2*(t9_press_iv * t9uhy1)/dt + (t9_vel1_iv * t9uhy2) - (t9_vel2_iv * t9uhy3) - (t9_vel1_iv * t9uhy4);
    t9_vel2_h = ( t9_vel2_iv ) - 2*(t9_press_iv * t9vhz1)/dt - (t9_vel1_iv * t9vhz2) + (t9_vel2_iv * t9vhz3) - (t9_vel2_iv * t9vhz4);
    
    t9_press_f = (t9_press_iv * t9pfx1) + (t9_vel1_iv * t9pfx2)/dt - (t9_vel2_iv * t9pfx3)/dt;
    t9_vel1_f = ( t9_vel1_iv ) + (t9_press_iv * t9ufy1)/dt + (t9_vel1_iv * t9ufy2) - (t9_vel2_iv * t9ufy3) - (t9_vel1_iv * t9ufy4);
    t9_vel2_f = ( t9_vel2_iv ) - (t9_press_iv * t9vfz1)/dt - (t9_vel1_iv * t9vfz2) + (t9_vel2_iv * t9vfz3) - (t9_vel2_iv * t9vfz4);
    
    
    % for type 10 point
    t10_press_h = (t10_press_iv * t10phx1) + 2*(t10_vel1_iv * t10phx2)/dt - 2*(t10_vel2_iv * t10phx3)/dt;
    t10_vel1_h = ( t10_vel1_iv ) + 2*(t10_press_iv * t10uhy1)/dt + (t10_vel1_iv * t10uhy2) - (t10_vel2_iv * t10uhy3) - (t10_vel1_iv * t10uhy4);
    t10_vel2_h = ( t10_vel2_iv ) - 2*(t10_press_iv * t10vhz1)/dt - (t10_vel1_iv * t10vhz2) + (t10_vel2_iv * t10vhz3) - (t10_vel2_iv * t10vhz4);
    
    t10_press_f = (t10_press_iv * t10pfx1) + (t10_vel1_iv * t10pfx2)/dt - (t10_vel2_iv * t10pfx3)/dt;
    t10_vel1_f = ( t10_vel1_iv ) + (t10_press_iv * t10ufy1)/dt + (t10_vel1_iv * t10ufy2) - (t10_vel2_iv * t10ufy3) - (t10_vel1_iv * t10ufy4);
    t10_vel2_f = ( t10_vel2_iv ) - (t10_press_iv * t10vfz1)/dt - (t10_vel1_iv * t10vfz2) + (t10_vel2_iv * t10vfz3) - (t10_vel2_iv * t10vfz4);
    
    
    % for type 11 point
    t11_press_h = (t11_press_iv * t11phx1) + 2*(t11_vel1_iv * t11phx2)/dt - 2*(t11_vel2_iv * t11phx3)/dt;
    t11_vel1_h = ( t11_vel1_iv ) + 2*(t11_press_iv * t11uhy1)/dt + (t11_vel1_iv * t11uhy2) - (t11_vel2_iv * t11uhy3) - (t11_vel1_iv * t11uhy4);
    t11_vel2_h = ( t11_vel2_iv ) - 2*(t11_press_iv * t11vhz1)/dt - (t11_vel1_iv * t11vhz2) + (t11_vel2_iv * t11vhz3) - (t11_vel2_iv * t11vhz4);
    
    t11_press_f = (t11_press_iv * t11pfx1) + (t11_vel1_iv * t11pfx2)/dt - (t11_vel2_iv * t11pfx3)/dt;
    t11_vel1_f = ( t11_vel1_iv ) + (t11_press_iv * t11ufy1)/dt + (t11_vel1_iv * t11ufy2) - (t11_vel2_iv * t11ufy3) - (t11_vel1_iv * t11ufy4);
    t11_vel2_f = ( t11_vel2_iv ) - (t11_press_iv * t11vfz1)/dt - (t11_vel1_iv * t11vfz2) + (t11_vel2_iv * t11vfz3) - (t11_vel2_iv * t11vfz4);
    
    
    % for type 12 point
    t12_press_h = (t12_press_iv * t12phx1) + 2*(t12_vel1_iv * t12phx2)/dt - 2*(t12_vel2_iv * t12phx3)/dt;
    t12_vel1_h = ( t12_vel1_iv ) + 2*(t12_press_iv * t12uhy1)/dt + (t12_vel1_iv * t12uhy2) - (t12_vel2_iv * t12uhy3) - (t12_vel1_iv * t12uhy4);
    t12_vel2_h = ( t12_vel2_iv ) - 2*(t12_press_iv * t12vhz1)/dt - (t12_vel1_iv * t12vhz2) + (t12_vel2_iv * t12vhz3) - (t12_vel2_iv * t12vhz4);
    
    t12_press_f = (t12_press_iv * t12pfx1) + (t12_vel1_iv * t12pfx2)/dt - (t12_vel2_iv * t12pfx3)/dt;
    t12_vel1_f = ( t12_vel1_iv ) + (t12_press_iv * t12ufy1)/dt + (t12_vel1_iv * t12ufy2) - (t12_vel2_iv * t12ufy3) - (t12_vel1_iv * t12ufy4);
    t12_vel2_f = ( t12_vel2_iv ) - (t12_press_iv * t12vfz1)/dt - (t12_vel1_iv * t12vfz2) + (t12_vel2_iv * t12vfz3) - (t12_vel2_iv * t12vfz4);
    
    
    % for type 13 point
    t13_press_h = (t13_press_iv * t13phx1) + 2*(t13_vel1_iv * t13phx2)/dt - 2*(t13_vel2_iv * t13phx3)/dt;
    t13_vel1_h = ( t13_vel1_iv ) + 2*(t13_press_iv * t13uhy1)/dt + (t13_vel1_iv * t13uhy2) - (t13_vel2_iv * t13uhy3) - (t13_vel1_iv * t13uhy4);
    t13_vel2_h = ( t13_vel2_iv ) - 2*(t13_press_iv * t13vhz1)/dt - (t13_vel1_iv * t13vhz2) + (t13_vel2_iv * t13vhz3) - (t13_vel2_iv * t13vhz4);
    
    t13_press_f = (t13_press_iv * t13pfx1) + (t13_vel1_iv * t13pfx2)/dt - (t13_vel2_iv * t13pfx3)/dt;
    t13_vel1_f = ( t13_vel1_iv ) + (t13_press_iv * t13ufy1)/dt + (t13_vel1_iv * t13ufy2) - (t13_vel2_iv * t13ufy3) - (t13_vel1_iv * t13ufy4);
    t13_vel2_f = ( t13_vel2_iv ) - (t13_press_iv * t13vfz1)/dt - (t13_vel1_iv * t13vfz2) + (t13_vel2_iv * t13vfz3) - (t13_vel2_iv * t13vfz4);
    
    
    % for type 14 point
    t14_press_h = (t14_press_iv * t14phx1) + 2*(t14_vel1_iv * t14phx2)/dt - 2*(t14_vel2_iv * t14phx3)/dt;
    t14_vel1_h = ( t14_vel1_iv ) + 2*(t14_press_iv * t14uhy1)/dt + (t14_vel1_iv * t14uhy2) - (t14_vel2_iv * t14uhy3) - (t14_vel1_iv * t14uhy4);
    t14_vel2_h = ( t14_vel2_iv ) - 2*(t14_press_iv * t14vhz1)/dt - (t14_vel1_iv * t14vhz2) + (t14_vel2_iv * t14vhz3) - (t14_vel2_iv * t14vhz4);
    
    t14_press_f = (t14_press_iv * t14pfx1) + (t14_vel1_iv * t14pfx2)/dt - (t14_vel2_iv * t1pfx3)/dt;
    t14_vel1_f = ( t14_vel1_iv ) + (t14_press_iv * t14ufy1)/dt + (t14_vel1_iv * t14ufy2) - (t14_vel2_iv * t14ufy3) - (t14_vel1_iv * t14ufy4);
    t14_vel2_f = ( t14_vel2_iv ) - (t14_press_iv * t14vfz1)/dt - (t14_vel1_iv * t14vfz2) + (t14_vel2_iv * t14vfz3) - (t14_vel2_iv * t14vfz4);
    
    
    % for type 15 point
    t15_press_h = (t15_press_iv * t15phx1) + 2*(t15_vel1_iv * t15phx2)/dt - 2*(t15_vel2_iv * t15phx3)/dt;
    t15_vel1_h = ( t15_vel1_iv ) + 2*(t15_press_iv * t15uhy1)/dt + (t15_vel1_iv * t15uhy2) - (t15_vel2_iv * t15uhy3) - (t15_vel1_iv * t15uhy4);
    t15_vel2_h = ( t15_vel2_iv ) - 2*(t15_press_iv * t15vhz1)/dt - (t15_vel1_iv * t15vhz2) + (t15_vel2_iv * t15vhz3) - (t15_vel2_iv * t15vhz4);
    
    t15_press_f = (t15_press_iv * t15pfx1) + (t15_vel1_iv * t15pfx2)/dt - (t15_vel2_iv * t15pfx3)/dt;
    t15_vel1_f = ( t15_vel1_iv ) + (t15_press_iv * t15ufy1)/dt + (t15_vel1_iv * t15ufy2) - (t15_vel2_iv * t15ufy3) - (t15_vel1_iv * t15ufy4);
    t15_vel2_f = ( t15_vel2_iv ) - (t15_press_iv * t15vfz1)/dt - (t15_vel1_iv * t15vfz2) + (t15_vel2_iv * t15vfz3) - (t15_vel2_iv * t15vfz4);
    
    % calculating for 8 points on a cell
    
    % for pressure
    
    % point 1
    press1_1(2:end,2:end) = t1_press_h;
    press1_1(1,2:end) = t2_press_h;
    press1_1(2:end,1) = t4_press_h;
    press1_1(1,1) = t6_press_h;
    
    press1_2(2:end,2:end) = t1_press_f;
    press1_2(1,2:end) = t2_press_f;
    press1_2(2:end,1) = t4_press_f;
    press1_2(1,1) = t6_press_f;
    
    % point 2
    
    press2_1(2:end,:) = t10_press_h;
    press2_1(1,:) = t11_press_h;
    
    press2_2(2:end,:) = t10_press_f;
    press2_2(1,:) = t11_press_f;
    
    % point 3
    
    press3_1(2:end,1:end-1) = t1_press_h;
    press3_1(1,1:end-1) = t2_press_h;
    press3_1(2:end,end) = t5_press_h;
    press3_1(1,end) = t7_press_h;
    
    press3_2(2:end,1:end-1) = t1_press_f;
    press3_2(1,1:end-1) = t2_press_f;
    press3_2(2:end,end) = t5_press_f;
    press3_2(1,end) = t7_press_f;
    
    % point 4
    
    press4_1(:,1:end-1) = t13_press_h;
    press4_1(:,end) = t15_press_h;
    
    press4_2(:,1:end-1) = t13_press_f;
    press4_2(:,end) = t15_press_f;
    
    % point 5
    
    press5_1(1:end-1,1:end-1) = t1_press_h;
    press5_1(end,1:end-1) = t3_press_h;
    press5_1(1:end-1,end) = t5_press_h;
    press5_1(end,end) = t9_press_h;
    
    press5_2(1:end-1,1:end-1) = t1_press_f;
    press5_2(end,1:end-1) = t3_press_f;
    press5_2(1:end-1,end) = t5_press_f;
    press5_2(end,end) = t9_press_f;
    
    % point 6
    
    press6_1(1:end-1,:) = t10_press_h;
    press6_1(end,:) = t12_press_h;
    
    press6_2(1:end-1,:) = t10_press_f;
    press6_2(end,:) = t12_press_f;
    
    % point 7
    
    press7_1(1:end-1,2:end) = t1_press_h;
    press7_1(1:end-1,1) = t4_press_h;
    press7_1(end,2:end) = t3_press_h;
    press7_1(end,1) = t8_press_h;
    
    press7_2(1:end-1,2:end) = t1_press_f;
    press7_2(1:end-1,1) = t4_press_f;
    press7_2(end,2:end) = t3_press_f;
    press7_2(end,1) = t8_press_f;
    
    % point 8
    
    press8_1(:,2:end) = t13_press_h;
    press8_1(:,1) = t14_press_h;
    
    press8_2(:,2:end) = t13_press_f;
    press8_2(:,1) = t14_press_f;
    
    
    
    % for velocity 1
    
    % point 1
    vel11_1(2:end,2:end) = t1_vel1_h;
    vel11_1(1,2:end) = t2_vel1_h;
    vel11_1(2:end,1) = t4_vel1_h;
    vel11_1(1,1) = t6_vel1_h;
    
    vel11_2(2:end,2:end) = t1_vel1_f;
    vel11_2(1,2:end) = t2_vel1_f;
    vel11_2(2:end,1) = t4_vel1_f;
    vel11_2(1,1) = t6_vel1_f;
    
    % point 2
    
    vel12_1(2:end,:) = t10_vel1_h;
    vel12_1(1,:) = t11_vel1_h;
    
    vel12_2(2:end,:) = t10_vel1_f;
    vel12_2(1,:) = t11_vel1_f;
    
    % point 3
    
    vel13_1(2:end,1:end-1) = t1_vel1_h;
    vel13_1(1,1:end-1) = t2_vel1_h;
    vel13_1(2:end,end) = t5_vel1_h;
    vel13_1(1,end) = t7_vel1_h;
    
    vel13_2(2:end,1:end-1) = t1_vel1_f;
    vel13_2(1,1:end-1) = t2_vel1_f;
    vel13_2(2:end,end) = t5_vel1_f;
    vel13_2(1,end) = t7_vel1_f;
    
    % point 4
    
    vel14_1(:,1:end-1) = t13_vel1_h;
    vel14_1(:,end) = t15_vel1_h;
    
    vel14_2(:,1:end-1) = t13_vel1_f;
    vel14_2(:,end) = t15_vel1_f;
    
    % point 5
    
    vel15_1(1:end-1,1:end-1) = t1_vel1_h;
    vel15_1(end,1:end-1) = t3_vel1_h;
    vel15_1(1:end-1,end) = t5_vel1_h;
    vel15_1(end,end) = t9_vel1_h;
    
    vel15_2(1:end-1,1:end-1) = t1_vel1_f;
    vel15_2(end,1:end-1) = t3_vel1_f;
    vel15_2(1:end-1,end) = t5_vel1_f;
    vel15_2(end,end) = t9_vel1_f;
    
    % point 6
    
    vel16_1(1:end-1,:) = t10_vel1_h;
    vel16_1(end,:) = t12_vel1_h;
    
    vel16_2(1:end-1,:) = t10_vel1_f;
    vel16_2(end,:) = t12_vel1_f;
    
    % point 7
    
    vel17_1(1:end-1,2:end) = t1_vel1_h;
    vel17_1(1:end-1,1) = t4_vel1_h;
    vel17_1(end,2:end) = t3_vel1_h;
    vel17_1(end,1) = t8_vel1_h;
    
    vel17_2(1:end-1,2:end) = t1_vel1_f;
    vel17_2(1:end-1,1) = t4_vel1_f;
    vel17_2(end,2:end) = t3_vel1_f;
    vel17_2(end,1) = t8_vel1_f;
    
    % point 8
    
    vel18_1(:,2:end) = t13_vel1_h;
    vel18_1(:,1) = t14_vel1_h;
    
    vel18_2(:,2:end) = t13_vel1_f;
    vel18_2(:,1) = t14_vel1_f;
    
    
    % for velocity 2
    
    % point 1
    vel21_1(2:end,2:end) = t1_vel2_h;
    vel21_1(1,2:end) = t2_vel2_h;
    vel21_1(2:end,1) = t4_vel2_h;
    vel21_1(1,1) = t6_vel2_h;
    
    vel21_2(2:end,2:end) = t1_vel2_f;
    vel21_2(1,2:end) = t2_vel2_f;
    vel21_2(2:end,1) = t4_vel2_f;
    vel21_2(1,1) = t6_vel2_f;
    
    % point 2
    
    vel22_1(2:end,:) = t10_vel2_h;
    vel22_1(1,:) = t11_vel2_h;
    
    vel22_2(2:end,:) = t10_vel2_f;
    vel22_2(1,:) = t11_vel2_f;
    
    % point 3
    
    vel23_1(2:end,1:end-1) = t1_vel2_h;
    vel23_1(1,1:end-1) = t2_vel2_h;
    vel23_1(2:end,end) = t5_vel2_h;
    vel23_1(1,end) = t7_vel2_h;
    
    vel23_2(2:end,1:end-1) = t1_vel2_f;
    vel23_2(1,1:end-1) = t2_vel2_f;
    vel23_2(2:end,end) = t5_vel2_f;
    vel23_2(1,end) = t7_vel2_f;
    
    % point 4
    
    vel24_1(:,1:end-1) = t13_vel2_h;
    vel24_1(:,end) = t15_vel2_h;
    
    vel24_2(:,1:end-1) = t13_vel2_f;
    vel24_2(:,end) = t15_vel2_f;
    
    % point 5
    
    vel25_1(1:end-1,1:end-1) = t1_vel2_h;
    vel25_1(end,1:end-1) = t3_vel2_h;
    vel25_1(1:end-1,end) = t5_vel2_h;
    vel25_1(end,end) = t9_vel2_h;
    
    vel25_2(1:end-1,1:end-1) = t1_vel2_f;
    vel25_2(end,1:end-1) = t3_vel2_f;
    vel25_2(1:end-1,end) = t5_vel2_f;
    vel25_2(end,end) = t9_vel2_f;
    
    % point 6
    
    vel26_1(1:end-1,:) = t10_vel2_h;
    vel26_1(end,:) = t12_vel2_h;
    
    vel26_2(1:end-1,:) = t10_vel2_f;
    vel26_2(end,:) = t12_vel2_f;
    
    % point 7
    
    vel27_1(1:end-1,2:end) = t1_vel2_h;
    vel27_1(1:end-1,1) = t4_vel2_h;
    vel27_1(end,2:end) = t3_vel2_h;
    vel27_1(end,1) = t8_vel2_h;
    
    vel27_2(1:end-1,2:end) = t1_vel2_f;
    vel27_2(1:end-1,1) = t4_vel2_f;
    vel27_2(end,2:end) = t3_vel2_f;
    vel27_2(end,1) = t8_vel2_f;
    
    % point 8
    
    vel28_1(:,2:end) = t13_vel2_h;
    vel28_1(:,1) = t14_vel2_h;
    
    vel28_2(:,2:end) = t13_vel2_f;
    vel28_2(:,1) = t14_vel2_f;
    
    % boundary conditions
    
    press1_1(:,1) = press3_1(:,end);
    press8_1(:,1) = press4_1(:,end);
    press7_1(:,1) = press5_1(:,end);
    press1_1(1,:) = press7_1(end,:);
    press2_1(1,:) = press6_1(end,:);
    press3_1(1,:) = press5_1(end,:);
    
    press1_2(:,1) = press3_2(:,end);
    press8_2(:,1) = press4_2(:,end);
    press7_2(:,1) = press5_2(:,end);
    press1_2(1,:) = press7_2(end,:);
    press2_2(1,:) = press6_2(end,:);
    press3_2(1,:) = press5_2(end,:);
    
    vel11_1(:,1) = vel13_1(:,end);
    vel18_1(:,1) = vel14_1(:,end);
    vel17_1(:,1) = vel15_1(:,end);
    vel11_1(1,:) = vel17_1(end,:);
    vel12_1(1,:) = vel16_1(end,:);
    vel13_1(1,:) = vel15_1(end,:);
    
    vel11_2(:,1) = vel13_2(:,end);
    vel18_2(:,1) = vel14_2(:,end);
    vel17_2(:,1) = vel15_2(:,end);
    vel11_2(1,:) = vel17_2(end,:);
    vel12_2(1,:) = vel16_2(end,:);
    vel13_2(1,:) = vel15_2(end,:);
    
    vel21_1(:,1) = vel23_1(:,end);
    vel28_1(:,1) = vel24_1(:,end);
    vel27_1(:,1) = vel25_1(:,end);
    vel21_1(1,:) = vel27_1(end,:);
    vel22_1(1,:) = vel26_1(end,:);
    vel23_1(1,:) = vel25_1(end,:);
    
    vel21_2(:,1) = vel23_2(:,end);
    vel28_2(:,1) = vel24_2(:,end);
    vel27_2(:,1) = vel25_2(:,end);
    vel21_2(1,:) = vel27_2(end,:);
    vel22_2(1,:) = vel26_2(end,:);
    vel23_2(1,:) = vel25_2(end,:);
    
    % bc ends
     
    % parameters for flux through face one

    pln1 = press1;
    pmn1 = press2;
    prn1 = press3;
    plnplushalf1 = press1_1; 
    pmnplushalf1 = press2_1;
    prnplushalf1 = press3_1;
    plnplusone1 = press1_2;
    pmnplusone1 = press2_2;
    prnplusone1 = press3_2;
    
    uln1 = vel11;
    umn1 = vel12;
    urn1 = vel13;
    ulnplushalf1 = vel11_1; 
    umnplushalf1 = vel12_1;
    urnplushalf1 = vel13_1;
    ulnplusone1 = vel11_2;
    umnplusone1 = vel12_2;
    urnplusone1 = vel13_2;
    
    vln1 = vel21;
    vmn1 = vel22;
    vrn1 = vel23;
    vlnplushalf1 = vel21_1; 
    vmnplushalf1 = vel22_1;
    vrnplushalf1 = vel23_1;
    vlnplusone1 = vel21_2;
    vmnplusone1 = vel22_2;
    vrnplusone1 = vel23_2;

    % parameters for flux through face two

    pln2 = press3;
    pmn2 = press4;
    prn2 = press5;
    plnplushalf2 = press3_1;
    pmnplushalf2 = press4_1;
    prnplushalf2 = press5_1;
    plnplusone2 = press3_2;
    pmnplusone2 = press4_2;
    prnplusone2 = press5_2;
    
    uln2 = vel13;
    umn2 = vel14;
    urn2 = vel15;
    ulnplushalf2 = vel13_1;
    umnplushalf2 = vel14_1;
    urnplushalf2 = vel15_1;
    ulnplusone2 = vel13_2;
    umnplusone2 = vel14_2;
    urnplusone2 = vel15_2;
    
    vln2 = vel23;
    vmn2 = vel24;
    vrn2 = vel25;
    vlnplushalf2 = vel23_1;
    vmnplushalf2 = vel24_1;
    vrnplushalf2 = vel25_1;
    vlnplusone2 = vel23_2;
    vmnplusone2 = vel24_2;
    vrnplusone2 = vel25_2;

    % parameters for flux through face three

    pln3 = press5;
    pmn3 = press6;
    prn3 = press7;
    plnplushalf3 = press5_1; 
    pmnplushalf3 = press6_1;
    prnplushalf3 = press7_1;
    plnplusone3 = press5_2;
    pmnplusone3 = press6_2;
    prnplusone3 = press7_2;
    
    uln3 = vel15;
    umn3 = vel16;
    urn3 = vel17;
    ulnplushalf3 = vel15_1; 
    umnplushalf3 = vel16_1;
    urnplushalf3 = vel17_1;
    ulnplusone3 = vel15_2;
    umnplusone3 = vel16_2;
    urnplusone3 = vel17_2;
    
    vln3 = vel25;
    vmn3 = vel26;
    vrn3 = vel27;
    vlnplushalf3 = vel25_1; 
    vmnplushalf3 = vel26_1;
    vrnplushalf3 = vel27_1;
    vlnplusone3 = vel25_2;
    vmnplusone3 = vel26_2;
    vrnplusone3 = vel27_2;

    % parameters for flux through face four

    pln4 = press7;
    pmn4 = press8;
    prn4 = press1;
    plnplushalf4 = press7_1;
    pmnplushalf4 = press8_1;
    prnplushalf4 = press1_1;
    plnplusone4 = press7_2;
    pmnplusone4 = press8_2;
    prnplusone4 = press1_2;
    
    uln4 = vel17;
    umn4 = vel18;
    urn4 = vel11;
    ulnplushalf4 = vel17_1;
    umnplushalf4 = vel18_1;
    urnplushalf4 = vel11_1;
    ulnplusone4 = vel17_2;
    umnplusone4 = vel18_2;
    urnplusone4 = vel11_2;
    
    vln4 = vel27;
    vmn4 = vel28;
    vrn4 = vel21;
    vlnplushalf4 = vel27_1;
    vmnplushalf4 = vel28_1;
    vrnplushalf4 = vel21_1;
    vlnplusone4 = vel27_2;
    vmnplusone4 = vel28_2;
    vrnplusone4 = vel21_2;
    
    % calculating R1, R2 and R3
    
    R1 = -(( vln1 + vrn1 + vlnplusone1 + vrnplusone1 ) + 4*( vmn1 + vmnplusone1 + vlnplushalf1 + vrnplushalf1 ) + 16*vmnplushalf1 )/36 + (( uln2 + urn2 + ulnplusone2 + urnplusone2 ) + 4*( umn2 + umnplusone2 + ulnplushalf2 + urnplushalf2 ) + 16*umnplushalf2 )/36 + (( vln3 + vrn3 + vlnplusone3 + vrnplusone3 ) + 4*( vmn3 + vmnplusone3 + vlnplushalf3 + vrnplushalf3 ) + 16*vmnplushalf3 )/36 - (( uln4 + urn4 + ulnplusone4 + urnplusone4 ) + 4*( umn4 + umnplusone4 + ulnplushalf4 + urnplushalf4 ) + 16*umnplushalf4 )/36;
    
    R2 = (( pln2 + prn2 + plnplusone2 + prnplusone2 ) + 4*( pmn2 + pmnplusone2 + plnplushalf2 + prnplushalf2 ) + 16*pmnplushalf2 )/36 - (( pln4 + prn4 + plnplusone4 + prnplusone4 ) + 4*( pmn4 + pmnplusone4 + plnplushalf4 + prnplushalf4 ) + 16*pmnplushalf4 )/36;
    
    R3 = -(( pln1 + prn1 + plnplusone1 + prnplusone1 ) + 4*( pmn1 + pmnplusone1 + plnplushalf1 + prnplushalf1 ) + 16*pmnplushalf1 )/36 + (( pln3 + prn3 + plnplusone3 + prnplusone3 ) + 4*( pmn3 + pmnplusone3 + plnplushalf3 + prnplushalf3 ) + 16*pmnplushalf3 )/36;
    
    for i=1:sizeu
        press9_2(i) = press9(i) -(dt/area(i))*R1(i);
    end
    
    for i=1:sizeu
        vel19_2(i) = vel19(i) -(dt/area(i))*R2(i);
    end
    
    for i=1:sizeu
        vel29_2(i) = vel29(i) -(dt/area(i))*R3(i);
    end
    
    
    t1_press_iv = press1_2(2:end,2:end);
    t1_vel1_iv = vel11_2(2:end,2:end);
    t1_vel2_iv = vel21_2(2:end,2:end);

    t2_press_iv = press1_2(1,2:end);
    t2_vel1_iv = vel11_2(1,2:end);
    t2_vel2_iv = vel21_2(1,2:end);

    t3_press_iv = press7_2(end,2:end);
    t3_vel1_iv = vel17_2(end,2:end);
    t3_vel2_iv = vel27_2(end,2:end);

    t4_press_iv = press1_2(2:end,1);
    t4_vel1_iv = vel11_2(2:end,1);
    t4_vel2_iv = vel21_2(2:end,1);

    t5_press_iv = press3_2(2:end,end);
    t5_vel1_iv = vel13_2(2:end,end);
    t5_vel2_iv = vel23_2(2:end,end);

    t6_press_iv = press1_2(1,1);
    t6_vel1_iv = vel11_2(1,1);
    t6_vel2_iv = vel21_2(1,1);

    t7_press_iv = press3_2(1,end);
    t7_vel1_iv = vel13_2(1,end);
    t7_vel2_iv = vel23_2(1,end);

    t8_press_iv = press7_2(end,1);
    t8_vel1_iv = vel17_2(end,1);
    t8_vel2_iv = vel27_2(end,1);

    t9_press_iv = press5_2(end,end);
    t9_vel1_iv = vel15_2(end,end);
    t9_vel2_iv = vel25_2(end,end);

    t10_press_iv = press2_2(2:end,:);
    t10_vel1_iv = vel12_2(2:end,:);
    t10_vel2_iv = vel22_2(2:end,:);

    t11_press_iv = press2_2(1,:);
    t11_vel1_iv = vel12_2(1,:);
    t11_vel2_iv = vel22_2(1,:);

    t12_press_iv = press6_2(end,:);
    t12_vel1_iv = vel16_2(end,:);
    t12_vel2_iv = vel26_2(end,:);

    t13_press_iv = press4_2(:,1:end-1);
    t13_vel1_iv = vel14_2(:,1:end-1);
    t13_vel2_iv = vel24_2(:,1:end-1);

    t14_press_iv = press8_2(:,1);
    t14_vel1_iv = vel18_2(:,1);
    t14_vel2_iv = vel28_2(:,1);

    t15_press_iv = press4_2(:,end);
    t15_vel1_iv = vel14_2(:,end);
    t15_vel2_iv = vel24_2(:,end);
    
    
    press9 = press9_2;
    press1 = press1_2;
    press2 = press2_2;
    press3 = press3_2;
    press4 = press4_2;
    press5 = press5_2;
    press6 = press6_2;
    press7 = press7_2;
    press8 = press8_2;
    
    vel19 = vel19_2;
    vel11 = vel11_2;
    vel12 = vel12_2;
    vel13 = vel13_2;
    vel14 = vel14_2;
    vel15 = vel15_2;
    vel16 = vel16_2;
    vel17 = vel17_2;
    vel18 = vel18_2;
    
    vel29 = vel29_2;
    vel21 = vel21_2;
    vel22 = vel22_2;
    vel23 = vel23_2;
    vel24 = vel24_2;
    vel25 = vel25_2;
    vel26 = vel26_2;
    vel27 = vel27_2;
    vel28 = vel28_2;
    
    figure(4)
    surf(X1,Y1,press9_2)
    %view(2)
    colorbar
    drawnow
    
    figure(5)
    surf(X1,Y1,vel19_2)
    %view(2)
    colorbar
    drawnow
    
    figure(6)
    surf(X1,Y1,vel29_2)
    %view(2)
    colorbar
    drawnow
    
end
    
for i=1:sizeu
    exactp(i) = cos(2*pi*tend)*[sin(2*pi*X1(i)) + sin(2*pi*Y1(i))]; 
end

for i=1:sizeu
    exactvel1(i) = sin(2*pi*tend)*cos(2*pi*X1(i)); 
end

for i=1:sizeu
    exactvel2(i) = sin(2*pi*tend)*cos(2*pi*Y1(i)); 
end

forl2 = Delta_x*Delta_y*dt*ones(lx1,lx1);

L2_press=sqrt((1/((xb-xa+2*Delta_x)*(yb-ya+2*Delta_y)*(dt))).*sum(sum(abs(press9_2.*forl2-exactp.*forl2).^2)))

L2_vel1=sqrt((1/((xb-xa+2*Delta_x)*(yb-ya+2*Delta_y)*(dt))).*sum(sum(abs(vel19_2.*forl2-exactvel1.*forl2).^2)))

L2_vel2=sqrt((1/((xb-xa+2*Delta_x)*(yb-ya+2*Delta_y)*(dt))).*sum(sum(abs(vel29_2.*forl2-exactvel2.*forl2).^2)))

ncell = sizeu;
nface = No_of_edges;
nvertex = sizev;
DOF = ncell + nface/2 + nvertex/4
h = 1/sqrt(DOF)    