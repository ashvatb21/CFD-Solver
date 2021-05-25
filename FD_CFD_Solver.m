clear all
clc
close all

% Input Parameters
nx      =   50;                             % Number of cells in x-direction
ny      =   50;                             % Number of cells in y-direction
Lx      =   2;                              % Length in x-direction
Ly      =   2;                              % Length in y-direction
nu      =   0.1;                            % Kinematic Viscosity
rho     =   1;                              % Density
t_final =   5;                              % Time Duration

% Boundary Conditions
u_top       =   1;                          % Velocity at (:, y = Ly)
u_bottom    =   0;                          % Velocity at (:, y = 0)
v_left      =   0;                          % Velocity at (x = Lx, :)
v_right     =   0;                          % Velocity at (x = Lx, :)

% Index Extends
i_min    =   2;                             % Min index in i-direction
i_max    =   i_min + nx - 1;                % Max index in i-direction
j_min    =   2;                             % Min index in j-direction
j_max    =   j_min + ny - 1;                % Max index in j-direction

% Create Mesh Grid
x(i_min:i_max+1)  =   linspace(0, Lx, nx + 1);                      % Left faces 
y(j_min:j_max+1)  =   linspace(0, Ly, ny + 1);                      % Bottom faces
xm(i_min:i_max)   =   0.5*(x(i_min:i_max) + x(i_min+1:i_max+1));    % i cell middle
ym(j_min:j_max)   =   0.5*(y(j_min:j_max) + y(j_min+1:j_max+1));    % j cell middle

% Create Mesh Sizes
dx      =   x(i_min+1) - x(i_min);             % Length of 1 cell
dy      =   y(j_min+1) - y(j_min);             % Height of 1 cell
dxi     =   1/dx;                          
dyi     =   1/dy;                          

% Declarations of other Parameters  
t       =   0;
dt      =   0.0025;
u       =   zeros(i_max+1, j_max+1);
v       =   zeros(i_max+1, j_max+1);
u_star  =   zeros(i_max+1, j_max+1);
v_star  =   zeros(i_max+1, j_max+1);
R       =   zeros(1, nx*nx);
pv      =   zeros(1, nx*nx);
p       =   zeros(nx+1, ny+1);

% Laplace Operator
L = zeros(nx*ny, nx*ny) ;

for j = 1:ny
    for i = 1:nx
        L(i+(j-1)*nx, i+(j-1)*nx) = 2*dxi^2 + 2*dyi^2;
        
        for ii = i-1: 2:i+1
            if(ii > 0 && ii <= nx) % Interior point(s)
                L(i+(j-1)*nx, ii+(j-1)*nx) = -dxi^2;
            else % Neumann conditions on boundary
                L(i+(j-1)*nx, i+(j-1)*nx) = L(i+(j-1)*nx, i+(j-1)*nx) - dxi^2;
            end
        end
        
        for jj = j-1: 2: j+1
            if(jj > 0 && jj <= ny) % Interior point(s)
                L(i+(j-1)*nx, i+(jj-1)*nx) = -dyi^2;
            else % Neumann conditions on boundary
                L(i+(j-1)*nx, i+(j-1)*nx) = L(i+(j-1)*nx, i+(j-1)*nx) - dyi^2;
            end
        end
        
    end
end

% Pressure in First Cell
L(1,:) = 0 ; L(1,1)= 1;

while t <= t_final
    
    t = t + dt;
    
    disp('Current Time')
    disp(t)
    
    % Boundary Conditions
    u(:,j_min-1) = u(:,j_min) - 2*(u(:,j_min) - u_bottom);
    u(:,j_max+1) = u(:,j_max) - 2*(u(:,j_max) - u_top);
    v(i_min-1,:) = v(i_min,:) - 2*(v(i_min,:) - v_left);
    v(i_max+1,:) = v(i_max,:) - 2*(v(i_max,:) - v_right);
    
    % u_star calc %
    for j = j_min:j_max
        for i = i_min+1:i_max
            v_curr = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j)+v(i,j+1));
            
            u_star(i,j) = u(i,j) + dt*(nu*(u(i-1,j) - 2*u(i,j)+u(i+1,j))*dxi^2 ...
                + nu*(u(i,j-1) - 2*u(i,j) + u(i,j+1))*dyi^2 ...
                - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi...
                - v_curr*(u(i,j+1) - u(i,j-1))*0.5* dyi);
        end
    end
    
    
    % v_star calc %
    for j = j_min+1:j_max
        for i = i_min:i_max
            u_curr = 0.25*(u(i,j-1) + u(i,j) + u(i+1,j-1) + u(i+1,j));
            
            v_star(i,j) = v(i,j) + dt*(nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi^2 ...
                + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi^2 ...
                - u_curr*(v(i+1,j) - v(i-1,j))*0.5*dxi...
                - v(i,j)*(v(i,j+1) - v(i,j-1))*0.5* dyi);
        end
    end
    
    % R calc
    n=0;
    
    for j = j_min:j_max
        for i = i_min:i_max
            n = n+1;
            R(n) = -(rho/dt)*((u_star(i+1,j) - u_star(i,j))*dxi + (v_star(i,j+1) - v_star(i,j))*dyi);
        end
    end
    
    % Pressure
    pv = L\(R');
    n = 0;
    p = zeros(i_max,j_max);
    
    for j = j_min:j_max
        for i = i_min:i_max
            n = n+1;
            p(i,j) = pv(n);
        end
    end
    
    % Corrector Step
    for j = j_min:j_max
        for i = i_min+1:i_max
            u(i,j) = u_star(i,j) - (dt/rho)*(p(i,j) - p(i-1,j))*dxi ;
        end
    end
    for j = j_min+1:j_max
        for i = i_min:i_max
            v(i,j) = v_star(i,j) - (dt/rho)*(p(i,j)-p(i,j-1))*dyi ;
        end
    end

end


%% Plotting
[X,Y] = meshgrid(x,y);

figure(1)
[c,h] = contourf(X,Y,u',25)
set(h, 'edgecolor','none');
title('Velocity Component, u')
h1 = colorbar;
set(get(h1,'label'),'string','u [m/s]'); %title the contour colorbar
ylabel('y')
xlabel('x')
zlabel('u')
colormap(jet)


figure(2)
[c,h] = contourf(X,Y,v',25)
set(h, 'edgecolor','none');
title('Velocity Component, v')
h1 = colorbar;
set(get(h1,'label'),'string','v [m/s]'); %title the contour colorbar
ylabel('y')
xlabel('x')
zlabel('v')
colormap(jet)

velocity = sqrt(u.^2 + v.^2);
figure(3)
[c,h] = contourf(X,Y,velocity',25)
set(h, 'edgecolor','none');
title('Net Velocity')
h1 = colorbar;
set(get(h1,'label'),'string','Velocity [m/s]'); %title the contour colorbar
ylabel('y')
xlabel('x')
zlabel('U')
colormap(jet)


[XM,YM] = meshgrid(xm,ym)

figure(4)
[c,h] = contourf(XM,YM,p',25)
set(h, 'edgecolor','none');
title('Pressure Distribution')
h1 = colorbar;
set(get(h1,'label'),'string','P [Pa]'); %title the contour colorbar
ylabel('y')
xlabel('x')
zlabel('U')
colormap(jet)
