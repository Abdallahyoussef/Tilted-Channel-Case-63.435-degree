%% This is general code for flow simulation with incompressible fluid using S-B equation
% All the units are in field units 
% constant presssure boundaries
% The pressure is defined at the first and last grid cell to use block
% centerd griding as way of approximation
close all
clear 
clc       
%% Reservoir parmters 
dx=0.5;% length by ft of each grid in x direction  
dy=0.25; % length by ft of each grid in y direction
dz=15;                                        % length by ft of each grid in z direction
lx=30;                                       % number of grids in x direction 
ly=150;                                        % number of grids in y direction 
[Kx_cave ] = Permeability( 48 );
N=lx*ly;                                 % number of grids in y direction
Ax=dy*dz;
Ay=dx*dz;
Kx=200*ones(ly,lx);
for i=1:lx
    if i==1
        Kx(1:length(Kx_cave),1)=Kx_cave;
    else 
        Kx(4*(i-1):4*(i-1)+length(Kx_cave)-1,i)=Kx_cave;
    end
end
Kx=Kx';Kx=Kx(:);
for i=1:length(Kx)
    if Kx(i)>1000
        Kx(i)=inf;
    end
end
Ky=Kx;
Bc=1.127E-3; B2c = 1.0624E-14;
%% boundary conditions 
p_left=1000; % pressure at left boundary
p_right=500; % pressure at right boundary
%% Calculations 
Tx=Bc*Kx*Ax/dx;  Ty=Bc*Ky*Ay/dy; % transmissibility
Cx=B2c*dx/(Bc*Ax); % shear term in X-direction
Cy=B2c*dy/(Bc*Ay); % shear term in Y-direction
A = spalloc(N+(lx+1)*ly+(ly-1)*lx,N+(lx+1)*ly+(ly-1)*lx,9*N + 7*ly*(lx+1) + 7*lx*(ly-1)); % pressure rate matrix;
RHS= zeros(N+(lx+1)*ly+(ly-1)*lx,1);
for m=1:N
    j=ceil(m/lx);
    i=m-(j-1)*lx;
    if i > 1
        A(m,N+i+(j-1)*(lx+1))=-1;  % Conservation of mass flow rate in from m-1 to m
    end
    if i<lx
        F=N+i+1+(j-1)*(lx+1);
        A(m,F)=1;       % Conservation of mass flow rate out from m to m+1 
        A(F,m)=-1;            % Conservation of momentuem pressure at grid m
        A(F,m+1)=1;           % Conservation of momentuem pressure at grid m+1
        A(F,F)=0.5*((1/Tx(m))+(1/Tx(m+1)))+2*Cx*((1/(dx^2))+(1/(dy^2)));%+(4/(dz^2))); % Conservation of momentuem Dary's part  
        A(F,F+1)=-Cx/(dx^2);
        A(F,F-1)=-Cx/(dx^2);
        if j==1
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F+(lx+1))=-(4/3)*Cx*(1/(dy^2));
        elseif j==ly
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F-(lx+1))=-(4/3)*Cx*(1/(dy^2));
        else
            A(F,F+(lx+1))=-Cx*(1/(dy^2));
            A(F,F-(lx+1))=-Cx*(1/(dy^2));
        end
    end
    if j>1
        A(m,N+ly*(lx+1)+m-lx)=-1;
    end
    if j<ly
        F=N+ly*(lx+1)+m;
        A(m,F)=1;
        A(F,m)=-1;
        A(F,m+lx)=1;
        A(F,F)=0.5*((1/Ty(m))+(1/Ty(m+lx)))+2*Cy*((1/(dx^2))+(1/(dy^2)));%+(4/(dz^2)));
        if j<(ly-1)
            A(F,F+lx)=-Cy/(dy^2);
        end
        if j>1
            A(F,F-lx)=-Cy/(dy^2);
        end
        if i==1
            A(F,F+1)=-(4/3)*Cy*(1/(dx^2));
            A(F,F)=A(F,F)+2*Cy*(1/(dx^2));
        elseif i==lx
            A(F,F-1)=-(4/3)*Cy*(1/(dx^2));
            A(F,F)=A(F,F)+2*Cy*(1/(dx^2));
        else
            A(F,F+1)=-Cy*(1/(dx^2));
            A(F,F-1)=-Cy*(1/(dx^2));
        end
    end
    if i==1  % S-B equation for left boundary (constant pressure and approximating flow rate)
        F=N+i+(lx+1)*(j-1);
        A(F,m)=2;   % pressure at i=1
        A(m,F)=-1; % rate at boundary
        RHS(F)=2*p_left; % constant pressure boundary
        A(F,F)=(1/Tx(m))+Cx*((1/(dx^2))+(2/(dy^2)));
        A(F,F+1)=-Cx/(dx^2);
        if j==1
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F+(lx+1))=-(4/3)*Cx*(1/(dy^2));
        elseif j==ly
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F-(lx+1))=-(4/3)*Cx*(1/(dy^2));
        else
           A(F,F+(lx+1))=-Cx*(1/(dy^2));
           A(F,F-(lx+1))=-Cx*(1/(dy^2));
        end
    end
    if i==lx  % S-B equation for right boundary (constant pressure and approximating flow rate)
        F=N+i+1+(lx+1)*(j-1);
        A(F,m)=-2;   % pressure at i=1
        A(m,F)=1;   % rate at boundary
        RHS(F)=-2*p_right; % constant pressure boundary
        A(F,F)=(1/Tx(m))+Cx*((1/(dx^2))+(2/(dy^2)));
        A(F,F-1)=-Cx/(dx^2);
        if j==1
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F+(lx+1))=-(4/3)*Cx*(1/(dy^2));
        elseif j==ly
            A(F,F)=A(F,F)+2*Cx*(1/(dy^2));
            A(F,F-(lx+1))=-(4/3)*Cx*(1/(dy^2));
        else
           A(F,F+(lx+1))=-Cx*(1/(dy^2));
           A(F,F-(lx+1))=-Cx*(1/(dy^2));
        end
    end
end
Solution=A\RHS;
Pressure=Solution(1:N); % pressure distribution
Fx=Solution(N+1:N+ly*(lx+1))*(3.8993E-3)/(dy*dz); % velocity in X direction 
Fy=Solution(N+ly*(lx+1)+1:end)*(3.8993E-3)/(dx*dz); % velocity in Y direction 
% pressure plot
% surf(vec2mat(Pressure,lx)); view(2);  shading interp % ploting pressure
% figure(2)
% surf(vec2mat(Fx,lx+1)); view(2);  shading interp % ploting rate in X direction 
% figure(3)
% surf(vec2mat(Fy,lx-1)); view(2);  shading interp % ploting rate in Y direction
Fx(1:lx+1:end)=[];
Fy=[Fy;zeros(lx,1)];
[x,y]=meshgrid(dx/2:dx:dx*(lx-0.5),dy/2:dy:dy*(ly-0.5));
figure(4)
contourf(x,y,vec2mat(Pressure,lx))
title('Pressure Map Stokes-Brinkman')
% y=y';y=y(:);x=x';x=x(:);
figure(5)
velocity=quiver(x,y,vec2mat(Fx,lx),vec2mat(Fy,lx));  % ploting net rate map
set(velocity,'linewidth',2,'color','r')
startx=[(dx/2)*ones(1,ly/10) 0.5*dx:dx:dx*(lx-0.5)]; starty=[dy/2:10*dy:dy*(ly-0.5) linspace(dy/2,30,lx)]; hlines=streamline(x,y,vec2mat(Fx,lx),vec2mat(Fy,lx),startx,starty);
set(hlines,'linewidth',1,'color','b')
title('Stramline Map Stokes-Brinkman')