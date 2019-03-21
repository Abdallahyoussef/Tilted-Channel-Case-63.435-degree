%% This is a code for Tilted Channel Flow Simulation using Gloabally Coupled Pressure Method
% The channel is tilted by angle theta and constant pressure boundaries at left and right boundaries
% Isothermal flow and ideal incopressible system
% This code uses GCP-2
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
Kx=200*ones(ly,lx);
for i=1:lx
    if i==1
        Kx(1:length(Kx_cave),1)=Kx_cave;
    else 
        Kx(4*(i-1):4*(i-1)+length(Kx_cave)-1,i)=Kx_cave;
    end
end
Kx=Kx';Kx=Kx(:);
Ky=Kx;
Theta=pi/2.8376; % inclination angel of channel by radian
Bc=1.127E-3; % conversion factor
P_left=1000; % pressure at left boundary
P_right=500; % pressure at right boundary
%% Calculations
COSIN=cos(Theta); SIN=sin(Theta); % sin and cos of rotation angle 
Ax=dy*dz;  % area in X direction
Ay=dx*dz;  % area in Y direction
Kxx=Kx*(COSIN^2)+Ky*(SIN^2); % Kxx permeability in full tensor
Kxy=(Kx-Ky)*COSIN*SIN; % Kxy and Kyx permeability in full tensor
Kyy=Ky*(COSIN^2)+Kx*(SIN^2); % Kyy permeability in full tensor
Nx=(lx-1)*ly; % number of boundaries in X direction
Nf=Nx+(ly+1)*lx; % number of axuilary unknowns 
B=spalloc(Nf,Nf,3*Nf);  % axuilary pressure matrix
C=spalloc(Nf,N,3*Nf); % primary pressure matrix
Boundary_pressure=zeros(Nf,1); % specifeid pressure matrix that is added to primary pressures
D= spalloc((lx+1)*ly+(ly-1)*lx,Nf,3*((lx+1)*ly+(ly-1)*lx));  % rate axuilary pressure matrix
E= spalloc((lx+1)*ly+(ly-1)*lx,N,3*((lx+1)*ly+(ly-1)*lx));   % rate primary pressure matrix
Boundary_rate=zeros((lx+1)*ly+(ly-1)*lx,1); % specifeid pressure matrix that is added to rate at boundaries
for m=1:N  % constructing primary and axuilary unknowns matrices
    j=ceil(m/lx); % index in Y direction
    i=m-(j-1)*lx; % index in X direction
    % start of forming rate matrices in X direction
    f_rate=i+(j-1)*(lx+1);  % index of rate at inteface between m and m-1
    E(f_rate,m)=-2*((Kxx(m)/dx)+(Kxy(m)/dy)); % rate primary pressure coefficient at m-1
    D(f_rate,m+ly*(lx-1))=(2/dy)*Kxy(m); % rate axuilary pressure coefficient at inteface between m and m-I
    if i > 1
        f=i-1+(j-1)*(lx-1); % index of intefrace between m and m-1 
        B(f,f)=(Kxx(m)+Kxx(m-1))*(2/dx);  % pressure coefficeint at interface between m and m-1
        B(f,ly*(lx-1)+i+(j-1)*lx)=-2*Kxy(m)/dy; % pressure coefficeint at interface between m and m+I
        B(f,ly*(lx-1)+i-1+j*lx)=2*Kxy(m-1)/dy; % pressure coefficeint at interface between m-1 and m-1+I
        C(f,m-1)=2*((Kxx(m-1)/dx)+(Kxy(m-1)/dy));  % pressure coefficeint m-1
        C(f,m)=2*((Kxx(m)/dx)+(Kxy(m)/dy));  % pressure coefficeint m
        D(f_rate,i-1+(j-1)*(lx-1))=(2/dx)*Kxx(m); % rate axuilary pressure coefficient at inteface between m and m-1
    end
    if j>1
        f=ly*(lx-1)+i+(j-1)*lx; % index of intefrace between m and m-I 
        C(f,m-lx)=2*((Kyy(m-lx)/dy)+(Kxy(m-lx)/dx));  % pressure coefficeint m-I
        C(f,m)=2*((Kyy(m)/dy)+(Kxy(m)/dx));  % pressure coefficeint m
        B(f,f)=(Kyy(m)+Kyy(m-lx))*(2/dy);  % pressure coefficeint at interface between m and m-I
        % dealing with specified pressure boundaries
        if i==1
            Boundary_pressure(f)=-2*Kxy(m)*P_left/dy; % specifed pressure at left boundary
            B(f,i+(j-2)*(lx-1))=2*Kxy(m-lx)/dy; % pressure coefficeint at interface between m-I and m-I+1
        elseif i==lx
            Boundary_pressure(f)=-2*Kxy(m-lx)*P_right/dy; % specifed pressure at right boundary
            B(f,i-1+(j-1)*(lx-1))=2*Kxy(m)/dy; % pressure coefficeint at interface between m and m-1
        else
            B(f,i+(j-2)*(lx-1))=2*Kxy(m-lx)/dy; % pressure coefficeint at interface between m-I and m-I+1
            B(f,i-1+(j-1)*(lx-1))=2*Kxy(m)/dy; % pressure coefficeint at interface between m and m-1
        end
    end
    if j<ly
        f_rate=(lx+1)*ly+m;  % index of rate at inteface between m and m+I
        E(f_rate,m)=2*((Kyy(m)/dy)+(Kxy(m)/dx)); % rate primary pressure coefficient at m-I
        D(f_rate,ly*(lx-1)+i+lx*j)=-(2/dy)*Kyy(m); % rate axuilary pressure coefficient at inteface between m-I and m-I+1
        if i==lx
            Boundary_rate(f_rate)=-2*Kxy(m)*P_right/dx; % specifed pressure at right boundary
        else
            D(f_rate,i+(j-1)*(lx-1))=-(2/dy)*Kxy(m); % rate axuilary pressure coefficient at inteface between m and m-1
        end
    end
    if j==1 % no flow boundary at South 
        f=(lx-1)*ly+i;  % index of South boundary pressure
        C(f,m)=2*((Kxy(m)/dx)+(Kyy(m)/dy));  % pressure coefficeint m
        B(f,f)=2*Kyy(m)/dy; % pressure coefficeint at South interface
        if i==1
            Boundary_pressure(f)=-2*Kxy(m)*P_left/dx; % specifed pressure at left boundary
        else
            B(f,i-1)=2*Kxy(m)/dx; % pressure coefficeint at interface between m and m-1
        end
    end
    if j==ly % no flow boundary at Norht 
        f=ly*(2*lx-1)+i;  % index of North boundary pressure
        C(f,m)=2*((Kxy(m)/dx)+(Kyy(m)/dy));  % pressure coefficeint m
        B(f,f)=2*Kyy(m)/dy; % pressure coefficeint at south interface
        if i==lx
            Boundary_pressure(f)=-2*Kxy(m)*P_right/dx; % specifed pressure at right boundary
        else
            B(f,i+(lx-1)*(ly-1))=2*Kxy(m)/dx; % pressure coefficeint at interface between m and m-1
        end
    end
    if i==lx % flow at right boundary
        f_rate=j*(lx+1); % rate index at right boundary
        Boundary_rate(f_rate)=-2*Kxx(m)*P_right/dx; % specifed pressure at right boundary
        E(f_rate,m)=2*((Kxx(m)/dx)+(Kxy(m)/dy));  % rate primary pressure coefficient at m
        D(f_rate,ly*(lx-1)+i+j*lx)=-(2/dy)*Kxy(m); % rate axuilary pressure coefficient at inteface between m and m+I        
    end
    if i==1 % flow at left boundary
        f_rate=i+(j-1)*(lx+1); % rate index at left boundary
        Boundary_rate(f_rate)=2*Kxx(m)*P_left/dx; % specifed pressure at left boundary
    end
end
%% Constructing Total Transmissiblity Matrix
T_p=D*(B\C)+E;  T_free=D*(B\Boundary_pressure)+Boundary_rate;  % two components of static transmissibility 
%% Construction of conservation of mass mtrix
AP= spalloc(N,(lx+1)*ly+(ly-1)*lx,4*N);  % rate coefficient matrix 
for m=1:N;  % constructing mass balance matric
    j=ceil(m/lx); % index in Y direction
    i=m-(j-1)*lx; % index in X direction
    % flow in X direction
    AP(m,i+(j-1)*(lx+1))=-1;  % Conservation of mass flow rate in from m-1 to m
    AP(m,i+1+(j-1)*(lx+1))=1;  % Conservation of mass flow rate in from m to m+1
    % flow in Y direction
    if j>1
        AP(m,ly*(lx+1)+m-lx)=-1;
    end
    if j<ly
        AP(m,ly*(lx+1)+m)=1;
    end
end
LHS=AP*T_p;   RHS=-AP*T_free;  % LHS and RHS for pressure coefficient matrices for pressure in conservation of mass
Pressure=LHS\RHS; % pressure at centre of grid blocks 
Flow_rate=T_p*Pressure+T_free;  % flow rate at intefaces
% figure(1)
% surf(vec2mat(Pressure,lx)); %view(2);
% shading interp % ploting pressure 
Fx=Flow_rate(1:(lx+1)*ly); Fy=Flow_rate(1+(lx+1)*ly:end);
% % figure(2)
% % surf(vec2mat(Fx,lx+1)); view(2);  shading interp % ploting rate in X direction 
% figure(3)
% surf(vec2mat(Fy,lx-1)); view(2);  shading interp % ploting rate in Y direction
Fx(1:lx+1:end)=[];
Fy=[Fy;zeros(lx,1)];
[x,y]=meshgrid(dx/2:dx:dx*(lx-0.5),dy/2:dy:dy*(ly-0.5));
figure(4)
contourf(x,y,vec2mat(Pressure,lx))
title('Pressure Map DMEPD')
% y=y';y=y(:);x=x';x=x(:);
figure(5)
velocity=quiver(x,y,vec2mat(Fx,lx),vec2mat(Fy,lx));  % ploting net rate map
set(velocity,'linewidth',2,'color','r')
startx=[(dx/2)*ones(1,ly/10) 0.5*dx:dx:dx*(lx-0.5)]; starty=[dy/2:10*dy:dy*(ly-0.5) linspace(dy/2,30,lx)]; hlines=streamline(x,y,vec2mat(Fx,lx),vec2mat(Fy,lx),startx,starty);
set(hlines,'linewidth',1,'color','b')
title('Stramline Map DMEPD')
