function [ Kx_grid ] = Permeability( N_grids )
%% This function calculates the apparent permeability  
Kx_grid=zeros(N_grids,1);
width_grid=0.25*cos(pi/2.8376);
for i=1:N_grids
    yi=width_grid*(i-1); yip1=width_grid*i;
    Kx_grid(i)=-(1/6)*((yip1^2)+(yi^2)+(yip1*yi))+1.342*(yi+yip1);
end
Kx_grid=Kx_grid*9.413092183E+13;
end

