function [div] = diver(U,V, nx, ny, dx, dy)
% % Calculate the divergence  of a vector field in 2D 

div=zeros(ny,nx);

div=0.5*(dx*(U(2:ny+1,2:nx+1)-U(2:ny+1,1:nx))+dy*(V(2:ny+1,2:nx+1)-V(1:ny,2:nx+1)));

