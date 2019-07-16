function [A,aW,aE,aN,aS] = coefficient_pressure_matrix(nx,ny,dx,dy)
% % Create the coefficient matrix (A)to solve the pressure 
%   A*p=div(velocity_star)/dt
%   the dimension of the matrix is (ny*nx,ny*nx)
%   A has five diagonals, one for each coefficient (aW,aE,aN,aS and aP)
%   in this case, as the mesh is uniform, the aW-aE and aN-aS coefficients
%   have the same value
%
%   Then:
%       aWE=dy/dx;  aNS=dx/dy;  aP=2*(aWE+aNS);

% % Create the empty matrix
A=spalloc(ny*nx,ny*nx,5*(ny*nx-2*nx));
% % Create diagonals 
aW=spalloc(ny,nx,ny*nx-2*nx);
aE=spalloc(ny,nx,ny*nx-2*nx);
aN=spalloc(ny,nx,ny*nx-2*nx);
aS=spalloc(ny,nx,ny*nx-2*nx);
for i=1:ny
    for j=1:nx
        aW(i,j)=dy/dx;
        aE(i,j)=dy/dx;
        aN(i,j)=dx/dy;      
        aS(i,j)=dx/dy;     
    end
end

aN(end,:)=0;
aS(1,:)=0;
aP=aW+aE+aN+aS;
aW=aW(:);
aE=aE(:);
aN=aN(:);
aS=aS(:);
aP=aP(:);

A=diag(-aP,0)+diag(aN(1:end-1),1)+diag(aS(2:end),-1)+diag(aE(1:end-ny),ny)+diag(aW(ny+1:end),-ny);

end

