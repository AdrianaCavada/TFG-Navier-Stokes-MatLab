function [Ru] = vector_Ru(U, V, Re, dx, dy, nx, ny)

% R(u)=aWu*UW+aEu*UE+aNu*UN+aSu*US-aPu*UP
% The matrices of coefficients have the same dimension as U
aWu=zeros(ny+2,nx+1);
aEu=zeros(ny+2,nx+1);
aNu=zeros(ny+2,nx+1);
aSu=zeros(ny+2,nx+1);
Ru=zeros(ny+2,nx+1);

for i=2:ny+1
    for j=2:nx
%         aWu(i,j)=(dy/Re/dx)+0.25*dy*((U(i,j)+U(i,j+1))*0.5+U(i,j-1));  % UP+UW
%         aEu(i,j)=(dy/Re/dx)-0.25*dy*((U(i,j)+U(i,j+1))*0.5+U(i,j+1)); % UP+UE
%         aNu(i,j)=(dx/Re/dy)-0.25*dx*((V(i,j)+V(i,j+1))*0.5+V(i,j+1)); % VP+VE ok
        aWu(i,j)=(dy/Re/dx)+0.25*dy.*(U(i,j)+U(i,j-1));  % UP+UW
        aEu(i,j)=(dy/Re/dx)-0.25*dy.*(U(i,j)+U(i,j+1)); % UP+UE
        aNu(i,j)=(dx/Re/dy)-0.25*dx.*(V(i,j+1)+V(i,j)); % VP+VE ok
        aSu(i,j)=(dx/Re/dy)+0.25*dx.*(V(i-1,j)+V(i-1,j+1)); % VS+VES ok  
        
%         aWpu(i,j)=(dy/Re/dx)-0.25*dy.*(U(i,j)+U(i,j-1));  % UP+UW
%         aEpu(i,j)=(dy/Re/dx)+0.25*dy.*(U(i,j)+U(i,j+1)); % UP+UE
%         aNpu(i,j)=(dx/Re/dy)+0.25*dx.*(V(i,j+1)+V(i,j)); % VP+VE ok
%         aSpu(i,j)=(dx/Re/dy)-0.25*dx.*(V(i-1,j)+V(i-1,j+1)); % VS+VES ok 
    end
end 
aPu=aWu+aEu+aNu+aSu;
% aPu=aWpu+aEpu+aNpu+aSpu;

% Multiply each coefficient by its corresponding U point (UW, UE...etc)
for i=2:ny+1
    for j=2:nx
        Ru(i,j)=aWu(i,j).*U(i,j-1)+aEu(i,j).*U(i,j+1)+aNu(i,j).*U(i+1,j)+aSu(i,j).*U(i-1,j)-aPu(i,j).*U(i,j);
    end
end

end

