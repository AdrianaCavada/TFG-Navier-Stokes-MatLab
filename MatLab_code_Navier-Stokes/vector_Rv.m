function [Rv] = vector_Rv(U, V, Re, dx, dy, nx, ny)

% R(v)=aWv*VW+aEv*VE+aNv*VN+aSv*VS-aPv*VP
% The matrices of coefficients have the same dimension as V
aWv=zeros(ny+1,nx+2);
aEv=zeros(ny+1,nx+2);
aNv=zeros(ny+1,nx+2);
aSv=zeros(ny+1,nx+2);
Rv=zeros(ny+1,nx+2);

for i=2:ny
    for j=2:nx
%         aNv(i,j)=(dx/Re/dy)-0.25*dx*((V(i,j)+V(i+1,j))*0.5+V(i+1,j)); % VP+VN
%         aSv(i,j)=(dx/Re/dy)+0.25*dx*((V(i,j)+V(i+1,j))*0.5+V(i-1,j)); % VP+VS
%         aEv(i,j)=(dy/Re/dx)-0.25*dy*((U(i,j)+U(i+1,j))*0.5+U(i+1,j)); % UP+UN
        aNv(i,j)=(dx/Re/dy)-0.25.*dx.*(V(i,j)+V(i+1,j)); % VP+VN
        aSv(i,j)=(dx/Re/dy)+0.25.*dx.*(V(i,j)+V(i-1,j)); % VP+VS
        aEv(i,j)=(dy/Re/dx)-0.25.*dy.*(U(i,j)+U(i+1,j)); % UN+UP
        aWv(i,j)=(dy/Re/dx)+0.25.*dy.*(U(i+1,j-1)+U(i,j-1)); % UWN+UW
        
%         aNpv(i,j)=(dx/Re/dy)+0.25.*dx.*(V(i,j)+V(i+1,j)); % VP+VN
%         aSpv(i,j)=(dx/Re/dy)-0.25.*dx.*(V(i,j)+V(i-1,j)); % VP+VS
%         aEpv(i,j)=(dy/Re/dx)+0.25.*dy.*(U(i,j)+U(i+1,j)); % UN+UP
%         aWpv(i,j)=(dy/Re/dx)-0.25.*dy.*(U(i+1,j-1)+U(i,j-1)); % UWN+UW
        
    end
end 
aPv=aWv+aEv+aNv+aSv;
% aPv=aWpv+aEpv+aNpv+aSpv;

for i=2:ny
    for j=2:nx
        Rv(i,j)=aWv(i,j).*V(i,j-1)+aEv(i,j).*V(i,j+1)+aNv(i,j).*V(i+1,j)+aSv(i,j).*V(i-1,j)-aPv(i,j).*V(i,j);
    end
end
end

