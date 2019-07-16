function [grad_p_U,grad_p_V, p,div_velocity_star] = pressure_correction( U_star, V_star, A, nx,ny, dx, dy, dt,tolerance_p,n_steps_p, p)
% % Calculate the pressure by solving the Poisson equation and
%   return the gradient of the calculated pressure 

% % Calculate the divergence of the velocity_star
div_velocity_star=diver(U_star,V_star,nx, ny, dx, dy);
%contourf(div_velocity_star)
%pause


p_prima=p(2:end-1,2:end-1);
% % Calculate the pressure

p_prima=bicgstab(A,div_velocity_star(:)/dt,tolerance_p,n_steps_p);

p_prima=reshape(p_prima,ny,nx);
% p=p_initial+p;

p(2:end-1,2:end-1)=p_prima;

p(:,1) = p(:,2); % If the first and the second columns have the same value, the derivative of the first column is zero and the same in the rest
p(:,end) = p(:,end-1);
 p(1,:) = p(2,:);
p(end,:) = p(end-1,:);

% Calculate the gradient of pressure
   grad_p_U=0.5*(dx*(p(1:end,2:end)-p(1:end,1:end-1)));
   grad_p_V=0.5*dy*(p(2:end,1:end)-p(1:end-1,1:end));


end

