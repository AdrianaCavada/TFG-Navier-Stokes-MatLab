function [U_star, V_star, Ru, Rv] = midway_velocity_AdamB(U, V, Re, dx, dy, nx, ny, dt, Ru_old, Rv_old)

Ru=vector_Ru(U,V,Re,dx,dy,nx,ny);
Rv=vector_Rv(U,V,Re,dx,dy,nx,ny);

U_star=U+dt*(1.5*Ru-0.5*Ru_old);
V_star=V+dt*(1.5*Rv-0.5*Rv_old);
% U_star=U+dt*Ru_old;
% V_star=V+dt*Rv_old;


end

