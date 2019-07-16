function [U_star, V_star, Ru, Rv] = midway_velocity_Euler(U, V, Re, dx, dy, nx, ny, dt)

Ru=vector_Ru(U,V,Re,dx,dy,nx,ny);
Rv=vector_Rv(U,V,Re,dx,dy,nx,ny);

U_star=U+dt*Ru;
V_star=V+dt*Rv;

end

