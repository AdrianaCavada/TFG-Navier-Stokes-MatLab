function [U,V,p_final,errorU,errorV] = navier_stokes(n_step_snst, nx, ny, dt) 
format long
tic
% % Physical variables
%       rho = 933; % density 
%       mu = 0.2390; % dinamic viscosity 
%       U_in = 0.1;  % modulus of initial velocity (m/s)
%       D_equivalent = (Lx*Ly)*4/(Lx+Lx+Ly+Ly); % equivalent diameter for the reynolds number
%       Re = rho*D_equivalent*v_in/mu; % Reynolds number
        Re = 100; % The reynolds number has to be less than 2300 for the laminar regime
        v0=1;        
        
% % Grid variables
%       nx = 20; % nx is the number of cells in x direction
%       ny = 20;  % ny is the number of cells in y direction
      Lx = 1;  % Lx is the domain length in x direction
      Ly = 1;  % Ly is the domain length in y direction
      dx = Lx/nx;  % dx is the cell size in x direcction
      dy = Ly/ny;  % dy is the cell size in y direcction
%       dx = 0.001;
%       dy = 0.001;
      
% % Time integration variables
% dt=5e-3; 
% n_steps=100000;
tolerance_p=0.001;
n_steps_p=3000;

% Plotting variables
makevideo=1;
% n_snapshot=100;
% n_step_snst is the number of steps between each snapshot
if (makevideo==1)
  h=figure('Position',[1 0 800 600],'MenuBar','none','ToolBar','none','resize','off','Visible','off'); % 800 x 600 hopefully
  vidobj = VideoWriter('fluid2.avi');
  open(vidobj);      
  axis tight;
  set(gca, 'CLim',[0,1],'position',[0 0 1 1],'Visible','off');
end 
% % Initial matrices
U=zeros(ny+2,nx+2); 
V=zeros(ny+2,nx+2);
p=zeros(ny+2,nx+2);

% % Bounsary conditions

U(:,1)= 0; V(:,1)=0; % east
U(:,end)= 0; V(:,end)=0; % west
U(1,:)= 0; V(1,:)=0; % south
U(end,:)= 1; V(end,:)=0; % north                 
  
% iter = 1 ---> in the first iteration the scheme for time discretization
% is explicit Euler (velocity_star = v_previous+dt*R(velocity))
%       R(velocity) is the discretizated vector of convecction and
%       difussion parts
% If iter > 1 ---> the time discretization scheme is Adam-Bashford
% (velocity_star = v_previous+dt*[1.5*R(velocity)- 0.5*R(velocity)])
iter = 1;
% % Calculate the coefficient matrix to solve the pressure
A=coefficient_pressure_matrix(nx,ny,dx,dy);
n_step=1;
errorU=1;
errorV=1;
while (errorV > 5e-6)
    
    % % Calculate the mid-way velocity

    if n_step == 1   % % Explicit Euler
        U=avg(U);
        V=avg(V')';
        [U_star, V_star, Ru_old, Rv_old] = midway_velocity_Euler(U, V, Re, dx, dy, nx, ny, dt);
       
        [grad_p_U,grad_p_V, p,div_velocity_star] = pressure_correction( U_star, V_star, A, nx,ny, dx, dy, dt,tolerance_p,n_steps_p, p);
        grad_p_U=grad_p_U.*dt;
        grad_p_V=grad_p_V.*dt;
        % % Calculate the velocity at step n+1 by correcting it with the
        % previous divergence of p 
        U = U_star - grad_p_U;
        V = V_star - grad_p_V;
        
        % % Bounsary conditions

        U(:,1)= 0; V(:,1)=0; % east
        U(:,end)= 0; V(:,end)=0; % west
        U(1,:)= 0; V(1,:)=0; % south
        U(end,:)= 1; V(end,:)=0; % north 
        
        n_step=n_step+1;
    else       % % Adam-Bashford
        U0=U;
        V0=V;
        [U_star, V_star, Ru_old, Rv_old] = midway_velocity_AdamB(U, V, Re, dx, dy, nx, ny, dt, Ru_old, Rv_old);

        % % Calculate the gradient of the pressure to correct the mid-way
        % velocity 
        [grad_p_U,grad_p_V, p,div_velocity_star] = pressure_correction( U_star, V_star, A, nx,ny, dx, dy, dt,tolerance_p,n_steps_p, p);
        grad_p_U=grad_p_U.*dt;
        grad_p_V=grad_p_V.*dt;
        % % Calculate the velocity at step n+1 by correcting it with the
        % previous divergence of p 
        U = U_star - grad_p_U;
        V = V_star - grad_p_V;

        % % Bounsary conditions

        U(:,1)= 0; V(:,1)=0; % east
        U(:,end)= 0; V(:,end)=0; % west
        U(1,:)= 0; V(1,:)=0; % south
        U(end,:)= 1; V(end,:)=0; % north        

        if (mod(n_step,n_step_snst)==0)
          if makevideo==1
            quiver(U(2:end,:),V(:,2:end))
            currframe=getframe;      
            writeVideo(vidobj,currframe);      
          end
    %       Plot and store figures
          filename = ['presion_' num2str(n_step/n_step_snst)];      
          h=dialog ( 'visible', 'off', 'windowstyle', 'normal');      
          ax=axes('parent', h, 'nextplot', 'add' );
          colormap(h,'Jet')
          contourf(ax,p(2:end-1,2:end-1),20);
          hh=colorbar('peer', ax);
          title('Contorno de presión');
          saveas ( ax, filename, 'png' )
          close(h)  
%           filename = ['velocidad_' num2str(n_step/n_step_snst)];      
%           v=dialog ( 'visible', 'off', 'windowstyle', 'normal'); 
%           ax=axes('parent', v, 'nextplot', 'add' );
%           colormap(v,'Jet')
%           quiver(ax, [ U zeros(1,ny+2)'],[zeros(1,nx+2);V])
%     %       contourf(ax,p);
%           title('Campo de velocidad');
%           saveas ( ax,filename, 'png' )
%           close(v)  
          filename = ['velocidadU_' num2str(n_step/n_step_snst)];      
          v=dialog ( 'visible', 'off', 'windowstyle', 'normal'); 
          ax=axes('parent', v, 'nextplot', 'add' );
          colormap(v,'Jet')
%           quiver(ax, U(1:end-1,:),V(:,1:end-1))
          contourf(ax,U,20);
          hh=colorbar('peer', ax);
          title('Campo de velocidad');
          saveas ( ax,filename, 'png' )
          close(v) 
          filename = ['velocidadV_' num2str(n_step/n_step_snst)];      
          v=dialog ( 'visible', 'off', 'windowstyle', 'normal'); 
          ax=axes('parent', v, 'nextplot', 'add' );
          colormap(v,'Jet')
%           quiver(ax, U(1:end-1,:),V(:,1:end-1))
          contourf(ax,V,20);
          hh=colorbar('peer', ax);
          title('Campo de velocidad');
          saveas ( ax,filename, 'png' )
          close(v) 
        elseif errorV<5e-6
            break
        end
                
        n_step=n_step+1
        errorU=norm(U-U0)/norm(U)
        errorV=norm(V-V0)/norm(V)
        p_final=p; 


    end 
    
end
time=toc;
print_matrix(U,V,p_final,errorU,errorV,n_step,time);
if (makevideo==1)
    close(vidobj);
end
    

end


