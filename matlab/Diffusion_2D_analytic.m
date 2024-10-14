% 2D Diffusion Equation: du/dt = D * (d^2u/dx^2 + d^2u/dy^2)

% Parameters
D = 0.1; % Diffusion coefficient
L = 10;  % Size of the square domain
T = 100;   % Total simulation time
Nx = 100; % Number of spatial grid points in x
Ny = 100; % Number of spatial grid points in y

% Spatial grid
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
[X, Y] = meshgrid(x, y);

% Initial condition
u0 = exp(-((X - L/2).^2 + (Y - L/2).^2)/2);
% u0=zeros(Nx,Ny);
%u0(Nx/2,Ny/2)=10;
%u0(:)=0.25

% Analytic solution
u_analytic = zeros(Nx, Ny, T);
total_field=[];
figure
tol=0.01;

for t = 1:T
    u_analytic(:,:,t) = u0 * exp(-D * ((pi/L)^2 + (pi/L)^2) * t);
%     
%      if t<3
%     u_analytic(Nx/2,Ny/2,t)=u_analytic(Nx/2,Ny/2,t)+10;
%      end
     
    total_field=[total_field sum(u_analytic(:))];
    
    
    %Visualization
         subplot(2,1,1),colormap
          mesh(x,y,u_analytic(:,:,t)) %,pcolor(x,y,u_analytic(:,:,t)), shading interp
          axis([0 L 0 L 0 max(u0(:))]);
          xlabel('Spatial co-ordinate (x) \rightarrow')
          ylabel('{\leftarrow} Spatial co-ordinate (y)')
          zlabel('Total field (u) \rightarrow')
          drawnow;
          
         %Visualization2
         subplot(2,1,2),colormap
          plot(total_field(1:t))
          xlabel('Time')
          ylabel('Total Field')
          drawnow;
          pause(0.001)
            
       
     
end


          


