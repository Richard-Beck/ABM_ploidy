
%%
%Specifying parameters
nx=200;                           %Number of steps in space(x)
ny=200;                           %Number of steps in space(y)       
nt=100;                           %Number of time steps 
dt=0.1;                         %Width of each time step
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
vis=0.01;                         %Diffusion coefficient/viscocity
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)


%Single source initial condition
u(nx/2,ny/2)=10;




%%
%B.C vector
bc=zeros(nx-2,ny-2);
bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs

bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
bc=vis*dt*bc;
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;
%%
%Calculating the field variable for each time step
i=2:nx-1;
j=2:ny-1;

FieldContainer=[];

totalField=sum(sum(u)); % Saves the sum of all values of the matrix at each time step

for it=0:nt
    un=u;
   

         
         %Visualization 
          mesh(x,y,u')
          axis([0 2 0 2 0 2])
          drawnow;
    
    %Uncomment as necessary
    %Implicit method:
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
  
%     
    %Dirichlet:
    u(1,:)=u(2,:);
    u(nx,:)=u(nx-1,:);
    u(:,1)=u(:,2);
    u(:,ny)=u(:,ny-1);

    %Uncomment this to allow for contiuous  source,  turned off at t=ts
    % if it<4
    %  u(nx/2,ny/2)=u(nx/2,ny/2)+dt*8;
    %  end
    % 
    
    
    
  
    FieldContainer(:,:,it+2)=u;
    totalField=[totalField sum(sum(u))];
    dot=[dot u(100,100)];

    
    
end

 plot(totalField);
 xlabel('Time')
 ylabel('Total Field')
 
