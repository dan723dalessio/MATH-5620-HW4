%x is the spatial vairable
%x_0 is the inial condition of the spatial variable 
%t is the time varialbe
%t_0 is the initial conidiotn of the time variable
%h is the spatial step
%k is the time step
%D is some heating constant
%g_0 is initial temperature at time t=0
%h_0 is some temerature at spacial location x=0
%h_1 is some temperature at spatial loation at x=1
%A is the coefficient matrix of the form Ax=b
%m is the number of spatial steps and time steps (assume spatial steps = time steps)

%the goal of this exercise to the find numerical approximaitons to the 1D
%heat equation with the given funciton 
function OneD = OneDBackwardEuler(x_0,t_0,h,k,D,g_0,g_1,m)
a = k/(h^2); %this reduces the complexity of our matrix

    %declare the function for the problem 
    function funct = f(x,t) 
        f = cos(t)*cos(2*pi*x)+4*pi*pi*D*sin(t)*cos(2*pi*x); %the function we've define for this problem
    end

A=zeros(m,m); %initialize the A matrix

%=====================================================================
%the following linees are used to construct my A matrix
A(1,1)=1+2*D*a; %define the first element of our matrix
A(1,2)=-D*a; %define the second element of our matrix
for i = 1:m %fill in the rows of the A matrix
    for j=2:m-1 %fill in the columns of the A matrix
        if i==j
            A(i,j)=2*D*a;
            A(i,j-1)=-D*a;
            A(i,j+1)=-D*a;
        end 
    end
end
A(m,m-1)=1+2*D*a; %define the second to last element of the array 
A(m,m)=-D*a; %define the last element of the array
%======================================================================
A_inv=inv(A); %define the inverse matrix

%the following lines are used to construct the matrix b
%here, we construct our solutions based off of the elements in the
%matrix A that we've constructed

b=zeros(1,m); %create a matrix for the b vector 
U_j=zeros(1,m); %create a matrix for the U_j vector
U_i_j=zeros(m,m); %create a matirx for the whole ``big" U vector
ut_0=zeros(1,m); %ut_0 is a time vector with the inital conndition for t=0
%this was given in the problem description 
ux_0=zeros(1,m); %ux_0 is a spatial vector with the initial condition for x=0
%fill in ux_0


U_j(1,1)=@(x,t) f(x_0,t_0);

for j = 1:m %iterate over all time steps 
    for i=2:(m-1) %iterate over all spatial steps
        ux_0(1,j) = sin(j); %this creates a vector of the initial conditions at each time step 
        U_j(1,i) = ut_0 + k*@(x,t) f(i,j) + D*a*ux_0;
        
    end 
end 



end