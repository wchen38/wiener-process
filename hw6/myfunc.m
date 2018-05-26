function xdot = myfunc(t, x)

c = 10263;
p_0 = 1220;
rho = p_0*(exp(-x(1)/c));
g = 9.8;
d = (rho*(x(2)^2))/ (2*x(3));


%loading the inital condition of covariance matrix P
%since the inital variance guess is independent from each other
%then the other covariance are zero.
xdot = zeros(9,1);
P(1,1)=x(4);        %variance of height
P(1,2)=x(5);
P(1,3)=x(6);
P(2,1) = x(7);
P(2,2)=x(8);        %variance of velocity
P(2,3)=x(9);
P(3,1)=x(10);
p(3,2) = x(11);
P(3,3)=x(12);       %variance of beta



%finding A(t) which is the jacobian
%dv/dx                 dv/dv                dv/dbeta
%d(d-g)/dx             d(d-g)/dv            d(d-g)/dbeta
%d(sqrt(1000))/dx       d(sqrt(1000))/dv     d(sqrt(1000))/dbeta 
%first row
A11 = 0;
A12 = 1;
A13 = 0;
%second row
A21 = -(610*exp(-x(1)/10263)*(x(2)^2)) / (10263*x(3));
A22 = (1220*exp(-x(1)/10263)*x(2) )/x(3);
A23 = -(610*exp(-x(1)/10263)*(x(2)^2))/(x(3)^2);

%third row
A31 = 0;
A32 = 0;
A33 = 0;

A = [A11 A12 A13; A21 A22 A23; A31 A32 A33];

b = [0 0 0;0 0 0;0 0 1000];

%calculativing the coveriance matrix
Pdot=A*P+P*A'+b;

%storing everything back to the return variable.
xdot(1) = x(2) ;    %height
xdot(2) = d - g;    %velociy
xdot(3) = 0;        %beta

xdot(4)=Pdot(1,1);
xdot(5)=Pdot(1,2);
xdot(6)=Pdot(1,3);
xdot(7)=Pdot(2,1);
xdot(8)=Pdot(2,2);
xdot(9)=Pdot(2,3);
xdot(10)=Pdot(3,1);
xdot(11)=Pdot(3,2);
xdot(12)=Pdot(3,3);
