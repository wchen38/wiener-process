%Wuyuan Chen
clear all
clc

filename1 = './problem/XYData_cm.csv'
filename2 = './problem/HeadingAngle_rad.csv'
position = csvread(filename1);
headingAngle = csvread(filename2);
xAxis = position(:,1);         %all the x axis data from csv file
yAxis = position(:,2);         %all the y axis data from csv file
dt = 1/3;

%inital guess of X
x0 = xAxis(3);
y0 = yAxis(3);
v0 = (xAxis(2) - xAxis(1)) / dt;
theta0 = headingAngle(1);
X = [x0; y0; v0; theta0];
vRec = X(3);
thetaRec = X(4);

%inital guess of P
varX0 = 1.3^2;
varY0 = 1.3^2;
varV0 = 2*varX0 / dt;
varTheta0 = (5*pi/180)^2;

P = [varX0  0    0    0; 
      0   varY0  0    0;
      0     0  varV0  0;
      0     0    0 varTheta0];
  
P11Rec = P(1,1);
p22Rec = P(2,2);
P33Rec = P(3,3);
P44Rec = P(4,4);

sigmaV = 0.1;
sigmaTheta = 0.2;

gamma = [0 0; 0 0; sigmaV*dt 0; 0 sigmaTheta*dt];
Q = [1 0; 0 1];

%-----------------Prediction step---------------

for k=1:length(xAxis)
    X = X + [vRec*cos(thetaRec)*dt; vRec*sin(thetaRec)*dt; 0; 0];
    
    phi = [1    0   dt*cos(thetaRec) -vRec*dt*sin(thetaRec); 
           0    1   dt*sin(thetaRec)  vRec*dt*cos(thetaRec);
           0    0           1                  1;
           0    0           0                  1;];

    P = phi*P*phi' + (gamma*Q*gamma');
    
     
%-----------Update step -----------------------
H = [1 0 0 0; 0 1 0 0];
R = [1.3^2 0; 0 1.3^2];
Z = [1 0 0 0; 0 1 0 0] * X + [1.3; 1.3];
K = P*H'*inv(H*P*H' + R);
X = X + K*(position(k)' - Z);
xRec(:, k) = X;
P = (eye(4) - K*H)*P;

end
