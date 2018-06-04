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
xRec(:, 1) = X;
vRec = X(3);
thetaRec = X(4);

%inital guess of P
varX0 = (1.3^2 +1.3^2+1.3^2)/3;
varY0 = (1.3^2 +1.3^2+1.3^2)/3;
varV0 = varX0 / dt;
varTheta0 = (1*pi/180)^2;

P = [varX0  0    0    0; 
      0   varY0  0    0;
      0     0  varV0  0;
      0     0    0 varTheta0];

P11Rec = P(1,1);
p22Rec = P(2,2);
P33Rec = P(3,3);
P44Rec = P(4,4);

P_Rec(:, 1) = [P11Rec; p22Rec; P33Rec; P44Rec];

sigmaV = 2;
sigmaTheta = 2;

gamma = [0 0; 0 0; sigmaV*sqrt(dt) 0; 0 sigmaTheta*sqrt(dt)];
Q = [1 0; 0 1];

%-----------------Prediction step---------------
index = 2;
for k=1:length(xAxis)
    X = X + [vRec*cos(headingAngle(k))*dt; vRec*sin(headingAngle(k))*dt; 0; 0];
    
    phi = [1    0   dt*cos(headingAngle(k)) -vRec*dt*sin(headingAngle(k)); 
           0    1   dt*sin(headingAngle(k))  vRec*dt*cos(headingAngle(k));
           0    0           1                  0;
           0    0           0                  1;];

    P = phi*P*phi' + (gamma*Q*gamma');
   
    P_Rec(:, index) = [P(1,1); P(2,2); P(3,3); P(4,4);];
    index = index +1;
%-----------Update step -----------------------
H = [1 0 0 0; 0 1 0 0];
R = [1.3^2 0; 0 1.3^2];
Z = [1 0 0 0; 0 1 0 0] * X;
K = (P*H')/(H*P*H' + R);
X = X + K*(position(k,:)' - Z);
xRec(:, k) = X;
P = (eye(4) - K*H)*P;
  P_Rec(:, index) = [P(1,1); P(2,2); P(3,3); P(4,4);]; 
  index = index +1;
end
plot(xAxis, yAxis); hold on
plot(xRec(1,:),xRec(2,:))
figure 
plot(headingAngle); hold on
plot(xRec(4,:));

%--------------------------2nd Order Filter--------------------------------
clear X
X = [x0; y0; v0; theta0];
F_k = [0 0 0 -v0*dt*cos(theta0); 
       0 0 0  v0*dt*sin(theta0);
       0 0 0 0
       0 0 0 0];
   
P_k = [varX0  0    0    0; 
      0   varY0  0    0;
      0     0  varV0  0;
      0     0    0 varTheta0];

G_k = [0 0; 0 0];

% for k = 1:length(xAxis)
%     b_k = b_k + (eye(length(xAxis)*());
% end 
len = length(xAxis);
%syms n
%S1 = symsum(eye(n)*(0.5*trace(F_k*P_k) + 0.5*trace(G_k*Q)), n, [0 4] );
%-----------------------Prediction Step---------------------------
% for k=1:length(xAxis)
%     X = X + [vRec*cos(headingAngle(k))*dt; vRec*sin(headingAngle(k))*dt; 0; 0];
% end 