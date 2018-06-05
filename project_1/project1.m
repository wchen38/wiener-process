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
theta0 = pi;%headingAngle(1);
X = [x0; y0; v0; theta0];
xRec = X;
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
    X = X + [vRec*cos(thetaRec)*dt; vRec*sin(thetaRec)*dt; 0; 0];
    
    phi = [1    0   dt*cos(thetaRec) -vRec*dt*sin(thetaRec); 
           0    1   dt*sin(thetaRec)  vRec*dt*cos(thetaRec);
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
xRec = [xRec, X];
vRec = X(3);
thetaRec = X(4);
xRec(4, k) = wrapToPi(X(4)); %wrap angles from pi to -pi
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
clear all
filename1 = './problem/XYData_cm.csv'
filename2 = './problem/HeadingAngle_rad.csv'
position = csvread(filename1);
headingAngle = csvread(filename2);
xAxis = position(:,1);         %all the x axis data from csv file
yAxis = position(:,2);         %all the y axis data from csv file
dt = 1/3;


sigmaV = 1.5;
sigmaTheta = 2.1;

x0 = xAxis(3);
y0 = yAxis(3);
v0 = (xAxis(2) - xAxis(1)) / dt;
theta0 = pi;%headingAngle(1);

%inital guess of P
varX0 = (1.3^2 +1.3^2+1.3^2)/9;
varY0 = (1.3^2 +1.3^2+1.3^2)/9;
varV0 = varX0 / dt;
varTheta0 = (2*pi/180)^2;

Xm = [x0; y0; v0; theta0];
XmRec = Xm;
Xhat = Xm;
vRec = Xm(3);
thetaRec = Xm(4);
gamma = [0 0; 0 0; sigmaV*sqrt(dt) 0; 0 sigmaTheta*sqrt(dt)];
Q = [1 0; 0 1];
   
Pm = [varX0  0    0    0; 
      0   varY0  0    0;
      0     0  varV0  0;
      0     0    0 varTheta0];

G_k = [0 0; 0 0];
b_k = 0;
H = [1 0 0 0; 0 1 0 0];
R = [1.3^2 0; 0 1.3^2];

% 2nd partial derivative of df(w, 0)/dx
df1 = [0 0 0 0; 0 0 0 0; 0 0 0 -dt*sin(theta0); 0 0 -dt*sin(theta0) -v0*dt*cos(theta0)];
df2 = [0 0 0 0; 0 0 0 0; 0 0 0 dt*cos(theta0); 0 0 dt*cos(theta0) -v0*dt*sin(theta0)];
df3 = zeros(4);
df4 = zeros(4);

b_k(1,1)=0.5*(trace(df1*Pm));
b_k(2,1)=0.5*(trace(df2*Pm));
b_k(3,1)=0.5*(trace(df3*Pm));
b_k(4,1)=0.5*(trace(df4*Pm));
 

 for k=1:length(xAxis)
    df1 = [0 0 0 0; 0 0 0 0; 0 0 0 -dt*sin(Xm(4)); 0 0 -dt*sin(Xm(4)) -Xm(3)*dt*cos(Xm(4))];
    df2 = [0 0 0 0; 0 0 0 0; 0 0 0 dt*cos(Xm(4)); 0 0 dt*cos(Xm(4)) -Xm(3)*dt*sin(Xm(4))];
    df3 = zeros(4);
    df4 = zeros(4);
    Z = [1 0 0 0; 0 1 0 0] * Xm;
    phi = [1    0   dt*cos(Xm(4)) -Xm(3)*dt*sin(Xm(4)); 
            0    1   dt*sin(Xm(4))  Xm(3)*dt*cos(Xm(4));
            0    0           1                  0;
            0    0           0                  1;];
%-----------------------Prediction Step---------------------------

    Xm = Xm + [Xm(3)*cos(Xm(4))*dt; Xm(3)*sin(Xm(4))*dt; 0; 0] + b_k;
%      b_k = 0.5*(e1*(trace(df1*Pm))...
%             + e2*(trace(df2*Pm))...
%                +e3*(trace(df3*Pm))...
%                 + e4*(trace(df4*Pm)));
    b_k(1,1)=0.5*(trace(df1*Pm));
    b_k(2,1)=0.5*(trace(df2*Pm));
    b_k(3,1)=0.5*(trace(df3*Pm));
    b_k(4,1)=0.5*(trace(df4*Pm));
      
    Pm = phi*Pm*phi' + (gamma*Q*gamma');

%----------------------Update Step--------------------------
     S = 0;
     Pm = Pm -  ((Pm*H')/(H*Pm*H' + R + S) * H*Pm);
     K = (Pm*H')/(R+S)
    % Svec = [0 0; 0 0; 0 0; 0 0];
     Xm = Xm + K*(position(k,:)' - Z)+S;
     XmRec = [XmRec, Xm];
     vRec = Xm(3);
     Xhat = Xm;
     thetaRec = Xm(4);
 end 
 
 figure
 plot(xAxis, yAxis); hold on
 plot(XmRec(1,:), XmRec(2,:));
 legend('Measure','Predict','Location','southwest')
 