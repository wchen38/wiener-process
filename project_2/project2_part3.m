%Wuyuan Chen 
clear all
clc
S = load('./dataPart3/dataSet1.mat');
% 
 plot(S.xm, S.ym); title('Evolution of Position'); hold on
% figure 
% plot(S.pitch); title('Evolution of Pitch');
% figure 
% plot(S.roll);  title('Evolution of Roll');
% figure 
% plot(S.thrust);  title('Evolution of Thrust');
% figure
% plot(S.yawm);  title('Evolution of Yawm');

pitch = S.pitch * (pi/180);
roll = S.roll * (pi/180);

Q = zeros(3,3);
H = zeros(3, 5);

c=0.1;
g = 9.8;
dt = 0.1;
varQ = [1 1 .001];

ctn = 0;
%------------------inital guess-------------------
xm0 = S.xm(44); ym0 = S.ym(44); vx0 = (S.xm(44)-S.xm(43))/dt;
vy0 = (S.ym(44)-S.ym(43))/dt;
yawm0 = S.yawm(44);

X = [xm0; ym0; vx0; vy0; yawm0];
Xrec = X;
varP0 = [0.01^2 0.01^2 2*(0.01^2/dt^2) 2*(0.01^2/dt^2) (pi/180)^2];
P = diag(varP0);

sigmaVx = 2*(0.01^2)/(dt^2);
sigmaVy = 2*(0.01^2)/(dt^2);
sigmaYaw = (pi/180)^2;

for k=45:length(S.xm)
    Q =  diag(varQ) * dt;
    gamma = [0 0 0;
             0 0 0;
            sigmaVx 0 0;
            0 sigmaVy 0;
            0 0 sigmaYaw];
    df4 = c*g*dt*(cos(X(5))*sin(roll(k-1)) - sin(X(5))*cos(roll(k-1))*sin(pitch(k-1)));
    df5 = c*g*dt*(cos(X(5))*cos(roll(k-1))*sin(pitch(k-1)) + sin(X(5))*sin(roll(k-1)));
    phi = [1 0 dt 0 0;
           0 1 0 dt 0;
           0 0 1 0 df4;
           0 0 0 1 df5;
           0 0 0 0 1 ];
%-------------------Prediction step-------------------------
    Lx = c*g*(sin(X(5))*sin(roll(k-1)) + cos(X(5))*cos(roll(k-1))*sin(pitch(k-1)))*dt;
    Ly = c*g*(sin(X(5))*cos(roll(k-1))*sin(pitch(k-1)) - cos(X(5))*sin(roll(k-1)))*dt;
    Xm = X + [X(3)*dt; X(4)*dt; Lx; Ly; 0];
    Pm = phi*P*phi' + gamma*Q*gamma';
       
%-----------------Update step---------------------------
 res = [S.xm(k); S.ym(k); S.yawm(k)] - [Xm(1); Xm(2); Xm(5)];
 if abs(res(3)) > pi/6
     ctn = ctn +1;
    H = [1 0 0 0 0; 0 1 0 0 0];
    R = diag([0.01^2 0.01^2]);
    K = Pm*H'/(H*Pm*H' + R);
    X = Xm + K*(res(1:2, :));
    P = (eye(5) - K*H)*Pm;
 else
    H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 1];
    R = diag([0.01^2, 0.01^2, (pi/180)^2]);
    K = Pm*H'/(H*Pm*H' + R);
    X = Xm + K*(res);
    P = (eye(5) - K*H)*Pm;
 end
 Xrec = [Xrec X];
end
plot(Xrec(1,:), Xrec(2,:))
