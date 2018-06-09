%Wuyuan Chen
clear all
clc
S = load('data3D.mat');

 plot3(S.xm, S.ym, S.zm); title('Evolution of position'); hold on
% figure
% plot(S.pitch);title('Evolution of Pitch')
% figure
% plot(S.roll);title('Evolution of Roll')
% figure
% plot(S.yawm);title('Evolution of yam')
% figure 
% plot(S.thrust); title('Evolution of thrust')
% figure 
% plot(S.ans);title('Evolution of angles')
% figure 
% plot(S.xm); title('Evolution of x')
% figure
% plot(S.ym); title('Evolution of y')
% figure
% plot(S.zm); title('Evolution of z')

len = length(S.xm);
sigmaVx = 2*0.01^2/(S.dtrec(44)^2);
sigmaVy = 2*0.01^2/(S.dtrec(44)^2);
sigmaVz = 2*0.01^2/(S.dtrec(44)^2);
sigmaYaw = (pi/180)^2;

sumRes= 0;
sumResRec = [];

 g = 9.8;
 %c = 1;

dt0 = S.dtrec(44);
varQ = [1 1 1 .001];


varR = [0.01^2, 0.01^2, 0.01^2, (pi/180)^2];
R = diag(varR);

 
%----------------------inital guess ----------------------------%
xm0 = S.xm(44); ym0 = S.ym(44); zm0 = S.zm(44); vx0 = (S.xm(44)-S.xm(43))/dt0;
vy0 = (S.ym(44)-S.ym(43))/dt0; vz0 = (S.zm(44)-S.zm(43))/dt0;
yawm0 = S.yawm(44);
pitch = S.pitch * (pi/180);
roll = S.roll * (pi/180);

Xm = [xm0; ym0; zm0; vx0; vy0; vz0; yawm0];
xRec = Xm;
varP0 = [0.01^2 0.01^2 0.01^2 2*0.01^2/(S.dtrec(44)^2) 2*0.01^2/S.dtrec(44)^2 2*0.01^2/S.dtrec(44)^2 (pi/180)^2];
Pm = diag(varP0);

for c0=2:0.01:3
    xRec = Xm;
    for k=45:length(S.xm)
        dt = S.dtrec(k-1);
        s = (255/6000) * S.thrust(k);
        c = (0.000409*(s^2) + 0.1405*s - 0.099)/c0;
        Q = diag(varQ) * dt;
        gamma = [0 0 0 0;
            0 0 0 0;
            0 0 0 0;
            sigmaVx 0 0 0;
            0 sigmaVy 0 0;
            0 0 sigmaVz 0;
            0 0 0 sigmaYaw];
        
        H = [1 0 0 0 0 0 0;
            0 1 0 0 0 0 0;
            0 0 1 0 0 0 0;
            0 0 0 0 0 0 1];
        % df4 = 0;
        %df5 = 0;
        df4 = c*g*S.dtrec(k-1)*(cos(Xm(7))*sin(roll(k-1)) - sin(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)));
        df5 = c*g*S.dtrec(k-1)*(cos(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)) + sin(Xm(7))*sin(roll(k-1)));
        phi = [1 0 0 dt 0 0 0;
            0 1 0 0 dt 0 0;
            0 0 1 0 0 dt 0;
            0 0 0 1 0 0 df4;
            0 0 0 0 1 0 df5;
            0 0 0 0 0 1 0;
            0 0 0 0 0 0 1 ];
        %-------------------Prediction step-------------------------
        Lx = c*g*(sin(Xm(7))*sin(roll(k-1)) + cos(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)))*dt;
        %Lx = 0;
        Ly = c*g*(sin(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)) - cos(Xm(7))*sin(roll(k-1)))*dt;
        %Ly = 0;
        Lz = (c*g*cos(roll(k-1))*cos(pitch(k-1)) - g)*dt;
        %Lz = 0;
        
        Xm = Xm + [Xm(4)*dt; Xm(5)*dt; Xm(6)*dt; Lx; Ly; Lz; 0];
        Pm = phi*Pm*phi' + gamma*Q*gamma';
        
        %-----------------Update step----------------------
        
        K = Pm*H'/(H*Pm*H' + R);
        res = [S.xm(k); S.ym(k); S.zm(k); S.yawm(k)] - [Xm(1); Xm(2); Xm(3); Xm(7)];
        Xm = Xm + K*(res);
        xRec = [xRec Xm];
        Pm = (eye(7) - K*H)*Pm;
        
        %sum of the residual
        sumRes = sumRes+(res(1)^2 + res(2)^2 +res(3)^2 +res(4)^2);
    end
    plot3(xRec(1,:), xRec(2,:), xRec(3,:)); hold on
    sumResRec = [sumResRec (1/(len-44))*sumRes];
    clear xRec
end 

