%Wuyuan Chen
S = load('data3D.mat');

% plot3(S.xm, S.ym, S.zm); title('Evolution of position')
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

sigmaVx = 1;
sigmaVy = 1;
sigmaVz = 1;
sigmaYaw = 1;

 g = 9.8;
 c = 1;

dt0 = S.dtrec(44);
varQ = [1 1 1 .001];
Q = diag(varQ) * dt0;

 


 
%----------------------inital guess ----------------------------%
vy0 = (S.ym(44)-S.ym(43))/dt0; vz0 = (S.zm(44)-S.zm(43))/dt0;
yawm0 = S.yawm(44);
pitch = S.pitch;
roll = S.roll;

Xm = [xm0; ym0; zm0; vx0; vy0; vz0; yawm0];

varP0 = [0.01^2 0.01^2 0.01^2 0.01^2/S.dtrec(44) 0.01^2/S.dtrec(44) 0.01^2/S.dtrec(44) (pi/180)^2];
Pm = diag(varP0);

for k=45:length(S.xm)
     dt = S.dtrec(k);
     gamma = [0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         sigmaVx*sqrt(dt) 0 0 0;
         0 sigmaVy*sqrt(dt) 0 0;
         0 0 sigmaVz*sqrt(dt) 0;
         0 0 0 sigmaYaw*sqrt(dt)];
     
     df4 = c*g*S.dtrec(k-1)*(cos(Xm(7))*sin(roll(k-1)) - sin(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)));
     df5 = c*g*S.dtrec(k-1)*(cos(Xm(7))*cos(roll(k-1))*sin(pitch(k-1)) + sin(Xm(7))*sin(roll(k-1)));
     phi = [1 0 0 dt 0 0 0;
        0 1 0 0 dt 0 0;
        0 0 1 0 0 dt 0;
        0 0 0 1 0 0 df4;
        0 0 0 0 1 0 df5;
        0 0 0 0 0 1 0;
        0 0 0 0 0 0 1;
        ];
    %-------------------Prediction step-------------------------
    Lx = c*g*(sin(Xm(7))*sin(roll(k)) + cos(Xm(7))*cos(roll(k))*sin(pitch(k)))*dt;
    Ly = c*g*(sin(Xm(7))*cos(roll(k))*sin(pitch(k)) - cos(Xm(7))*sin(roll(k)))*dt;
    Lz = (c*g*cos(roll(k))*cos(pitch(k)) - g)*dt;
   
    Xm = Xm + [Xm(4)*dt; Xm(5)*dt; Xm(6)*dt; Lx; Ly; Lz; Xm(7)];
    Pm = phi*Pm*phi' + gamma*Q*gamma';
end 


