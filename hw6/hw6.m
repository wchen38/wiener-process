clear all
clc 


%Wuyuan Chen
%
global r1 r2;
r1=1000; r2=500; varR=5;
y_measure = [9055; 8560; 7963; 7467; 7000; 6378; 5885; 5400; 4928; 4503];

x_hat(:, 1) = [10000; -500; 6*(10^7)];
P_vector(:, 1) = [50; 0; 0; 0; 200; 0; 0; 0; 2*(10^12)];
P_matrix = [50 0 0; 0 200 0; 0 0 2*(10^12)];
init_value = [10000; -500; 6*(10^7); 50; 0; 0; 0; 200; 0; 0; 0; 2*(10^12)];
index = 2;
for k=2:1:10
    %-------------------Step 1 Prediction ----------------------%
    %using ode45 to solve for x_dot_hat(t) and A(t) in the prediction step 
    [t, x] = ode45('myfunc', [0 1], init_value); 
    t = t';                             %turn it into a col vector
    x = x';                             %turn it into a col vector
    %getting the last column of values in each ode45 run
    lastX = x(:,end);
    x_hat(:, k) = [lastX(1); lastX(2); lastX(3)];
    x_before_update = x_hat(:, k);
    P_vector(:, index) = [lastX(4); lastX(5); lastX(6)
                      lastX(7); lastX(8); lastX(9)
                      lastX(10); lastX(11); lastX(12)];
    index = index + 1;             
    P_matrix = [lastX(4) lastX(5) lastX(6);
                lastX(7) lastX(8) lastX(9);
                lastX(10) lastX(11) lastX(12)];
    
         
    %-------------------Step 3 Update ----------------------%
%find the partial derivative of h = sqrt(r1^2 + (x_k-r2)^2  )
% [dh/dx;    dh/dv;   dh/dbeta]
dhdx = (x_hat(1)-r2)/sqrt( (x_hat(1)^2)-(r1*x_hat(1)) + 1250000);
H = [dhdx 0 0];
%y_vector = [zkvk_func(x_before_update(1)) zk_func(x_before_update(2)) zk_func(x_before_update(3))];
%obs_h = [zkvk_func(x_before_update(1)) zk_func(x_before_update(2)) zk_func(x_before_update(3))];
y_obs = zkvk_func(x_before_update(1));
h_obs = zkvk_func(x_before_update(1));
%find the kalman gain Version 1
K = P_matrix*H'*inv(H*P_matrix*H' + varR);
x_hat(:, 1) = x_before_update + K*(y_measure(1) - h_obs);
P_matrix = (eye(3) - K*H)*P_matrix;
P_vector(:, index) = [P_matrix(1); P_matrix(4); P_matrix(7);
                   P_matrix(2); P_matrix(5); P_matrix(8);
                   P_matrix(3); P_matrix(6); P_matrix(9)];
index = index + 1;    
init_value = [x_hat(1); x_hat(2); x_hat(3);
    P_matrix(1); P_matrix(4); P_matrix(7); 
    P_matrix(2); P_matrix(5); P_matrix(8);
    P_matrix(3); P_matrix(6); P_matrix(9)];
end

function zk = zkvk_func(x)
global r1 r2
    vk = sqrt(5)*randn;
    zk = sqrt(r1^2 + (x - r2)^2) + vk;
end
function zk = zk_func(x)
global r1 r2
    vk = sqrt(5)*randn;
    zk = sqrt(r1^2 + (x - r2)^2);
end

