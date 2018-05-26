clear all
clc 


%Wuyuan Chen
%
x_hat(:, 1) = [10000; -500; 6*(10^7)];
P_matrix(:, 1) = [50; 0; 0; 0; 200; 0; 0; 0; 2*(10^12)];

init_value = [10000; -500; 6*(10^7); 50; 0; 0; 0; 200; 0; 0; 0; 2*(10^12)];

for k=1:1:10
    %-------------------Step 1 Prediction ----------------------%
    %using ode45 to solve for x_dot_hat(t) and A(t) in the prediction step 
    [t, x] = ode45('myfunc', [0 1], init_value); 
    t = t';                             %turn it into a col vector
    x = x';                             %turn it into a col vector
    %getting the last column of values in each ode45 run
    lastX = x(:,end);
    x_hat(:, k) = [lastX(1); lastX(2); lastX(3)];
    P_matrix(:, k) = [lastX(4); lastX(5); lastX(6)
                      lastX(7); lastX(8); lastX(9)
                      lastX(10); lastX(11); lastX(12)];
    
         
    
    
    
   
    
    
    
    
end


