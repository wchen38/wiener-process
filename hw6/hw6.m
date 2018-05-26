clear all
clc 


%Wuyuan Chen
%
init_value = [10000 -500 6*(10^7) 50 0 0 0 200 0 0 0 2*(10^12)];
for k=1:1:10
    %using ode45 to solve for x_dot_hat(t) and A(t) in the prediction step 
    [t, x] = ode45('myfunc', [0 1], init_value); 
    t = t';
    x = x';
    
end


