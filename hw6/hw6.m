clear all
clc 


%Wuyuan Chen
%
init_value = [10000 -500 6*(10^7) 50 0 0 0 200 0 0 0 2*(10^12)];
%init_value = [10000; -500; 6*(10^7); 50; 0; 0; 200; 0; 2*(10^12)];
%init_value = [1 1];
[t, x] = ode45('myfunc', [0 1], init_value);
t = t';
x = x';

