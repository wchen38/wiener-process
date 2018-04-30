phi = [1.5 1; -0.7 0];
gamma = [1; 0.5];
t = [0:0.05 10];
P = [0 0; 0 0];
oldP = [0 0; 0 0];
for k=2:1:50
    P = (phi*oldP*phi') + (gamma*gamma');
    oldP = P;
    
    plot(t,P(1), 'o');
    hold on
end
display(P);