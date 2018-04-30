%k = 1.500;

t = 1:1:500;
N = 500;
phi = [1.5 1; -0.7 0];
gamma = [1; 0.5];

x = [];
for row=1:1:5
    x(:,1) = [0;0];
    for k=2:1:N
        x(:,k) = (phi*x(:,k-1)) + (gamma*rand);
    end
     plot(t, x(2,:)), xlabel('t'), ylabel('Output x1'), title('5 Sequences')
     hold on
end



