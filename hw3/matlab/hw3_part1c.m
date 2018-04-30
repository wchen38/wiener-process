%k = 1.500;

t = 1:1:50;
N = 50;
phi = [1.5 1; -0.7 0];
gamma = [1; 0.5];

x = []; stdx1 = []; stdx2 = [];
for row=1:1:N
    x(:,1) = [0;0];
    for k=2:1:N
        x(:,k) = (phi*x(:,k-1)) + (gamma*randn);
    end
    hold on
    plot(t, x(2,:));
    stdx1(row,:) = x(1,:);
    stdx2(row,:) = x(2,:);
end
std_x1 = std(stdx1(:,2))
std_x2 = std(stdx2(:,2))


