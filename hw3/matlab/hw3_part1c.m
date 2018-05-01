%k = 1.500;

t = 1:1:50;
N = 50;
phi = [1.5 1; -0.7 0];
gamma = [1; 0.5];

x = []; stdx1 = []; stdx2 = [];
for row=1:1:N
    x(:,1) = [0;0];
    for k=2:1:500
        x(:,k) = (phi*x(:,k-1)) + (gamma*randn);
    end
    %storing the x1 and x2 in seperate vector
    stdx1(row,:) = x(1,:);
    stdx2(row,:) = x(2,:);

end

%plotting the std of all the columns
for k=1:1:500
    std_x1(k) = std(stdx1(:,k));
    std_x2(k) = std(stdx2(:,k));
end



plot(std_x2), xlabel('t'), ylabel('std of x1'), title('500 Sequences')
mean(std_x2)
