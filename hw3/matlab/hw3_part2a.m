expectW= [];
varW = [];
delta_t = 0.05;
t1 = 1;
t2 = 2;
t3 = -1;
%when t equal to 2, the index will be 22
index = 21;
x(1) = 1;
x(index) = 2;
n = 2-0.05;
%wiener process of w1 and w2 is just a normal distribution
w(1) = rand*sqrt(1);
w(index)= rand*sqrt(2-1);

%expected value of given variables
expectW(1) = 0;
expectW(index) = 0;

%variance of given variables
varW(1) = 1;
varW(index) = 1;

%calcuates the expected value from t1 to t2
%each expected value is stored in expectW(k)

    k=1;
    for t=1.05:0.05:n
        k = k+1;
        t_prev=t-delta_t;
        A = ( (t2-t)/(t2-t_prev) )*w(k);
        B = ( (t-t_prev)/(t2-t_prev) )*w(index);
        expectW(k) = A+B;
        
    end
    
   % plot(expectW);
   % hold on

%calcuates the variance from t1 to t2
%each variance value is stored in varW(k)
    k=1;
    for t=1.05:0.05:n
        k = k+1;
        t_prev=t-delta_t;
        A = (t2-t) * (t2-t_prev);
        B = t2-t_prev;
        varW(k) = A/B;
        
    end
    %plot(varW);
    
%calculate x(k)

k=1;
for t=2:1:index-1
    rand_val = expectW(k) + (rand*sqrt(varW(k)));
    x(t) = x(t-1) + rand_val;
    k =k+1;
end

plot(x);


