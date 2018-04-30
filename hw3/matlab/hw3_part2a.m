expectW= [];
varW = [];
delta_t = 0.05*2;
t1 = 1;
t2 = 2;
t3 = 3;
%when t equal to 2, the index will be 22
index = 21;
t3_index = 41;
x(1) = 1;
x(index) = 2;
x(t3_index) = -1;
n = t2-0.05;
m = t3-0.05;


%expected value of given variables
expectW(1) = 0;
expectW(index) = 0;
expectW(t3_index) = 0;
%variance of given variables
varW(1) = 1;
varW(index) = 1;
varW(t3_index) = 1;


%run the sequence 100 times
for b=1:1:100
    k=1;
    %calcuates the expected value from t1 to t2
    %each expected value is stored in expectW(k)
    for t=1.05:0.05:n
        t_prev=t-delta_t;
        
        %expected value
        A = ( (t2-t)/(t2-t_prev) )*x(k);
        B = ( (t-t_prev)/(t2-t_prev) )*x(index);
        expectW(k+1) = A+B;
        
        %variance
        A = (t2-t) * (t-t_prev);
        B = t2-t_prev;
        varW(k+1) = A/B;
        
        %x(k)
        x(k+1) = (rand*varW(k+1)) + expectW(k+1);
        k = k+1;
        
    end
    
    %calcuates the expected value from t2 to t3
    %each expected value is stored in expectW(k)
    k=index;
    for t=2.05:0.05:m
        
        t_prev=t-delta_t;
        
        %expected value
        A = ( (t3-t)/(t3-t_prev) )*x(k);
        B = ( (t-t_prev)/(t3-t_prev) )*x(t3_index);
        expectW(k+1) = A+B;
        
        %variance
        A = (t3-t) * (t-t_prev);
        B = t3-t_prev;
        varW(k+1) = A/B;
        
        %x(k)
        x(k+1) = (rand*varW(k+1)) + expectW(k+1);
        k = k+1;
        
    end
    plot(x);
    hold on
    
end




