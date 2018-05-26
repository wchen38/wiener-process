function xdot=func(t, x)

g_L = 0.1;
xdot(1)=x(2);
xdot(2)= g_L*(sin(x(1)));
xdot=xdot';
