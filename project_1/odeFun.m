function dX=odeFun(t,X)
  
  x=X(1,1);
  dotx=X(2,1);
  beta=X(3,1);
  P(1,1)=X(4,1);
  P(1,2)=X(5,1);
  P(1,3)=X(6,1);
  P(2,2)=X(7,1);
  P(2,3)=X(8,1);
  P(3,3)=X(9,1);
  P(2,1)=P(1,2);
  P(3,1)=P(1,3);
  P(3,2)=P(2,3);

  rho_0 = 1220 ;
  c = 10263 ;
  g = 9.8;
  d = rho_0*exp(-x/c)*dotx^2/(2*beta);
  
  A21 = -rho_0*exp(-x/c)*dotx^2/(2*c*beta);
  A22 =  rho_0*exp(-x/c)*dotx/beta;
  A23 = -rho_0*exp(-x/c)*dotx^2/(2*beta^2); 
  A=[0 1 0 ; A21 A22 A23;0 0 0];

  dX(1,1)=dotx;
  dX(2,1)=d-g;
  dX(3,1)=0;
  
  Gamma=[0; 0; sqrt(1000)];
  dotP=A*P+P*A'+Gamma*Gamma';
   
  dX(4,1)=dotP(1,1);
  dX(5,1)=dotP(1,2);
  dX(6,1)=dotP(1,3);
  dX(7,1)=dotP(2,2);
  dX(8,1)=dotP(2,3);
  dX(9,1)=dotP(3,3);
  
end