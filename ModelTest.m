%% Model testing
clc
clear all
close all



M=1;
tau=.1;
L=1;

A1 = [0 1 0 0 0 0
      0 0 1 0 0 0
      L^3 0 0 L^3 L^2 L
      0 0 0 0 1 0
      0 0 0 0 0 1
      0 0 0 0 0 -1/tau];

A2 = [L 0 0 L 0 0
      0 L 0 0 L 0
      0 0 L 0 0 L
      0 0 0 0 1 0
      0 0 0 0 0 1
      0 0 0 0 0 -1/tau];


A3 = [0 1 0 0 0 0
      0 0 1 0 0 0
      0 0 L 0 0 L
      0 0 0 0 1 0
      0 0 0 0 0 1
      0 0 0 0 0 -1/tau];

  
A4 = [0 1 0 0 0
      L^2 0 L^2 L 0
      0 0 0 1 0
      0 0 0 0 1
      0 0 0 0 -1/tau] ; 
  
B = [0 0 0 0 0 1/(tau*M)]';
B4 = [0 0 0 0 1/(tau*M)]';

H  = [1 0 0 0 0 0 ];

H4 = [1 1 0 0 0];
% Observability
ob1 = rank(obsv(A1,H))
ob2 = rank(obsv(A2,H))
ob3 = rank(obsv(A3,H))
ob4 = rank(obsv(A4,H4))

eg1 = eig(A1);
eg2 = eig(A2);
eg3 = eig(A3);
eg4 = eig(A4);
 
 
 
R = 1;
Q = diag([1e5 0 0 0 0 0]);
N = [1 0 0 0 0 0]';

%% LQR
try
    [K1,S1,CLP1] = lqr(A1,B,Q,R,N);
catch
    Flag1=0;
end

try
    [K2,S2,CLP2] = lqr(A2,B,Q,R,N);
catch
    Flag2=0;
end

try
    [K3,S3,CLP3] = lqr(A3,B,Q,R,N);
catch
    Flag3=0;
end







