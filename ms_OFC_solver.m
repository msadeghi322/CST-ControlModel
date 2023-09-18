function sout = ms_OFC_solver(simdata,C,Ke)
Flag = 1; % 1: solve ofc, 0: use already solved ofc
if nargin>1
    Flag=0;
end
% sout = minmaxfc_pointMass(xinit,xfinal,simdata)
%
% Calculates a trajectory with initial condition, final target and
% parameters defined in the input structure simdata.
%   XINIT: Initial State
%   XFINAL: Target State
%   Input data structure must contain:
%   SIMDATA.delta       Discretization step
%          .delay       Hard temporal delay in the closed loop system
%          .pert        1x2 vector with step force magnitude along x and y axes
%          .time        Time horizon
%          .gamma       1x2 with Parameter for optimal disturbance
%                       rejection level. The second entry (1 or 0) indicates
%                       whether the routine should optimize this value.
%          .nStep       Number of time steps
%          .noise       1x2 vector of scaling parameters for noise matrices
%                       (default: [1 1])
%          .ralpha      matrix with one row per state variable and one column per time
%                       step with the cost of the corresponding state and time
%          .nsimu       number of simulation runs.
%
%
%   SOUT: output data structure with the following fields:
%   sout.L              Series of optimal robust control gains
%       .C              Series of optimal LQG control gains
%       .x              State - Robust control
%       .xest           State Estimate - Robust control
%       .z              State - LQG
%       .zest           State Estimate - LQG
%       .u              Series of Control Vector - Robust control
%       .v              Series of Control Vector - LQG
%       .minlambda      Minimum eigen value (optimized or used). Must be > 0.
%       .cost           1x2 vector with movement cost (1: Robust, 2: LQG)
%       .gammaopt       Optimal or used gamma parameter
%
%
%
%   Uses: > AugRobustControl
%         > extLQG
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances.
%   Crevecoeur F., Scott S. H., Cluff T.
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

delta = simdata.delta;
delay = simdata.delay;
ralpha = simdata.ralpha;
nStep = simdata.nStep;
statedim = size(simdata.xinit,1);

% Mapping all final targets on 0 to ensure positive definiteness of the
% cost matrices
xinit = simdata.xinit-simdata.xfinal;

% System Matrices
A = simdata.A;
B = simdata.B;
Aest = A;
DA = (A-Aest)*delta; % Used when there is a model error
A = eye(size(A))+delta*A;
Aest = eye(size(Aest))+delta*Aest;
B = delta*B;

% Observability Matrix
H = simdata.H;

% Definition of the cost matrices:
Q = zeros(size(A,1),size(A,2),nStep);
Id = eye(statedim);

%Filling in the cost matrices
for j = 1:nStep
    for i = 1:statedim
        
        Q(:,:,j) = Q(:,:,j) + ralpha(i,j)*Id(:,i)*Id(:,i)';
        
    end
end



% Augment the System to Include the Feedback Delay
A0 = A;
DA0 = DA;
Aest0 = Aest;
B0 = B;
Q0 = Q;
H0 = H;
[A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta);
[Aest,~,~,~,~] = AugRobustControl(Aest0,DA0,B0,Q0,H0,delay,delta);


%Signal Dependent Noise
nc = size(B,2);
Csdn = zeros(size(B,1),nc,nc);
for i = 1:nc
    Csdn(:,i,i) = simdata.noise(2)*B(:,i);
end


D = zeros(size(A));
%D(1:8,1:8) = eye(8);

%--------------------------------------------------------------------------
statedim = size(A,1);

%Forward Simulation of the System Trajectory
h = max(0,floor(delay/delta))+1;

% Parallel Simulation for LQG control
currentZ = kron(ones(h,1),xinit);
currentZEst = currentZ;
z = zeros(nStep,statedim); z(1,:) = currentZ;
zest = z;
v = zeros(nStep-1,size(B,2));
Oxi = simdata.noise(1)*B*B';
%Omega = eye(size(H0,1))*simdata.noise(2);
Omega = simdata.Omega; % eye(size(H0,1))*simdata.noise(2);

%--------------------------------------------------------------------------
% Extended LQG
RLQG = zeros(size(B0,2),size(B0,2),nStep-1);
%RLQG = zeros(2,2,nStep-1);
for i = 1:nStep-1
    RLQG(:,:,i) = simdata.effort*eye(size(B0,2));
end
Cstate = eye(statedim)*10^-2;

%Q(1:8:end,1:8:end,:)=0;
%Q(3:8:end,3:8:end,:)=Q(3:8:end,3:8:end,:)/1300;
%Q(1,1,:)=0*Q(1,1,:);
if Flag
[C,Ke,~,~,~,~,~] = extLQG(Aest,B,Csdn,0*H,H,Q,RLQG,Oxi,Omega,0*A,Cstate,Cstate,currentZ);
end
%--------------------------------------------------------------------------

% Compute the total cost
cost = zeros(1,1);

for i = 1:nStep-1
    
    sensoryNoise = mvnrnd(zeros(size(Omega,1),1),Omega)';
    motorNoise   = mvnrnd(zeros(size(Oxi,1),1),Oxi)';
    
    
    %LQG CONTROL ----------------------------------------------------------
    yz = H*currentZ + sensoryNoise;
    v(i,:) = (-C(:,:,i)*currentZEst)';
    K = Ke(:,:,i);
    currentZEst = Aest*currentZEst + B*v(i,:)' + K*(yz-H*currentZEst);
    
    % Cost LQG
    cost(1) = cost(1) + currentZ'*Q(:,:,i)*currentZ + v(i,:)*v(i,:)';
    
    %Signal Dependent Noise - LQG
    sdn = 0;
    for isdn = 1:nc
        sdn = sdn + normrnd(0,1)*Csdn(:,:,isdn)*v(i,:)';
    end
    
    wz = DA*currentZ;
    currentZ = Aest*currentZ + B*v(i,:)' + D*wz + motorNoise + sdn;
    z(i+1,:) = currentZ(1:statedim)';
    zest(i+1,:) = currentZEst(1:statedim)';
    
end

% Add the final cost
cost(1) = cost(1) + currentZ'*Q(:,:,end)*currentZ;

% Output structure
sout.C = C;                      % Series of optimal LQG control gains
sout.z = z(:,1:length(H0));               % State - LQG
sout.zest = zest(:,1:length(H0));         % State Estimate - LQG
sout.v = v;                      % Series of Control Vector - LQG
sout.cost = cost;                % Movement cost (1: Robust, 2: LQG)
sout.K    = Ke;

end




% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





function [L,K,newCost,Sx,Se,s,JC] = extLQG(A,B,C,D,H,Q,R,cxi,comega,ceta,ce,cxhat,xInit)

% [L,K,newCost,Sx,Se,s,sall,JC] = extLQG(A,B,C,D,H,Q,R,cxi,comega,ceta,ce,cxhat,cxe,xInit)
%
% Calculates a series of feedback gains and Kalman gains for extended LQG
% control. The routine follows (Todorov, 2005, Neural Comp, 17, 1084-1108).
%   [A, B, C, D, H]     matrices for the state pace representation
%   [Q, R]              cost matrices
%   CXI, COMEGA, CETA   motor, sensory and internal noise
%                       covariance matrices
%   CE, CXHAT           initial covariance matrices for the error
%                       and the estimation
%   XINIT               initial state vector
%
% The output matrices are:
%   L                   time series of optimal feedback gains
%   K                   time series of non adaptive Kalmen gains
%   NEWCOST             expected cost after stopping the iterations
%   SX, SE              series of cost matrices determined by the backwards
%                       recurrences for the state (SX) and error (SE) terms
%   S                   scalar component of the total expected cost
%   JC                  covariance matrix
%
%
%   Uses: > computeOFC
%         > computeExtKalman
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances.
%   Crevecoeur F., Scott S. H., Cluff T.
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019


n = size(A,1);
m = size(B,2);
p = size(H,1);
c = size(C,3);
step = size(R,3);

K = zeros(n,p,step);

tol = 10^-14;
current = 10^6;
itmax = 100;
count = 0;
found = false;


% The optimal control and Kalman gains are calculated iteratively (no more
% than itmax times if it does not converge)
while ~found && count < itmax
    
    [L,Sx,Se,s] = computeOFC(A,B,C,D,H,Q,R,K,cxi,comega,ceta);
    [K,JC] = computeExtKalman(A,B,C,D,H,cxi,comega,ceta,L,cxhat,ce);
    
    % Expected cost
    newCost = xInit'*Sx*xInit + trace((Sx+Se)*ce)+ s;
    
    % Relative improvement
    dCost = abs(current - newCost)/newCost;
    %fprintf('dCost = %.14f \n',(current - newCost)/newCost)
    current = newCost;
    
    % Check: of the relative improvement is small, the solution is found
    if dCost > tol
        found = false;
    else
        found = true;
    end
    
    count = count + 1;
    
    
end

fprintf('Number of iterations (ext LQG): %d\n', count); % Number of iterations

end




% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




function [A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta)

% [A,DA,B,Q,H] = AugRobustControl(A0,DA0,B0,Q0,H0,delay,delta)
%
% Augment the system matrices to take the feedback delay into account.
%   [A0, DA0, B0, Q0, H0] are the state space representation matrices
%   without delay. The output matrices include the delay.
%
%   DELAY: hard time shift in the closed loop system in ms.
%   DELTA: discretization step.
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances.
%   Crevecoeur F., Scott S. H., Cluff T.
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

h = floor(delay/delta); %Feedback delay in number of sample times

n = size(A0,1);
m = size(B0,2);
t = size(Q0,3);
p = size(H0,1);

A  = zeros((h+1)*n,(h+1)*n);
DA = zeros((h+1)*n,(h+1)*n);
B  = zeros((h+1)*n,m);
Q  = zeros((h+1)*n,(h+1)*n,t);
H  = zeros(p,(h+1)*n);

A(1:n,1:n) = A0;
A(n+1:end,1:end-n) = eye(h*n);
DA(1:n,1:n) = DA0;
B(1:n,:) = B0;
H(:,end-n+1:end) = H0;

% Adding h times the constraint Q1:
Qaug = zeros(n,n,t+h);
for i = 1:h
    Qaug(:,:,i) = Q0(:,:,1);
end
for i = 1:t
    Qaug(:,:,i+h) = Q0(:,:,i);
end

%Filling the diagonal Q matrices
for time = 1:t
    for ii = 0:h
        Q(ii*n+1:(ii+1)*n,ii*n+1:(ii+1)*n,time) = Qaug(:,:,time+h-ii)/(h+1);
    end
end

end




% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



function [L,Sx,Se,s] = computeOFC(A,B,C,D,H,Q,R,K,oXi,oOmega,oEta)

%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances.
%   Crevecoeur F., Scott S. H., Cluff T.
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

n = size(A,1);
m = size(B,2);
c = size(C,3);
d = size(D,3);
step = size(R,3);
L = zeros(m,n,step);

currSx = Q(:,:,end);
currSe = 0;
currs = 0;

for i = step:-1:1
    
    sdn = 0;
    
    for j = 1:c
        
        sdn = sdn + C(:,:,j)'*(currSx + currSe)*C(:,:,j);
        
    end
    
    statedn = 0;
    
    for j = 1:d
        
        statedn = statedn + D(:,:,j)'*K(:,:,i)'*currSe*K(:,:,i)*D(:,:,j);
        
    end
    
    L(:,:,i) = (R(:,:,step) + B'*currSx*B + sdn)\(B'*currSx*A);
    currSxTemp = currSx;
    currSx = Q(:,:,i) + A'*currSx*(A-B*L(:,:,i)) + statedn;
    currSeTemp = currSe;
    currSe = A'*currSxTemp*B*L(:,:,i)+...
        (A-K(:,:,i)*H)'*currSeTemp*(A-K(:,:,i)*H);
    currs = trace(currSxTemp*oXi+currSeTemp*(...
        oXi+oEta+K(:,:,i)*oOmega*K(:,:,i)'))+currs;
    
end

Sx = currSx;
Se = currSe;
s = currs;

end





% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




function [K,JC] = computeExtKalman(A,B,C,D,H,oXi,oOmega,oEta,L,cxhat,ce)

%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances.
%   Crevecoeur F., Scott S. H., Cluff T.
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019

n = size(A,1);
k = size(H,1);
d = size(D,3);
c = size(C,3);
step = size(L,3);
sigmaE = ce;
sigmaX = cxhat;
sigmaEX = zeros(n);

K = zeros(n,k,step);
JC = zeros(step,3);

for i = 1:step
    
    statedn = 0;
    sTemp = (sigmaE+sigmaX+sigmaEX+sigmaEX');
    
    for j = 1:d
        
        statedn = statedn + D(:,:,j)*sTemp*D(:,:,j)';
        
    end
    
    sdn = 0;
    
    for j = 1:c
        
        sdn = sdn + C(:,:,j)*L(:,:,i)*sigmaX*L(:,:,i)'*C(:,:,j)';
        
    end
    
    K(:,:,i) = (A*sigmaE*H') / (H*sigmaE*H'+oOmega+statedn);
    
    
    
    sigmaETemp = sigmaE;
    sigmaE = oXi+oEta+(A-K(:,:,i)*H)*sigmaE*A'+sdn;
    term = (A-B*L(:,:,i))*sigmaEX*H'*K(:,:,i)';
    sigmaX = oEta + K(:,:,i)*H*sigmaETemp*A'+(A-B*L(:,:,i))*sigmaX*(A-B*L(:,:,i))'...
        + term + term';
    sigmaEX = (A-B*L(:,:,i))*sigmaEX*(A-K(:,:,i)*H)'-oEta;
    
    JC(i,:) = [sigmaE(1,1),sTemp(1,1),0];
    
end

end














































