%% Simulating the CST based on a 5-State model


%% Simulating experiments with OFC
clc
clear
close all

% There are three different A matrices in the model_1 function. 


FileName = 'sim';
FilePath = 'Reduced data/CST_R1.2_test5';
PerturbFlag = 0; % random initial condition for cursor on or off

%% Fixed parameters
Lambda_List = (.5:.4:6.5)';
Obs_List = [1 1 1 1 1];
%Q_List = [0   1e10 0 0 0 ]; % Let's fix things on velocity control
Q_List = [1e5 0 0 0 0 
          0   1e10 0 0 0];
      
      
%% Independent parameters 
DelayList      = [20, 50, 100];
EffortCostList = [10, 100 , 1000];

%% Tunning parameters:
% These parameters are tuned to result in a fixed success rate when
% changing the independent parameters
MotorNoiseList = [.1  2.35
                  .1  2.05];
SensoryNoise   = diag([1e-6 , 1e-6 , 1e-6 , 1e-6 , 1e-6]);              
simdata.tau    = .07;
simdata.mass   = 1;
simdata.delta  = .01;        % Discretization step: 10ms
simdata.delay  = .05;        % feedback loop delay, xx time steps
simdata.time   = 8;         % Reach time
simdata.nStep  = simdata.time/simdata.delta+1;         % Number of time steps corresponding to reach time (600ms), plus terminal step
Time           = (0:simdata.delta:simdata.time);
Trials         = 200; % Number of simulation runs.

% feedback type (Matrix H): cursor(c)/hand(h); pos(p)/vel(v)/acc(a): Assuming only one
% observation
Obs = licols( diag(Obs_List(1,:)) )';
Omega = licols(Obs*SensoryNoise);
simdata.Omega = Omega;



c=1;
for i=1:size(Q_List,1)
    % Q matrix: which state should be penalized more
    qq = Q_List(i,:)';
    simdata.noise  = MotorNoiseList(i,:);   % Motor, and Signal dependent noise, standard values.

    for j=1:size(EffortCostList,2)
        simdata.effort = EffortCostList(j);
    
        for k=1:length(Lambda_List)
            L = Lambda_List(k);
            simdata.Lambda = L;
            % Initial conditions
            hp0  = 0;
            hv0  = 0;
            cp0  = 0.00;
            cv0  = L*(cp0+hp0);
            F0   = 0;
            Init = [cp0,cv0,hp0,hv0,F0]';
            Finl = [0 0 0 0 0]';

            clear Sim
            flp = sprintf('%s/Q%d_H%d',FilePath,i,j);
            fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
            if ~exist(flp,'file')
               mkdir(flp) 
            end
            
            if exist(fln,'file')
                fprintf('%s already exists \n',fln);
                c=c+Trials;
                continue
            end
            
            Flag1 = 1; % solve
            
            % Solve the OFC and check if solution is valid (converged). If
            % not, perturb the Q matrix a bit and try again.
            Flag_NumericalError = 1; % check if the solution is valid
            while Flag_NumericalError==1
                Sim1 = model_1(Flag1,simdata,Init,Finl,Obs,qq,[],[]);
                if max(abs(Sim1.C_p))>1e4 && L<8
                    qq = qq*( 1+ .1*randn(1) );
                else
                    Flag_NumericalError=0;
                end
            end
            
            C = Sim1.C;
            Ke = Sim1.K;
            for n=1:Trials
                tic
                Flag2 = 0;
                X0_rand = Init + PerturbFlag*[.01*sign(randn(1)) , 0,0,0,0]';
                Sim(n) = model_1(Flag2,simdata,X0_rand,Finl,Obs,qq,C,Ke);
                if n>1 % we only need the gains for first trial, the rest is repeated
                    Sim(n).C=[];
                    Sim(n).K=[];
                end
                
                fprintf('Trial %d \n',n)
                c=c+1;
                
            end
            Sim(1).SimulationMetaData = simdata;
            Sim(1).LambdaIncrements = Lambda_List;
            Sim(1).ControlPolicyList = Q_List;
            Sim(1).SensoryNoise = SensoryNoise;
            Sim(1).MotorNoiseList = MotorNoiseList;
            Sim(1).ObsMatrix = Obs_List;
            save(fln,'Sim','-v7.3');
            fprintf('%s is saved \n',fln);
            
        end
    end
end



%%
D   = load(fln);
Sim = D.Sim;
L = simdata.Lambda;
for i = 1:Trials
    
    C_p = Sim(i).C_p;
    C_v = Sim(i).C_v;
    H_p = Sim(i).H_p;
    H_v = Sim(i).H_v;
    H_f = Sim(i).H_f;
    
    if i>10
       continue 
    end
    
    
    figure(i)
    sb1 = 3;
    sb2 = 2;
    clf
    % Puts hold on to add simulation to the figure
    subplot(sb1,sb2,1)
    hold all
    plot([0,Time(end)],[0,0],':k')
    H1 = plot(Time,C_p,'b','linewidth',2);
    H2 = plot(Time,H_p,'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    legend([H1,H2],'Cursor','Hand','location','best')
    ylim([-.05 .05])
    
    subplot(sb1,sb2,2)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], max(abs(C_v))*[-1,1] ,':k');
    plot(C_p,C_v,'color',[.6,.6,.6])
    hh = scatter(C_p,C_v,50,copper(length(C_v))); hh.Marker = '.';
    xlabel('Cursor P')
    ylabel('Cursor V')
    xlim([-.05 .05])
    
    
    subplot(sb1,sb2,3)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], .05*[-1,1] ,'-k','linewidth',1);
    plot(.05*[1,-1] , .05*[-1,1] ,'-k','linewidth',1);
    hh = scatter(C_p,H_p,50,copper(length(C_p))); hh.Marker = '.'; 
    ylabel('Hand P')
    xlabel('Cursor P')
    ylim([-.05 .05])
    xlim([-.05 .05])
    axis equal
    
    
    subplot(sb1,sb2,4)
    hold all
    plot(max(abs(C_v))*[-1,1] , [0,0],':k');
    plot([0,0], max(abs(H_v))*[-1,1] ,':k');
    plot(max(abs(C_v))*[1,-1] , max(abs(H_v))*[-1,1] ,'-k','linewidth',1);
    hh = scatter(C_v,H_v,50,copper(length(C_p))); hh.Marker = '.';
    ylabel('Hand V')
    xlabel('Cursor V')
    axis square
    
    subplot(sb1,sb2,5)
    hold all
    plot(max(abs(H_p))*[-1,1] , [0,0],':k');
    plot([0,0], max(abs(H_v))*[-1,1] ,':k');
    hh = scatter(H_p,H_v,50,copper(length(H_p))); hh.Marker = '.';
    ylabel('Hand P')
    xlabel('Hand V')
    axis square

    
    
    
end












%% 

function Sim = model_1(Flag, simdata,Init,Final,Obs,qq,C,Ke)

% Flag = 1: solve for optimal controller
% Flag = 0: generate trajectories based on solved controller


% System 
tau = simdata.tau;
M   = simdata.mass;
L   = simdata.Lambda;

simdata.B = [0 0 0 0 1/(tau)]';
simdata.A = [0 1 0 0 0
             L^2 0 L^2 L 0
             0 0 0 1 0
             0 0 0 0 1/M
             0 0 0 0 -1/tau];         

simdata.H = Obs;

% Boundry conditions
simdata.xinit  = Init;  
simdata.xfinal = Final;


% Shape the time-dependent function for Q matrix

runningalpha = zeros(size(simdata.A,1),simdata.nStep);
runningalpha(:,1) = qq;
for i = 2:simdata.nStep
    runningalpha(:,i) = qq;
end
simdata.ralpha = runningalpha;


% Run the minmax control simulation
if Flag
    simul = ms_OFC_solver(simdata);
else
    simul = ms_OFC_solver(simdata,C,Ke);
end


% Output
Sim        = struct;
Sim.QList  = qq;
Sim.Lambda = L;
Sim.Obs    = Obs;

% Actual states of the system
Sim.C_p = simul.z(:,1);
Sim.C_v = simul.z(:,2);
Sim.H_p = simul.z(:,3);
Sim.H_v = simul.z(:,4);
Sim.H_f = simul.z(:,5);


% Control command
Sim.u = simul.v;

% Control gain
Sim.C = simul.C;

% Kalman Gain
Sim.K = simul.K;


end










