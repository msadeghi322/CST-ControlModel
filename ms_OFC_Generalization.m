
% Here the OFC is solved for a given lambda and then used for a range of
% other lambdas that it has never seen



clear 
clc
close all
set(0,'DefaultFigureWindowStyle','docked') 
set(0,'defaultAxesFontSize',13)
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
set(groot,'defaultAxesBox','off')
set(0, 'DefaultFigureRenderer', 'painters');


%%% LQG 
%   scritp_minmax_PointMass
%
%   Writtent by F. Crevecoeur - Spet. 6, 2019
%   Used in: Robust control in human reaching movements: a model free
%   strategy to compensate for unpredictable disturbances. 
%   Crevecoeur F., Scott S. H., Cluff T. 
%   DOI: https://doi.org/10.1523/JNEUROSCI.0770-19.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


simdata.tau    = .06;    % muscle parameter
simdata.mass   = 1;
simdata.delta  = .01;        % Discretization step: 10ms
simdata.delay  = .05;        % feedback loop delay, xx time steps
simdata.time   = 8;         % Reach time
simdata.nStep  = simdata.time/simdata.delta+1;         % Number of time steps corresponding to reach time (600ms), plus terminal step
simdata.noise  = [.01 .4];   % Motor, and Signal dependent noise, standard values.
%simdata.noise = [.4 1 1.5];       % Motor, Sensory and Signal dependent noise, standard values.
simdata.effort = 200;
Time           = (0:simdata.delta:simdata.time);
simdata.nsimu  = 150; % Number of simulation runs.

EffectiveTime  = 6;


% feedback type (Matrix H): cursor(c)/hand(h); pos(p)/vel(v)/acc(a): Assuming only one
% observation
fb_cp = 1;
fb_cv = 1;
fb_ca = 0;
fb_hp = 0;
fb_hv = 0;
fb_ha = 0;
Obs = licols( diag([fb_cp fb_cv fb_ca fb_hp fb_hv fb_ha]) )';
SensoryNoise = 10*diag([1e-6 , 1e-5 , 1e-3 , 1e-3 , 1e-3 , 1e-2]);
Omega = licols(Obs*SensoryNoise);
if isempty(Omega)
    Omega=0;
end
simdata.Omega = Omega;


% Q matrix: which state should be penalized more
Q_cp = 1e5;
Q_cv = 1e-1;
Q_ca = 1e-1;
Q_hp = 1e-1;
Q_hv = 1e-1;
Q_ha = 1e-1;
qq = [Q_cp , Q_cv , Q_ca , Q_hp , Q_hv , Q_ha]';



%% Cross validation of lambda
Lambda_List_Train = [3,4,5,6,7,8];
Lambda_List_Test  = [2:.5:9];
Cpolicy = qq>10; Cpolicy=double(Cpolicy);

FilePath = 'CrossValidSim';
flp = sprintf('%s/LRes%d',FilePath,length(Lambda_List_Test));
fln = sprintf('%s/Trials%d_P%d_V%d_D%d.mat',flp,simdata.nsimu,Cpolicy(1),Cpolicy(2),simdata.delay*1000);
if ~exist(flp,'file')
    mkdir(flp)
end

if exist(fln,'file')
    fprintf('%s already exists \n',fln);
    D = load(fln);
    Sim = D.Sim;
else
    
    for i=1:length(Lambda_List_Train)
        L1 = Lambda_List_Train(i); % Training lambda
        simdata.Lambda = L1;
        
        % Initial conditions (assumed zero)
        hp0  = .0;
        hv0  = 0;
        ha0  = 0;
        cp0  = 0.0;
        cv0  = L1*(cp0+hp0);
        ca0  = L1*(cv0+hv0);
        Init = [cp0,cv0,ca0,hp0,hv0,ha0]';
        Finl = [0 0 0 0 0 0]';
        
        %%%%% solve 1 trial and generate many trials
        Flag1 = 1; % solve
        Sim1 = model_1(Flag1,simdata,Init,Finl,Obs,qq,[],[]);
        C = Sim1.C;
        Ke = Sim1.K;
        
        for j=1:length(Lambda_List_Test)
            L2 = Lambda_List_Test(j); % Testing lambda
            simdata.Lambda = L2;
            for n=1:simdata.nsimu
                Flag2 = 0;
                Sim(i,j,n) = model_1(Flag2,simdata,Init,Finl,Obs,qq,C,Ke);
                fprintf('Trial %d \n',n)
            end
            
        end
    end
    save(fln,'Sim','-v7.3');
    
    
end


%% Success%
Success = zeros(simdata.nsimu,length(Lambda_List_Train),length(Lambda_List_Test));
for i=1:length(Lambda_List_Train)
    for j=1:length(Lambda_List_Test)
        for n=1:simdata.nsimu
            C_p  = Sim(i,j,n).C_p;
            ii = abs(C_p)>.05;
            if sum(ii)==0
                Success(n,i,j)=1;
            end
            
        end
    end
end

SuccessRate = squeeze(mean(Success,1));

clf
for i=1:length(Lambda_List_Train)
subplot(2,3,i)
hold all
plot(Lambda_List_Test,SuccessRate(i,:),'.-','linewidth',1.5,'markersize',15)
kk = Lambda_List_Test==Lambda_List_Train(i);
plot(Lambda_List_Train(i),SuccessRate(i,kk),'.','linewidth',1.5,'markersize',20)
ylim([-.2 1.2])
text(2.5, 1.1, sprintf('Training \\lambda = %.1f',Lambda_List_Train(i)),'fontsize',12)
%set(gca,'xtick',Lambda_List_Test(2:2:end))
grid on
xlabel('Testing \lambda')
ylabel('Success Rate ')
end

%%
i = 4; % training condition
figure(100)
sb1=5;
sb2=ceil(length(Lambda_List_Test)/3);
clf
pp=0;
for j=3:3:length(Lambda_List_Test)
    pp = pp+1;
    for n = 1:sb1
        
        C_p  = Sim(i,j,n).C_p;
        C_v  = Sim(i,j,n).C_v;
        C_a  = Sim(i,j,n).C_a;
        H_p  = Sim(i,j,n).H_p;
        H_v  = Sim(i,j,n).H_v;
        H_a  = Sim(i,j,n).H_a;
        
        ind = (1:EffectiveTime/simdata.delta)';
        % Puts hold on to add simulation to the figure
        subplot(sb1,sb2,pp+(n-1)*sb2)
        hold all
        plot([0,EffectiveTime],[0,0],':k')
        H1 = plot(Time(ind),C_p(ind),'b','linewidth',2);
        H2 = plot(Time(ind),H_p(ind),'r','linewidth',2);
        %xlabel('Time (s)')
        %ylabel('Position')
        ylim([-.05 .05])
        set(gca,'fontsize',11)
        if n==1
            title(sprintf('L = %.2f',Lambda_List_Test(j)))
        end
    end
    
end
legend([H1,H2],'C','H')






return
%% 
for i = 1:simdata.nsimu
    
    C_p  = Sim(i).C_p;
    C_v  = Sim(i).C_v;
    C_a  = Sim(i).C_a;
    H_p  = Sim(i).H_p;
    H_v  = Sim(i).H_v;
    H_a  = Sim(i).H_a;
    CG   = squeeze(Sim(i).C);
    U    = Sim(i).u;
    KG   = Sim(i).K;
    CG_all(i,:) = nanmean( CG(1:6,:) , 2); 
    
    ii = C_p>.05;
    if sum(ii)>0
        Success(i)=0;
    else
        Success(i)=1;
    end
    
    if i>10
       continue
    end
    figure(i)
    sb1 = 3;
    sb2 = 4;
    clf
    ind = (1:EffectiveTime/simdata.delta)';
    % Puts hold on to add simulation to the figure
    subplot(sb1,sb2,1)
    hold all
    plot([0,EffectiveTime],[0,0],':k')
    H1 = plot(Time(ind),C_p(ind),'b','linewidth',2);
    H2 = plot(Time(ind),H_p(ind),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    ylim([-.05 .05])
    
    subplot(sb1,sb2,2)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], .06*[-1,1] ,':k');
    plot(C_p(ind),C_v(ind),'color',[.6,.6,.6])
    hh = scatter(C_p(ind),C_v(ind),50,copper(length(C_v(ind)))); hh.Marker = '.';
    xlabel('Cursor P')
    ylabel('Cursor V')
    xlim([-.03 .03])
    ylim([-.06 .06])
    
    subplot(sb1,sb2,3)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], .05*[-1,1] ,'-k','linewidth',1);
    plot(.05*[1,-1] , .05*[-1,1] ,'-k','linewidth',1);
    hh = scatter(C_p(ind),H_p(ind),50,copper(length(C_p(ind)))); hh.Marker = '.'; 
    ylabel('Hand P')
    xlabel('Cursor P')
    ylim([-.05 .05])
    xlim([-.05 .05])
    axis equal
    
    
    subplot(sb1,sb2,4)
    hold all
    plot(max(abs(H_p(ind)))*[-1,1] , [0,0],':k');
    plot([0,0], max(abs(H_v(ind)))*[-1,1] ,':k');
    hh = scatter(H_p(ind),H_v(ind),50,copper(length(H_p(ind)))); hh.Marker = '.';
    xlabel('Hand P')
    ylabel('Hand V')
    axis square

    
    subplot(sb1,sb2,5)
    hold all
    [yy,xx] = xcorr(-H_p(ind),C_p(ind)); xx = xx*10; % to ms
    R(i) = corr(C_p(ind),H_p(ind));
    [mm,imx]=max(abs(yy));
    plot(xx,yy,'linewidth',2);
    plot([xx(1),xx(end)],[0,0],':k');
    plot([0 0],mm*[-2,2],':k');
    plot(xx(imx)*[1 1],mm*[0 1],'-k')
    text(2000,mm,sprintf('lag = %.2f \n Corr = %.2f',xx(imx),R(i)));
    xlabel('Lag ms')
    ylabel('XCorr')
    ylim(1.2*mm*[-.5,1])
    Lag(i) = xx(imx);
    
    
    subplot(sb1,sb2,6)
    hold all
    plot([0,Time(ind(end))],[0,0],':k')
    H1 = plot(Time(ind),C_v(ind),'b','linewidth',2);
    H2 = plot(Time(ind),L*(H_p(ind)+C_p(ind)),'--r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Vel.')
    legend([H1,H2],'dx','L*(x+u)','location','best')
    
    
    subplot(sb1,sb2,7)
    hold all
    H1 = plot(Time(ind),CG(1:3,ind)','linewidth',2);
    plot([0,EffectiveTime],[0,0],':k')
    xlabel('Time (s)')
    ylabel('Gain')
    legend(H1,'Pos','Vel','Acc','location','best')
    
    
    subplot(sb1,sb2,8)
    hold all
    plot([0,EffectiveTime],[0,0],':k')
    H1 = plot(Time(ind),C_v(ind),'b','linewidth',2);
    H2 = plot(Time(ind),H_v(ind),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Vel')
    legend([H1,H2],'Cursor','Hand','location','best')
    ylim([-.1 .1])
    
    
    
    subplot(sb1,sb2,9)
    hold all
    F = simdata.mass*H_a;
    H1 = plot(Time(ind),U(ind),'linewidth',1.5);
    plot(Time(ind) , F(ind),'r','linewidth',2)
    legend('u','F')
    plot([0,EffectiveTime],[0,0],'k')
    xlabel('Time (s)')
    ylabel('Motor command')
    
    
    
end





%%

function Sim = model_1(Flag, simdata,Init,Final,Obs,qq,C,Ke)

% Flag = 1: solve for optimal controller
% Flag = 0: generate trajectories based on solved controller


Sim = struct;

% System 
tau = simdata.tau;
M   = simdata.mass;
L   = simdata.Lambda;


% observability rank = Full (6)
% Controllability rank = 4
simdata.A = [0 1 0 0 0 0
             0 0 1 0 0 0
             L^3 0 0 L^3 L^2 L
             0 0 0 0 1 0
             0 0 0 0 0 1
             0 0 0 0 0 -1/tau];   

% % observability rank = 5          
% % Controllability rank = 4
% simdata.A = [L 0 0 L 0 0
%              0 L 0 0 L 0
%              0 0 L 0 0 L
%              0 0 0 0 1 0
%              0 0 0 0 0 1
%              0 0 0 0 0 -1/tau];
%          
% % observability rank = 4          
% % Controllability rank = 4
% simdata.A = [0 1 0 0 0 0
%              0 0 1 0 0 0
%              0 0 L 0 0 L
%              0 0 0 0 1 0
%              0 0 0 0 0 1
%              0 0 0 0 0 -1/tau];   
         
         
         
simdata.B = [0 0 0 0 0 1/(tau*M)]';
simdata.H = Obs;

% Boundry conditions
simdata.xinit  = Init;  
simdata.xfinal = Final;


% Shape the time-dependent function for Q matrix

runningalpha = zeros(size(simdata.A,1),simdata.nStep);
runningalpha(:,1) = qq;
for i = 2:simdata.nStep
    runningalpha(:,i) = qq;
    %runningalpha(:,i) = (1-exp((-i/10)))*qq;
    %runningalpha(:,i) = runningalpha(:,i-1)+exp((-i/10))*qq;
    %runningalpha(:,i) = exp((-(i-300)^2)/1200)*qq;
end
simdata.ralpha = runningalpha;


% Run the minmax control simulation
if Flag
    simul = ms_OFC_solver(simdata);
else
    simul = ms_OFC_solver(simdata,C,Ke);
end

Sim.QList  = qq;
Sim.Lambda = L;
Sim.Obs    = Obs;

% Actual states of the system
Sim.C_p = simul.z(:,1);
Sim.C_v = simul.z(:,2);
Sim.C_a = simul.z(:,3);
Sim.H_p = simul.z(:,4);
Sim.H_v = simul.z(:,5);
Sim.H_a = simul.z(:,6);

% Estimated states of the system using LQG
Sim.C_p_e = simul.zest(:,1);
Sim.C_v_e = simul.zest(:,2);
Sim.C_a_e = simul.zest(:,3);
Sim.H_p_e = simul.zest(:,4);
Sim.H_v_e = simul.zest(:,5);
Sim.H_a_e = simul.zest(:,6);

% Control command
Sim.u = simul.v;

% Control gain
Sim.C = simul.C;

% Kalman Gain
Sim.K = simul.K;


end



%%










