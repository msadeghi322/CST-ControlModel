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

          L    = 4.0;
simdata.tau    = .07;    % muscle parameter
simdata.Lambda = L;
simdata.mass   = 1;
simdata.delta  = .01;        % Discretization step: 10ms
simdata.delay  = .05;        % feedback loop delay, xx time steps
simdata.time   = 8;         % Reach time
simdata.nStep  = simdata.time/simdata.delta+1;         % Number of time steps corresponding to reach time (600ms), plus terminal step
simdata.noise  = [.1 2.05];   % Motor, and Signal dependent noise, standard values.
simdata.effort = 10;
Time           = (0:simdata.delta:simdata.time);
simdata.nsimu  = 100; % Number of simulation runs.

EffectiveTime  = 6;

% Initial conditions
hp0  = .0;
hv0  = 0;
cp0  = 0.0;
cv0  = L*(cp0+hp0); 
F0   = 0; 
Init = [cp0,cv0,hp0,hv0,F0]';
Finl = [0 0 0 0 0]';

% feedback type (Matrix H): cursor(c)/hand(h); pos(p)/vel(v)/acc(a): Assuming only one
% observation
fb_cp = 1;
fb_cv = 1;
fb_hp = 1;
fb_hv = 1;
fb_f  = 1;
Obs = licols( diag([fb_cp fb_cv fb_hp fb_hv fb_f]) )';

SensoryNoise = 1*diag([1e-6 , 1e-6 , 1e-6 , 1e-6  , 1e-6]);
Omega = licols(Obs*SensoryNoise);
if isempty(Omega)
    Omega=0;
end
simdata.Omega = Omega;


% Q matrix: which state should be penalized more
Q_cp = 0e5;
Q_cv = 1e10;
Q_hp = 0;
Q_hv = 0;
Q_f  = 0;
qq = [Q_cp , Q_cv , Q_hp , Q_hv , Q_f]';



%% Control at muscle activity level

%%%%% solve 1 trial and generate many trials
Flag1 = 1; % solve
Sim1 = model_1(Flag1,simdata,Init,Finl,Obs,qq,[],[]);
C = Sim1.C;
Ke = Sim1.K;
for i=1:simdata.nsimu
    Flag2 = 0;
    Sim(i) = model_1(Flag2,simdata,Init,Finl,Obs,qq,C,Ke);
    fprintf('Trial %d \n',i)
end


%%%%% Solve all the trials
% for i=1:simdata.nsimu
%     Flag = 1;
%     Sim(i) = model_1(Flag,simdata,Init,Finl,Obs,qq,[],[]);
%     fprintf('Trial %d \n',i)
% end


%% Success%
Success = zeros(simdata.nsimu,1);
for i = 1:simdata.nsimu

    C_p  = Sim(i).C_p;
    ii = abs(C_p)>.05;
    if sum(ii)==0
        Success(i,1)=1;
    end
    

end

fprintf('Success rate: %.1f \n',sum(Success)/length(Success)*100);




%%
figure(100)
sb1=5;
sb2=5;
clf
for i = 1:min(sb1*sb2,simdata.nsimu)
    
    C_p  = Sim(i).C_p;
    C_v  = Sim(i).C_v;
    H_p  = Sim(i).H_p;
    H_v  = Sim(i).H_v;
    H_f  = Sim(i).H_f;

    
    CG   = squeeze(Sim(i).C);
    U    = Sim(i).u;
    KG   = Sim(i).K;
    CG_all(i,:) = nanmean( CG(1:6,:) , 2); 
    
    ind = (1:EffectiveTime/simdata.delta)';
    % Puts hold on to add simulation to the figure
    subplot(sb1,sb2,i)
    hold all
    plot([0,EffectiveTime],[0,0],':k')
    H1 = plot(Time(ind),C_p(ind),'b','linewidth',2);
    H2 = plot(Time(ind),H_p(ind),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    ylim([-.05 .05])
    set(gca,'fontsize',11)
    if i==1
       title(sprintf('L = %.2f',L)) 
    end
end
legend([H1,H2],'C','H')


%%

figure(200)
sb1=4;
sb2=4;
clf
for i = 1:min(sb1*sb2,simdata.nsimu)/2
    
    C_p  = Sim(i).C_p;
    C_v  = Sim(i).C_v;
    H_p  = Sim(i).H_p;
    H_v  = Sim(i).H_v;
    H_f  = Sim(i).H_f;
    CG   = squeeze(Sim(i).C);
    U    = Sim(i).u;
    KG   = Sim(i).K;
    CG_all(i,:) = nanmean( CG(1:6,:) , 2); 
    
    ind = (1:EffectiveTime/simdata.delta)';
    % Puts hold on to add simulation to the figure
    
    subplot(sb1,sb2,2*i-1)
    hold all
    plot([0,EffectiveTime],[0,0],':k')
    H1 = plot(Time(ind),C_p(ind),'b','linewidth',2);
    H2 = plot(Time(ind),H_p(ind),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    ylim([-.05 .05])
    set(gca,'fontsize',11)
    if i==1
       title(sprintf('L = %.2f',L)) 
    end
    
    subplot(sb1,sb2,2*i)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], .06*[-1,1] ,':k');
    plot(C_p(ind),C_v(ind),'color',[.6,.6,.6])
    hh = scatter(C_p(ind),C_v(ind),50,copper(length(C_v(ind)))); hh.Marker = '.';
    xlabel('C_P (m)')
    ylabel('C_V (m/s)')
    xlim([-.05 .05])
    ylim([-.06 .06])
    
    
end
legend([H1,H2],'C','H')




%% 
return
for i = 1:simdata.nsimu
    
    C_p  = Sim(i).C_p;
    C_v  = Sim(i).C_v;
    H_p  = Sim(i).H_p;
    H_v  = Sim(i).H_v;
    H_f  = Sim(i).H_f;
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
    xlim([-.05 .05])
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
% simdata.A = [0 1 0 0 0 0
%              0 0 1 0 0 0
%              L^3 0 0 L^3 L^2 L
%              0 0 0 0 1 0
%              0 0 0 0 0 1
%              0 0 0 0 0 -1/tau];   

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

Sim.QList  = qq;
Sim.Lambda = L;
Sim.Obs    = Obs;

% Actual states of the system
Sim.C_p = simul.z(:,1);
Sim.C_v = simul.z(:,2);
Sim.H_p = simul.z(:,3);
Sim.H_v = simul.z(:,4);
Sim.H_f = simul.z(:,5);

% Estimated states of the system using LQG
Sim.C_p_e = simul.zest(:,1);
Sim.C_v_e = simul.zest(:,2);
Sim.H_p_e = simul.zest(:,3);
Sim.H_v_e = simul.zest(:,4);
Sim.H_f_e = simul.zest(:,5);

% Control command
Sim.u = simul.v;

% Control gain
Sim.C = simul.C;

% Kalman Gain
Sim.K = simul.K;


end



%%










