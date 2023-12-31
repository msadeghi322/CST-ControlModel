
clear 
clc
close all
set(0,'DefaultFigureWindowStyle','docked') 
set(0,'defaultAxesFontSize',16)
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
set(groot,'defaultAxesBox','off')
set(0, 'DefaultFigureRenderer', 'painters');

%mutual information path
addpath('D:\OneDrive - Northeastern University\Action Lab\01 Projects\External-Matlab-Toolbox\MIToolbox\MI2022\mi');


%%

% Lambda case: plotting agaist:
%_____ Lambda:                     Case=1
%_____ Lambda/Lambda_c:            Case=2
%_____ (Lambda-Lambda_c)/Sigma_c:  Case=3;
%_____ Success rate:               Case=4

LambdaCase = 1;
LCLabels = {'L','LC','LN','Success'};

StrategyList = [3,6,9]; % 1:QpFp, 2:QpFv, 3:QpFpv , ... , 6:QvFpv, 9:QpvFpv


SaveFigures = 0;

%% Analysis

FileName = 'sim';


Case = 1;

switch Case
    
    case 1
        %%%%  Full model : R = 1;
        %FilePath = 'Reduced data/OFC_R1';
        %FilePath = 'Reduced data/OFC_R1_SfN_A3';
        %FilePath = 'Reduced data/OFC_R1_SfN_A1';
        %FilePath = 'Reduced data/OFC_R1_SfN_A1_V02';
        %FilePath = 'Reduced data/OFC_R1_SfN_A1_V04';
        %FilePath  = 'Reduced data/OFC_R1_SfN_A1_V05';
        %FilePath = 'Reduced data/OFC_R1_SfN_A1_V06';
        %FilePath = 'Reduced data/OFC_R1_SfN_A1_V07';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V01';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V02';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V03';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V04';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V05';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V06';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V09';
        %FilePath  = 'Reduced data/OFC_R1_COSYNE_A1_V10'; 
        FilePath  = 'Reduced data/CST_Paper_V01'; 
        
        
        
        Obs_List = [1 0 0 0 0 0
            0 1 0 0 0 0
            1 1 0 0 0 0];
        Q_List   = [1e3 1 1 1 1 1 ; 1 1e3 1 1 1 1 ; 1e3 1e3 1 1 1 1];

        
        
        
        
        
        fpth = sprintf('%s/R1',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
    case 2
        %%%%  Full model : R = 1e5;
        FilePath = 'Reduced data/OFC_R1e5';
        fpth = sprintf('%s/R1e5',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
    case 3
        %%%%  Full model : R = 1e-5;
        FilePath = 'Reduced data/OFC_R1e_5';
        fpth = sprintf('%s/R1e_5',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
        
    case 4
        %%%%  No signal dependent noise : R = 1;
        FilePath = 'Reduced data/OFC_R1_lqg';
        fpth = sprintf('%s/R1_lqg',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .008;
        
    case 5
        %%%%  Full model : R = 1;
        FilePath = 'Reduced data/OFC_R1e_5_lqg';
        fpth = sprintf('%s/R1e_5_lqg',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .008;
        
    case 6
        %%%%  Full model : R = 1;
        FilePath = 'Reduced data/OFC_R1e5_lqg';
        fpth = sprintf('%s/R1e_5_lqg',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
    case 7
        %%%%  Full model : R = 1;
        FilePath = 'Reduced data/OFC_R1e_5_lqg_D0';
        fpth = sprintf('%s/R1e_5_lqg_D0',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .008;
        
    case 8
        %%%%  Full model : R = 1;
        FilePath = 'Reduced data/OFC_R1_Delay0';
        fpth = sprintf('%s/R1_Delay0',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
    case 9
        %%%%  Full model : R = 1;
        FilePath = 'Reduced data/OFC_R1_D100';
        fpth = sprintf('%s/R1_D100',FilePath);
        Lambda_List = [.5:.4:10]';
        duration = 8;
        Effective_Duration = 6; % this is what is analysed
        FailBound = .05;
        
        
        
end



Conditions = {'C_p,fb_p','C_p,fb_v','C_p','C_v,fb_p','C_v,fb_v','C_v','C_p_v,fb_p','C_p_v,fb_v','C_p_v'};
delta    = .01;
Time     = (0:delta:duration)';

if exist(sprintf('%s.mat',fpth),'file')
    S = load(sprintf('%s.mat',fpth));
    D = S.D;
    clear S;
    fprintf('%s Loaded. \n',fpth)
    
else
    
    CondN = 9;
    TRN   = 500;
    Lag_P       = NaN(TRN,length(Lambda_List), CondN);
    Lag_V       = NaN(TRN,length(Lambda_List), CondN);
    Lag_A       = NaN(TRN,length(Lambda_List), CondN);
    Corr_P      = NaN(TRN,length(Lambda_List), CondN);
    Corr_V      = NaN(TRN,length(Lambda_List), CondN);
    Corr_A      = NaN(TRN,length(Lambda_List), CondN);
    MI_P        = NaN(TRN,length(Lambda_List), CondN);
    MI_V        = NaN(TRN,length(Lambda_List), CondN);
    MI_A        = NaN(TRN,length(Lambda_List), CondN);
    HandRMS_P   = NaN(TRN,length(Lambda_List), CondN);
    HandRMS_V   = NaN(TRN,length(Lambda_List), CondN);
    CursRMS_P   = NaN(TRN,length(Lambda_List), CondN);
    CursRMS_V   = NaN(TRN,length(Lambda_List), CondN);
    HandCursRMS_P = HandRMS_P;  % RMS of the diff between hand and cursor: RMS(diff), NOT diff(RMS) 
    HandCursDiff_P = HandCursRMS_P ; % diff between H and C
    Success       = NaN(TRN,length(Lambda_List), CondN);
    ControlGain   = NaN(TRN,6,length(Lambda_List), CondN);
    ControlGain_I = NaN(TRN,6,length(Lambda_List), CondN);
    ControlGain_F = NaN(TRN,6,length(Lambda_List), CondN);
    KalmanGain    = NaN(TRN,6,length(Lambda_List), CondN);
    SuccessT      = NaN(TRN,length(Lambda_List), CondN);
    ControlComm   = NaN(TRN,length(Lambda_List), CondN);
    HC_XpassCount = ControlComm;
    HandCursRMS_v = HandCursRMS_P;
    CursorMaxVel  = HandCursRMS_P;
    CursorMeanPos = HandCursRMS_P;
    CursorMeanVel = HandCursRMS_P;
    
    Cond     = 0;
    for i=1:size(Q_List,1)
        for j=1:size(Obs_List,1)
            Cond = Cond+1;
            
            for k=1:length(Lambda_List)
                
                flp = sprintf('%s/Q%d_H%d',FilePath,i,j);
                fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
                fprintf('%s \n',fln)
                D = load(fln);
                Sim = D.Sim;
                clear D;
                FN = fieldnames(Sim);
                
                % loop over trials
                for n=1:length(Sim)
                    for fn=1:length(FN)
                        eval(sprintf('%s = Sim(n).%s;',FN{fn},FN{fn}));
                    end
                    
                    % calculate lag
                    ii = (1:length(C_p))'<=Effective_Duration/delta;
                    [yy,xx] = xcorr(-(H_p(ii)-mean(H_p(ii))),C_p(ii)-mean(C_p(ii))); xx = xx*10; % to ms
                    [~,ind]=max(abs(yy));
                    Lag_P(n,k,Cond) = xx(ind);
                    
                    
                    [yy,xx] = xcorr(-(H_v(ii)-mean(H_v(ii))),C_v(ii)-mean(C_v(ii))); xx = xx*10; % to ms
                    [~,ind]=max(abs(yy));
                    Lag_V(n,k,Cond) = xx(ind);
                    
                    [yy,xx] = xcorr(-(H_a(ii)-mean(H_a(ii))),C_a(ii)-mean(C_a(ii))); xx = xx*10; % to ms
                    [~,ind]=max(abs(yy));
                    Lag_A(n,k,Cond) = xx(ind);
                    
                    
                    % calculate correlation
                    %Corr_P(n,k,Cond) = corr(C_p(ii),H_p(ii));
                    %Corr_V(n,k,Cond) = corr(C_v(ii),H_v(ii));
                    %Corr_A(n,k,Cond) = corr(C_a(ii),H_a(ii));
                    rp = xcorr(C_p(ii)-mean(C_p(ii)),H_p(ii)-mean(H_p(ii)),'coeff');
                    Corr_P(n,k,Cond) = min(rp);
                    rv = xcorr(C_v(ii)-mean(C_v(ii)),H_v(ii)-mean(H_v(ii)),'coeff');
                    Corr_V(n,k,Cond) = min(rv);
                    ra = xcorr(C_a(ii)-mean(C_a(ii)),H_a(ii)-mean(H_a(ii)),'coeff');
                    Corr_A(n,k,Cond) = min(ra);
                    
                    % calculate Mutual information
                    kk = abs(C_p(ii))>2*FailBound;
                    if sum(kk)==0
                        ccp = (C_p(ii)-mean(C_p(ii)))*1000; % has to be in mm to be consistent with data- also for discritization
                        hhp = (H_p(ii)-mean(H_p(ii)))*1000;
                        MI_P(n,k,Cond) = mutualinfo(round(ccp),round(hhp));
                        
                        ccv = (C_v(ii)-mean(C_v(ii)))*1000;
                        hhv = (H_v(ii)-mean(H_v(ii)))*1000;
                        MI_V(n,k,Cond) = mutualinfo(round(ccv),round(hhv));
                        
                        %cca = (C_a(ii)-mean(C_a(ii)))*1000;
                        %hha = (H_a(ii)-mean(H_a(ii)))*1000;
                        %MI_A(n,k,Cond) =mutualinfo(round(cca),round(hha));
                    else
                        fprintf('Trial %d  ---  L%d  --- Cond%d ---- MI problem \n',n,k,Cond);
                    end
                    
                    % Calculate deviation (only first 5 seconds to make sure
                    % the termination effects are not involved.
                    jj = abs(C_p)>FailBound;
                    if sum(jj)==0
                        ind = find(ii==1);
                        %ind = find(ii.*jj==1);
                        HandRMS_P(n,k,Cond) = rms(H_p(ind));
                        CursRMS_P(n,k,Cond) = rms(C_p(ind));
                        HandRMS_V(n,k,Cond) = rms(H_v(ind));
                        CursRMS_V(n,k,Cond) = rms(C_v(ind));
                        HandCursRMS_P(n,k,Cond) = rms( H_p(ind) - C_p(ind) );
                        HandCursDiff_P(n,k,Cond) = abs( sum(H_p(ind)) - sum(C_p(ind)) );
                        HandCursRMS_v(n,k,Cond) = rms( H_v(ind) - C_v(ind) );
                        CursorMeanPos(n,k,Cond) = mean(C_p(ind));
                        CursorMeanVel(n,k,Cond) = mean(C_v(ind));
                        
                    end
                    
                    
                    % Control gain average 
                    ControlGain(n,:,k,Cond)   = squeeze( nanmean(Sim(1).C(1,1:6,ii),3) );
                    ControlGain_I(n,:,k,Cond) = squeeze( nanmean(Sim(1).C(1,1:6,1:100),3) );
                    ControlGain_F(n,:,k,Cond) = squeeze( nanmean(Sim(1).C(1,1:6,500:600),3) );
                    
                    
                    % Kalman Gain
                    KalmanGain(n,:,k,Cond) = nanmean(Sim(1).K(1:6,1,ii),3);
                    
                    
                    
                    % Percent success (Fail after 6 s doesn't count)
                    kk = abs(C_p(ii))>FailBound;
                    if sum(kk)==0
                        Success(n,k,Cond) = 1;
                    else
                        Success(n,k,Cond) = 0;
                    end
                    
                    
                    % Time of fail: % of duration the cursor was
                    % successfully kept within the boundaries
                    jj = abs(C_p(ii))>FailBound;
                    if sum(jj)==0
                        SuccessT(n,k,Cond) = 100;
                    else
                        ind_f = find(jj==1,1,'first');
                        SuccessT(n,k,Cond) = 100*ind_f/length(jj); 
                    end
                    
                    
                    % Control command
                    jj = abs(C_p(ii))>FailBound;
                    if sum(jj)==0
                        ControlComm(n,k,Cond) = trapz(delta,u(ii).^2);
                    end
                    
                    
                    %\____________________________________________________/
                    %Measures that could help identify control policy: 
                    %Vel vs Pos control
                    %/^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\
                    
                    % Times cursor and hand crossed passed: this is a
                    % measure of determining pos vs vel control: if hand
                    % and cursor cross pass a lot -> pos control, if not,
                    % vel control.
                    if k==0
                        S = H_p(ii);
                        S2 = C_p(ii);
                        Y = fft(S);
                        Y2 = fft(S2);
                        L = length(S);
                        Fs = 1/.01;
                        f = Fs*(0:(L/2))/L;
                        P2 = abs(Y/L);
                        P22 = abs(Y2/L);
                        P1 = P2(1:L/2+1);
                        P12 = P22(1:L/2+1);
                        P1(2:end-1) = 2*P1(2:end-1);
                        P12(2:end-1) = 2*P12(2:end-1);
                        
                        subplot(2,2,1)
                        plot(C_p(ii)); hold on
                        plot(H_p(ii)); hold off
                        
                        subplot(2,2,2)
                        h1 = plot(f,P12,'.-'); hold on
                        h2 = plot(f,P1,'.-'); hold on
                        title('Single-Sided Amplitude Spectrum of S(t)')
                        xlabel('f (Hz)')
                        ylabel('|P1(f)|')
                        xlim([0 10])
                        legend([h1,h2],'Cursor','Hand')
                        hold off
                    end
                    
                    kk = abs(C_p(ii))>FailBound;
                    if sum(kk)==0
                        Crss = -abs(H_p(ii) - C_p(ii));
                        [pks,lpks] = findpeaks(Crss, 'minPeakProminence', .001);
                        jj = find(pks<-.01); pks(jj)=[]; lpks(jj)=[];
                        HC_XpassCount(n,k,Cond) = length(pks);
                        % \__ notes: hyper parameters are the bounds for peak and for exclusion
                    end
                    
                    
                    
                    % Peak velocity
                    kk = abs(C_p(ii))>FailBound;
                    if sum(kk)==0
                        mxv = max(abs(C_v(ii)));
                        CursorMaxVel(n,k,Cond) = mxv;
                        % \__ notes: hyper parameters are the bounds for peak and for exclusion
                    end
                    
                    
                    
                    
                    
                end
            end
        end
    end
    
    D.Lag_P       = Lag_P;
    D.Lag_V       = Lag_V;
    D.Lag_A       = Lag_A;
    D.Corr_P      = Corr_P;
    D.Corr_V      = Corr_V;
    D.Corr_A      = Corr_A;
    D.MI_P        = MI_P;
    D.MI_V        = MI_V;
    D.MI_A        = MI_A;
    D.HandRMS_P   = HandRMS_P;
    D.HandRMS_V   = HandRMS_V;
    D.CursRMS_P   = CursRMS_P;
    D.CursRMS_V   = CursRMS_V;
    D.HandCursRMS_P = HandCursRMS_P;
    D.HandCursDiff_P  = HandCursDiff_P;
    D.Success     = Success;
    D.ControlGain = ControlGain;
    D.ControlGain_I = ControlGain_I;
    D.ControlGain_F = ControlGain_F;
    D.KalmanGain  = KalmanGain;
    D.SuccessT    = SuccessT;
    D.ControlComm = ControlComm;
    D.HC_XpassCount = HC_XpassCount;
    D.HandCursRMS_v = HandCursRMS_v;
    D.CursorMaxVel  = CursorMaxVel;
    D.CursorMeanPos = CursorMeanPos;
    D.CursorMeanVel = CursorMeanVel;
    
    
    save(fpth,'D','-v7.3');
    
end



%% plot settings
CL = [.8,.5,.5
      .6,0,0
      .6,.5,0
      .5,.8,.5;
      0,.4,0
      0,.7,.6
      .5,.5,.8
      0,0,.7
      .6,0,.9
      ];
%CL = get(groot,'DefaultAxesColorOrder'); CL=[CL;nthroot(flip(CL),4)];



%% 

figure(8)
clf
sb1=3;
sb2=8;
pp=0;

lmbd_idx = 1:length(Lambda_List);
Lambda_c = NaN(9,1);
Sigma_c = NaN(9,1);

% Success
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
X  = Lambda_List(lmbd_idx); 
y  = squeeze( nanmean(D.Success,1) );  y = y(lmbd_idx,:); % mean over trials
ft=fittype('50*( 1-erf((x - mu)/sqrt(2)/sigma) )', 'independent', 'x', 'dependent', 'y' );
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    B0 = [-2,4,4];
    FF = fit(X,100*y(:,i),ft, 'StartPoint',[3.5,0.3],'Robust','Bisquare');
    l = 0:.01:Lambda_List(end);
    yf = FF(l);
    %[~,ind] = min(abs(yf-.5)); % find critical lambda
    Lambda_c(i) = FF.mu;
    Sigma_c(i)  = FF.sigma;
    hh(i) = plot(l,yf,'color',CL(i,:),'linewidth',1.5);
    plot(X,100*y(:,i),'.','color',CL(i,:),'linewidth',2,'markersize',15)
    a=1;
end
plot([0,Lambda_List(end)],[50,50],':k')
ylim([0 100])
%xlim([0,Lambda_List(end)+1])
xlim([1,6])
ylabel('% Success')
xlabel('\lambda')
set(gca,'fontsize',12)
%legend(hh,Conditions,'location','best')




% Critical lambda based on % success
pp=pp+1;
subplot(sb1,sb2,sb2+1)
hold all
for kk=1:length(StrategyList)
   i = StrategyList(kk);
   bb = bar(2*kk-1,Lambda_c(i),.7);
   bb.FaceColor = CL(i,:);
end
set(gca,'xtick',2*[1:9]-1,'xticklabel',Conditions(StrategyList))
xtickangle(70)
ylim([0,Lambda_List(end)+1]);
ylabel('\lambda_C')
%legend(Conditions(StrategyList),'location','best','orientation','horizontal')
grid 
set(gca,'fontsize',12)


if 0
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
ZZ = reshape(Lambda_c,3,3);
cc=0;
for j=1:size(ZZ,2)
    for i=1:size(ZZ,1)
        cc=cc+1;
        aa = zeros(size(ZZ));
        aa(i,j)=1;
        Z = reshape(Lambda_c,3,3);
        B3 = bar3(ZZ.*aa);
        B3(1).FaceColor= CL(cc,:);
        B3(2).FaceColor= CL(cc,:);
        B3(3).FaceColor= CL(cc,:);
    end
end
set(gca,'xtick',[1:3],'xticklabel',{'C_P', 'C_V','C_P_V'})
set(gca,'ytick',[1:3],'yticklabel',{'fb_P', 'fb_V','fb_P_V'})
grid 
zlabel('\lambda_C')
zlim([min(Lambda_c)-1,8])
view(25,15)


pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y  = squeeze( nanmean(D.SuccessT,1) );  y = y(lmbd_idx,:); % mean over trials
ft=fittype('50*( 1-erf((x - mu)/sqrt(2)/sigma) )', 'independent', 'x', 'dependent', 'y' );
for i=1:size(y,2)
    B0 = [-2,4,4,5];
    FF = fit(X,y(:,i),ft,'StartPoint',[3.5,0.3],'Robust','Bisquare');
    l = 0:.01:Lambda_List(end);
    yf = FF(l);
    Lambda_c2(i) = FF.mu;
    Sigma_c2(i)  = FF.sigma;
    hh(i) = plot(l,yf,'color',CL(i,:),'linewidth',1.5);
    plot(X,y(:,i),'.','color',CL(i,:),'linewidth',2,'markersize',15)
    a=1;
end
plot([0,Lambda_List(end)],[50,50],':k')
ylim([0 100])
xlim([0,Lambda_List(end)+1])
ylabel('%Duration of success')
xlabel('\lambda')



% Critical lambda based on % success duration
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
for i=1:length(Lambda_c)
   bb = bar(2*i-1,Lambda_c2(i),.7);
   bb.FaceColor = CL(i,:);
end
set(gca,'xtick',2*[1:9]-1,'xticklabel',Conditions)
xtickangle(70)
ylim([0,Lambda_List(end)+1]);
ylabel('\lambda_C_2')
grid 


pp=pp+1;
subplot(sb1,sb2,pp)
hold all
ZZ = reshape(Lambda_c2,3,3);
cc=0;
for j=1:size(ZZ,2)
    for i=1:size(ZZ,1)
        cc=cc+1;
        aa = zeros(size(ZZ));
        aa(i,j)=1;
        Z = reshape(Lambda_c,3,3);
        B3 = bar3(ZZ.*aa);
        B3(1).FaceColor= CL(cc,:);
        B3(2).FaceColor= CL(cc,:);
        B3(3).FaceColor= CL(cc,:);
    end
end
set(gca,'xtick',[1:3],'xticklabel',{'C_P', 'C_V','C_P_V'})
set(gca,'ytick',[1:3],'yticklabel',{'fb_P', 'fb_V','fb_P_V'})
grid 
zlabel('\lambda_C_2')
zlim([min(Lambda_c)-1,8])
view(25,15)

end


if SaveFigures
    FLN = sprintf('ModelSuccessRate');
    print(FLN,'-dpng')
end









%%


sb1 = 3;
sb2 = 5;
y1 = D.CursRMS_P; y1 = y1(:,lmbd_idx,:); % mean over trials
y2 = D.CursRMS_V; y2 = y2(:,lmbd_idx,:); % mean over trials

figure(15)
clf
DiffLvl = {'Easy','Moderate','Hard'};

%____________________________________/ Regression without offset
PhaseSlope = [];
for kk=3:-1:1%length(StrategyList)
    
    i = StrategyList(kk);
    % Difficulty levels
    Thrsh = .75; % threshold for difficulty
    Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<1;
    DL(3).ii = Lambda_n1>1;
    
    for lm = 1:length(DL)  % loop over lambda
    
        subplot(sb1,sb2,kk+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:2:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:2:end);
        B = regress( Y , X );
        PhaseSlope(kk,lm) = B(1);
        Y_h =X*B; 
        h1(kk) = plot(X,Y,'.','color',CL(i,:));
        plot(X,Y_h,'-','color',CL(i,:)/2);
        ylim([0,10])
        xlim([0 5])
        text(.5,9,sprintf('a = %.2f',B(1)),'fontsize',11);
        if kk==1
        xlabel('Cursor pos. RMS (cm)')
        ylabel(sprintf('%s \n Cursor vel. RMS (cm/s)',DiffLvl{lm}))
        end
        set(gca,'fontsize',12)
        
    end    
end
legend(h1,'Pos','Vel','location','southeast')

%____________________________________/ Regression with offset

if 0
for kk=3:-1:1%length(StrategyList)
    
    i = StrategyList(kk);
    % Difficulty levels
    Thrsh = .75; % threshold for difficulty
    Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<1;
    DL(3).ii = Lambda_n1>1;
    
    for lm = 1:length(DL)  % loop over lambda
    
        subplot(sb1,sb2,kk+3+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:2:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:2:end);
        B = regress( Y , [ones(size(X)),X] );
        Y_h =[ones(size(X)),X]*B; 
        h1(kk) = plot(X,Y,'.','color',CL(i,:));
        plot(X,Y_h,'-','color',CL(i,:)/2);
        ylim([0,10])
        xlim([0 5])
        text(.5,9,sprintf('a = %.2f',B(2)),'fontsize',11);
        if kk==1
        xlabel('RMS_P (cm)')
        ylabel(sprintf('%s \n RMS_V (cm/s)',DiffLvl{lm}))
        end
        set(gca,'fontsize',11)
        
    end    
end
legend(h1,'Pos','Vel','location','southeast')
end

if SaveFigures
    FLN = sprintf('ModelPhaseRMS');
    print(FLN,'-dpng')
end

%%
figure(151)
clf
sb1 = 2;
sb2=12;
for dl=1:3
    subplot(sb1,sb2,2*dl-1)
    hold all
    
    for gr=1:2
        i = StrategyList(gr);
        DD = PhaseSlope(gr,dl);
        Y = mean(DD,1);
        S = std(DD,[],1)/sqrt(size(DD,1));
        X = gr;
        %errorbar(X,Y,S,'.k','capsize',.5)
        B = bar(X,Y,.6);
        B.EdgeColor = CL(i,:);
        B.FaceColor = CL(i,:);
        
    end
    if dl==1;  ylabel('Slope (s^{-1})');   end
    xlim([.5 2.5])
    ylim([0,2.5])
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb1 = 2;
sb2 = 2;

y1 = D.CursRMS_P; y1 = y1(:,lmbd_idx,:); % mean over trials
y2 = D.CursRMS_V; y2 = y2(:,lmbd_idx,:); % mean over trials


i1 = StrategyList(1);
i2 = StrategyList(2);

% Difficulty levels
Thrsh = .75; % threshold for difficulty
Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i1);
DL1.ii = Lambda_n1<=Thrsh;
Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i2);
DL2.ii = Lambda_n1<=Thrsh;


X  = y1(:,DL1.ii,i1)*100; X=X(:); X=X(1:2:end);
Y  = y2(:,DL1.ii,i1)*100; Y=Y(:); Y=Y(1:2:end);
X2 = y1(:,DL2.ii,i2)*100; X2=X2(:); X2=X2(1:2:end);
Y2 = y2(:,DL2.ii,i2)*100; Y2=Y2(:); Y2=Y2(1:2:end);

Data  = [X , Y ; X2 , Y2];
Label = [zeros(size(X)) ; ones(size(X2))];
SVMModel = fitcsvm(Data,Label,'BoxConstraint',20); 

subplot(sb1,sb2,3)
sv = SVMModel.SupportVectors;
H = gscatter(Data(:,1),Data(:,2),Label);
hold on
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('Pos','Vel','Support Vector')
hold off

%---------------------

[PredictedLabel,Score] = predict(SVMModel,Data);
subplot(sb1,sb2,4)
H = gscatter(Data(:,1),Data(:,2),PredictedLabel);
hold on
hold off
legend('Pos','Vel')





%%
% 
% B = regress( Y , X );
% B2 = regress( Y2 , X2 );
% Y_h =X*B;
% Y_h2 =X2*B2;
% 
% 
% subplot(sb1,sb2,1)
% hold all
% h1 = plot(X,Y,'.','color',CL(i1,:));
% plot(X,Y_h,'-','color',CL(i1,:)/2);
% 
% h2 = plot(X2,Y2,'.','color',CL(i2,:));
% plot(X2,Y_h2,'-','color',CL(i2,:)/2);
% 
% ylim([0,10])
% xlim([0 5])
% text(.5,9,sprintf('a = %.2f',B(1)),'fontsize',11);
% text(4,9,sprintf('a = %.2f',B2(1)),'fontsize',11);
% 
% xlabel('Cursor pos. RMS (cm)')
% ylabel(sprintf('%s \n Cursor vel. RMS (cm/s)',DiffLvl{lm}))
% 
% set(gca,'fontsize',12)
% 




 %%

figure(14)
clf
sb1 = 3;
sb2 = 5;
y1 = D.CursorMeanPos; y1 = y1(:,lmbd_idx,:); % mean over trials
y2 = D.CursorMeanVel; y2 = y2(:,lmbd_idx,:); % mean over trials

[ RcorrModel , PhaseEigRatio ] = deal(NaN(3,3)); % model x difficulty 

txtLbl = {'R_p','R_v','R_p_v'};
for kk=3:-1:1%length(StrategyList)
    
    i = StrategyList(kk);
    % Difficulty levels
    Thrsh = .75; % threshold for difficulty
    Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<1;
    DL(3).ii = Lambda_n1>1;
    
    
    for lm = 1:length(DL) % loop over lambda
        subplot(sb1,sb2,kk+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:2:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:2:end);
        h1(kk) = plot(X,Y,'.','color',CL(i,:));
        plot([-5,5],[0,0],':k');
        plot([0,0],[-2 2],':k');
        ylim([-2,2])
        xlim([-5 5])
        set(gca,'fontsize',12)
        jjj = ~isnan(X);
        if sum(jjj)~=0
            crl = corr(X(jjj),Y(jjj));
            RcorrModel(kk,lm) = crl;
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca([X(jjj),Y(jjj)]);
            PhaseEigRatio(kk,lm) = max(LATENT)/min(LATENT);
            text(-4,1.5,sprintf('%s = %.2f',txtLbl{kk},crl),'fontsize',11);
            text(-4,-1.5,sprintf('ratio = %.2f',max(LATENT)/min(LATENT)),'fontsize',11);
        end
        if kk==1
        xlabel('Mean cursor pos. (cm)')
        ylabel(sprintf('%s \n Mean cursor vel. (cm/s)',DiffLvl{lm}))
        end
    end
    
end
legend(h1,'Pos','Vel','location','southeast')




for kk=2:-1:1%length(StrategyList)
    
    i = StrategyList(kk);
    % Difficulty levels
    Thrsh = .75; % threshold for difficulty
    Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<1;
    DL(3).ii = Lambda_n1>1;
    
    
    for lm = 1:length(DL) % loop over lambda
        subplot(sb1,sb2,5+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:2:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:2:end);
        h1(kk) = plot(X,Y,'.','color',CL(i,:));
        plot([-5,5],[0,0],':k');
        plot([0,0],[-2 2],':k');
        ylim([-2,2])
        xlim([-5 5])
        set(gca,'fontsize',12)
        jjj = ~isnan(X);
        if sum(jjj)~=0
        crl = corr(X(jjj),Y(jjj));
        text(-4,sign(1.5-kk)*1.5,sprintf('%s = %.2f',txtLbl{kk},crl),'fontsize',11);
        end
        if kk==1
        xlabel('Mean cursor pos. (cm)')
        ylabel(sprintf('%s \n Mean cursor vel. (cm/s)',DiffLvl{lm}))
        end
    end
    
end


if SaveFigures
    FLN = sprintf('ModelPhaseMean');
    print(FLN,'-dpng')
end


%%
figure(141)
clf
sb1 = 2;
sb2=12;
for dl=1:3
    subplot(sb1,sb2,2*dl-1)
    hold all
    
    for gr=1:2
        i = StrategyList(gr);
        DD = PhaseEigRatio(gr,dl);
        Y = mean(DD,1);
        S = std(DD,[],1)/sqrt(size(DD,1));
        X = gr;
        %errorbar(X,Y,S,'.k','capsize',.5)
        B = bar(X,Y,.6);
        B.EdgeColor = CL(i,:);
        B.FaceColor = CL(i,:);
        
    end
    if dl==1;  ylabel('PC1 / PC2');   end
    xlim([.5 2.5])
    ylim([0,100])
    
end

for dl=1:3
    subplot(sb1,sb2,sb2+2*dl-1)
    hold all
    
    for gr=1:2
        i = StrategyList(gr);
        DD = RcorrModel(gr,dl);
        Y = mean(DD,1);
        S = std(DD,[],1)/sqrt(size(DD,1));
        X = gr;
        %errorbar(X,Y,S,'.k','capsize',.5)
        B = bar(X,Y,.6);
        B.EdgeColor = CL(i,:);
        B.FaceColor = CL(i,:);
    end
    if dl==1;  ylabel('R');   end
    xlim([.5 2.5])
    ylim([0 1])
end






%%
sb1 = 4;
sb2 = 6;
y1 = D.HandCursDiff_P; y1 = y1(:,lmbd_idx,:); % mean over trials
figure(11)
clf
for kk=1:2%length(StrategyList)
    i = StrategyList(kk);
    for lm = 1:size(y1,2) % loop over lambda
        subplot(sb1,sb2,lm)
        hold all
        h = histogram(y1(:,lm,i),[0:1:40]);
        h.LineStyle = 'none';
        title(sprintf('\\lambda = %.2f',Lambda_List(lm) ) , 'fontweight','normal')
        ylim([0,100])
        xlim([0 40])
        text(20,90-10*kk,sprintf('total: %d \n',sum(h.Values)));
        set(gca,'fontsize',11)    
    end
    
end




sb1 = 4;
sb2 = 6;
y1 = D.CursorMaxVel; y1 = y1(:,lmbd_idx,:); % mean over trials
figure(12)
clf
for kk=1:2%length(StrategyList)
    i = StrategyList(kk);
    for lm = 1:size(y1,2) % loop over lambda
        subplot(sb1,sb2,lm)
        hold all
        h = histogram(y1(:,lm,i),[0:.003:.15]);
        h.LineStyle = 'none';
        title(sprintf('\\lambda = %.2f',Lambda_List(lm) ) , 'fontweight','normal')
        ylim([0,100])
        xlim([0 .15])
        text(20,90-10*kk,sprintf('total: %d \n',sum(h.Values)));
        set(gca,'fontsize',11)
    end
    
end



%%
figure(10)
clf
sb1 = 4;
sb2 = length(StrategyList)+2;
pp = 0;



pp=pp+1;
hold all
y1  = squeeze( nanmean(D.HandRMS_V,1) );  y1 = y1(lmbd_idx,:);% mean over trials
y2  = squeeze( nanmean(D.CursRMS_V,1) );  y2 = y2(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,10],'color',[.7,.7,.7],'linewidth',2)
    hh1 = plot(X,y1(:,i),'.-','color',[.8,.2,.2],'linewidth',1.5,'markersize',10);
    hh2 = plot(X,y2(:,i),'.-','color',[.2,.2,.8],'linewidth',1.5,'markersize',10);
    ylim([0 .1])
    xlim(Limits)
    title(sprintf('%s',Conditions{i}));
    if kk==1; ylabel('RMS_V (m/s)'); end
    xlabel(Label)
    set(gca,'fontsize',11)
end
legend([hh1,hh2],'Hand','Cursor','Locatoin','NorthWest')





y1 = squeeze( nanmean(D.HandCursRMS_v,1) ); y1 = y1(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk+2*sb2)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,10],'color',[.7,.7,.7],'linewidth',2)
    plot(X,y1(:,i),'.-','color',CL(i,:)/1.2,'linewidth',2,'markersize',15)
    ylim([0 .2])
    if kk==1; ylabel('RMS(H-C) Vel'); end
    xlabel(Label)
    xlim(Limits)
    set(gca,'fontsize',11)
end



y1 = squeeze( nanmean(D.HandRMS_V-D.CursRMS_V,1) ); y1 = y1(lmbd_idx,:);% mean over trials
y2 = squeeze( nanmean(D.HandRMS_V,1) ); y2 = y2(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk+3*sb2)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,10],'color',[.7,.7,.7],'linewidth',2)
    H(kk) = plot(X,y1(:,i)./y2(:,i),'.-','color',CL(i,:)/1.2,'linewidth',2,'markersize',15);
    if kk==1; ylabel('1 - RMS[ C / H ]  Vel'); end
    xlabel(Label)
    xlim(Limits)
    ylim([0 1])
    set(gca,'fontsize',11)

end


if SaveFigures
    FLN = sprintf('ModelRMS_V_%s',LCLabels{LambdaCase});
    print(FLN,'-dpng')
end


%%


figure(9)
clf
sb1 = 4;
sb2 = length(StrategyList)+2;
pp = 0;



pp=pp+1;
hold all
y1  = 1000*squeeze( nanmean(D.HandRMS_P,1) );  y1 = y1(lmbd_idx,:); % mean over trials
y2  = 1000*squeeze( nanmean(D.CursRMS_P,1) );  y2 = y2(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
    hh1 = plot(X,y1(:,i),'.-','color',[.8,.2,.2],'linewidth',1.5,'markersize',10);
    hh2 = plot(X,y2(:,i),'.-','color',[.2,.2,.8],'linewidth',1.5,'markersize',10);
    ylim([0 30])
    xlabel(Label)
    xlim(Limits)
    title(sprintf('%s',Conditions{i}));
    if kk==1; ylabel('Pos. RMS (mm)'); end
    set(gca,'fontsize',11)
end
legend([hh1,hh2],'Hand','Cursor','location','best')



y  = squeeze( nanmean(D.HC_XpassCount,1) );  y = y(lmbd_idx,:);% mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk+sb2)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
    plot(X,y(:,i),'.-','color',CL(i,:)/1.2,'linewidth',2,'markersize',15);
    ylim([0 15])
    xlabel(Label)
    xlim(Limits)
    if kk==1; ylabel('XCross#'); end
    set(gca,'fontsize',11)
end



y1 = squeeze( nanmean(D.HandCursRMS_P,1) ); y1 = y1(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk+2*sb2)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
    plot(X,y1(:,i),'.-','color',CL(i,:)/1.2,'linewidth',2,'markersize',15)
    if kk==1; ylabel('RMS(H-C)'); end
    xlabel(Label)
    xlim(Limits)
    ylim([0 .05])
    set(gca,'fontsize',11)
end




y1 = squeeze( nanmean(D.HandRMS_P-D.CursRMS_P,1) ); y1 = y1(lmbd_idx,:); % mean over trials
y2 = squeeze( nanmean(D.HandRMS_P,1) );  y2 = y2(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    subplot(sb1,sb2,kk+3*sb2)
    hold all
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
    H(kk) = plot(X,y1(:,i)./y2(:,i),'.-','color',CL(i,:)/1.2,'linewidth',2,'markersize',15);
    if kk==1; ylabel('1 - RMS [C / H]'); end
    xlabel(Label)
    xlim(Limits)
    ylim([0,1])
    set(gca,'fontsize',11)
end



if SaveFigures
    FLN = sprintf('ModelRMS_P_%s',LCLabels{LambdaCase});
    print(FLN,'-dpng')
end

%%

figure(7)
clf
sb1=2;
sb2=3;
pp=0;

% RMS Cursor Pos (should either remove failed trials or only take the RMS prior to fail)
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y1 = squeeze( nanmean(D.CursRMS_P,1) );  y1 = y1(lmbd_idx,:);% mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y1(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('RMS Cursor Pos. (m)');
ylim([0,.04])



% RMS Cursor Pos (should either remove failed trials or only take the RMS prior to fail)
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y2 = squeeze( nanmean(D.HandRMS_P,1) );  y2 = y2(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('RMS Hand Pos. (m)');
ylim([0,.04])




% RMS difference
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y2 = squeeze( nanmean(D.HandRMS_P./D.CursRMS_P,1) );  y2 = y2(lmbd_idx,:);% mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('Gain ( H/C )');
ylim([0,5])




pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y1 = squeeze( nanmean(D.CursRMS_V,1) );  y1 = y1(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y1(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('RMS Cursor Vel. (m)');
ylim([0 max(y1(:))])



% RMS Hand Vel (should either remove failed trials or only take the RMS prior to fail)
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y1 = squeeze( nanmean(D.HandRMS_V,1) );y1 = y1(lmbd_idx,:);% mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y1(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('RMS Hand Vel. (m)');
ylim([0 max(y1(:))])



% Effort
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
X = Lambda_List; X(1)=[];
y1 = squeeze( nanmean(D.ControlComm,1) ); y1 = y1(lmbd_idx,:);% mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y1(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,100],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylabel('u^Tu');
ylim([0 max(y1(:))])








%%
figure(5)
clf
sb1=3;
sb2=8;
pp=0;



% Lag Position
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
indL = lmbd_idx(3:end-10);
y = squeeze( nanmean(D.Lag_P,1) );  y = y(indL,:); % mean over trials
s = squeeze( nanstd(D.Lag_P,0,1) ); s = s(indL,:);
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',2,'markersize',15)
end
plot([PerfLevel PerfLevel],[0,1000],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([0 400])
ylabel('Lag Pos. (ms)');
set(gca,'fontsize',12)



% Lag Velocity
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Lag_V,1) );  y = y(indL,:); % mean over trials
s = squeeze( nanstd(D.Lag_V,0,1) ); s = s(indL,:);
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[0,1000],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([0 400])
%ylabel('Lag Vel. (ms)');
set(gca,'fontsize',12)
    
    




% Lag Acc
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Lag_A,1) );  y = y(indL,:); % mean over trials
s = squeeze( nanstd(D.Lag_A,0,1) ); s = s(indL,:);
clear HH
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    HH(kk) = plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10);
end
plot([PerfLevel PerfLevel],[0,1000],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([0 200])
%ylabel('Lag Acc. (ms)');
legend(HH,Conditions,'location','best')
set(gca,'fontsize',12)


% Correlation Pos
pp=1+sb2;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Corr_P,1) );  y = y(indL,:); % mean over trials
clear HH
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    HH(kk) = plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',2,'markersize',15);
end
plot([PerfLevel PerfLevel],[-10,10],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([-1 -.7])
ylabel('Correlation');
set(gca, 'YDir','reverse')
set(gca,'fontsize',12)



% Correlation Vel
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Corr_V,1) );  y = y(indL,:); % mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[-10,10],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([-1 0])
%ylabel('Correlation (Vel.)');
set(gca, 'YDir','reverse')
set(gca,'fontsize',11)


% Correlation Acc
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Corr_A,1) ); y = y(indL,:); % mean over trials
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
end
plot([PerfLevel PerfLevel],[-10,10],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([-1 0])
%ylabel('Correlation (Acc.)');
set(gca, 'YDir','reverse')
set(gca,'fontsize',11)


% MI Pos
pp=1+2*sb2;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.MI_P,1) );  y = y(indL,:); % mean over trials
clear HH
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(indL), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    HH(kk) = plot(X,y(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10);
end
plot([PerfLevel PerfLevel],[0,1000],'color',[.7,.7,.7],'linewidth',2)
xlabel(Label)
xlim(Limits)
ylim([0,6])
ylabel('MI');
set(gca,'fontsize',11)




if SaveFigures
    FLN = sprintf('ModelLagCorMI_%s',LCLabels{LambdaCase});
    print(FLN,'-dpng')
end







%% Individual trial

qq  = 2; % control: 1:pos, 2:vel, 3:pv 
fb  = 3; % feedback: 1:pos, 2:vel, 3:pv
lm = length(Lambda_List)*[4/length(Lambda_List),8/length(Lambda_List),11/length(Lambda_List)]; 
tr = 13;


figure(2)
clf
sb1 = 3;
sb2 = 4;
pp=0;
for i=1:length(lm)
    k = lm(i);
    flp = sprintf('%s/Q%d_H%d',FilePath,qq,fb);
    fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
    fprintf('%s \n',fln)
    DD = load(fln);
    Sim = DD.Sim;
    
    C_p = Sim(tr).C_p;
    C_v = Sim(tr).C_v;
    C_a = Sim(tr).C_a;
    H_p = Sim(tr).H_p;
    H_v = Sim(tr).H_v;
    H_a = Sim(tr).H_a;
    C   = squeeze( Sim(1).C(1,1:6,:));
    
    ii = (1:length(C_p))'<Effective_Duration/delta;
    
    
    
    % Puts hold on to add simulation to the figure
    subplot(sb1,sb2,(i-1)*sb2+1)
    hold all
    plot([0,Time(end)],[0,0],':k')
    H1 = plot(Time(ii),C_p(ii),'b','linewidth',2);
    H2 = plot(Time(ii),H_p(ii),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    if i==1
        legend([H1,H2],'Cursor','Hand','location','best')
    end
    ylim([-.058 .058])
    text(.2,.05,sprintf('\\lambda = %.2f',Lambda_List(k)),'FontWeight','bold')
    xlim([0,Effective_Duration]);
    set(gca,'fontsize',11)
    
    
    subplot(sb1,sb2,(i-1)*sb2+2)
    hold all
    [yy,xx] = xcorr(-H_p(ii),C_p(ii)); xx = xx*10; % to ms
    R(i) = corr(C_p(ii),H_p(ii));
    [mm,ind]=max(abs(yy));
    plot(xx,yy,'linewidth',2);
    plot([xx(1),xx(end)],[0,0],':k');
    plot([0 0],mm*[-2,2],':k');
    plot(xx(ind)*[1 1],mm*[0 1],'-k')
    text(2000,mm,sprintf('lag = %.2f \n Corr = %.2f',xx(ind),R(i)));
    xlabel('Lag ms')
    ylabel('XCorr')
    ylim(1.2*mm*[-.5,1])
    Lag_P(i) = xx(ind);
    set(gca,'fontsize',11)
    
    subplot(sb1,sb2,(i-1)*sb2+3)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], max(abs(C_v(ii)))*[-1,1] ,':k');
    plot(C_p(ii),C_v(ii),'color',[.6,.6,.6])
    hh = scatter(C_p(ii),C_v(ii),50,copper(length(C_v(ii)))); hh.Marker = '.';
    colormap(copper)
    colorbar
    xlabel('Cursor P')
    ylabel('Cursor V')
    xlim(.05*[-1.05 1.05])
    set(gca,'fontsize',11)
    
    subplot(sb1,sb2,(i-1)*sb2+4)
    hold all
    plot(.05*[-1,1] , [0,0],':k');
    plot([0,0], .05*[-1,1] ,'-k','linewidth',1);
    plot(.05*[1,-1] , .05*[-1,1] ,'-k','linewidth',1);
    hh = scatter(C_p(ii),H_p(ii),50,copper(length(C_p(ii)))); hh.Marker = '.';
    ylabel('Hand P')
    xlabel('Cursor P')
    ylim(max(abs(H_p(ii)))*[-1 1])
    xlim(max(abs(C_p(ii)))*[-1 1])
    axis equal
    set(gca,'fontsize',11)

end




%%
figure(1)
clf
sb1 = 3;
sb2 = 4;
pp=0;
for i=1:length(lm)
    k = lm(i);
    flp = sprintf('%s/Q%d_H%d',FilePath,qq,fb);
    fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
    fprintf('%s \n',fln)
    DD = load(fln);
    Sim = DD.Sim;
    
    C_p = Sim(tr).C_p;
    C_v = Sim(tr).C_v;
    C_a = Sim(tr).C_a;
    H_p = Sim(tr).H_p;
    H_v = Sim(tr).H_v;
    H_a = Sim(tr).H_a;
    C   = squeeze( Sim(1).C(1,1:6,:));
    K   = Sim(1).K(1:6,:);
    u   = Sim(tr).u;
    
    ii = (1:length(C_p))'<Effective_Duration/delta;
    
    % Puts hold on to add simulation to the figure
    subplot(sb1,sb2,(i-1)*sb2+1)
    hold all
    plot([0,Time(end)],[0,0],':k')
    H1 = plot(Time(ii),C_p(ii),'b','linewidth',2);
    H2 = plot(Time(ii),H_p(ii),'r','linewidth',2);
    xlabel('Time (s)')
    ylabel('Position')
    if i==1
        legend([H1,H2],'Cursor','(-) Hand','location','best')
    end
    ylim([-.058 .058])
    text(.2,.05,sprintf('\\lambda = %.2f',Lambda_List(k)),'FontWeight','bold')
    xlim([0,Effective_Duration]);
    set(gca,'fontsize',11)
    
    subplot(sb1,sb2,(i-1)*sb2+2)
    hold all
    plot(Time(ii),u(ii),'linewidth',2)
    plot([0,Time(end)],[0,0],':k')
    ylabel('Control input: u')
    xlabel('Time (s)')
    xlim([0,Effective_Duration]);
    set(gca,'fontsize',11)
    
    subplot(sb1,sb2,(i-1)*sb2+3)
    hold all
    plot(Time(ii),C(1:2,(ii))','linewidth',2);
    ylabel('Control Gain: L')
    xlabel('Time (s)')
    if i==1
    legend('Pos.','Vel','location','best')
    end
    xlim([0,Effective_Duration]);
    set(gca,'fontsize',11)
    
    subplot(sb1,sb2,(i-1)*sb2+4)
    hold all
    plot(Time(ii),K(1:2,ii)','linewidth',2);
    ylabel('Kalman Gain: K')
    xlabel('Time (s)')
    if i==1
    legend('Pos.','Vel','location','best')
    end
    xlim([0,Effective_Duration]);
    set(gca,'fontsize',11)
end



%% 

for qq =1:2
    figure(1000+qq)
    clf
    sb1 = 4;
    sb2 = 5;
    fb = 3; % feedback type, both pos and vel
    k = 8; % lambda index
    
    flp = sprintf('%s/Q%d_H%d',FilePath,qq,fb);
    fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
    fprintf('%s \n',fln)
    DD = load(fln);
    Sim = DD.Sim;
    
    pp=0;
    for tr=1:500
        
        C_p = Sim(tr).C_p * 100;
        H_p = Sim(tr).H_p * 100;
        ii = (1:length(C_p))'<Effective_Duration/delta;
        
        % check for success
        jjj = abs(C_p(ii))>5;
        if sum(jjj)~=0
            continue
        end
        
        % Puts hold on to add simulation to the figure
        pp = pp+1;
        if pp>sb1*sb2; break; end
        subplot(sb1,sb2,pp)
        hold all
        plot([0,Time(end)],[0,0],':k')
        H1 = plot(Time(ii),C_p(ii),'b','linewidth',2);
        H2 = plot(Time(ii),H_p(ii),'r','linewidth',2);
        xlabel('Time (s)')
        ylabel('Cursor Pos (cm)')
        if i==1
            legend([H1,H2],'Cursor','Hand','location','best')
        end
        ylim([-5 5])
        text(.2,5,sprintf('\\lambda = %.2f',Lambda_List(k)),'FontWeight','bold')
        text(.2,-4,sprintf('C_{mn} = %.2f',mean(C_p(ii))))
        text(3,-4,sprintf('C_{rms} = %.2f',rms(C_p(ii))))
        xlim([0,Effective_Duration]);
        set(gca,'fontsize',9)
        
    end
   
end




%%

for qq =1:2
    figure(1100+qq)
    clf
    sb1 = 4;
    sb2 = 5;
    fb = 3; % feedback type, both pos and vel
    
    k = 8; % lambda index
    flp = sprintf('%s/Q%d_H%d',FilePath,qq,fb);
    fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
    fprintf('%s \n',fln)
    DD = load(fln);
    Sim = DD.Sim;
    
    pp=0;
    for tr=1:200
        
        
        
        C_p = Sim(tr).C_p * 100;
        H_p = Sim(tr).H_p * 100;
        
        C_v = Sim(tr).C_v * 100;
        H_v = Sim(tr).H_v * 100;
        
        
        ii = (1:length(C_p))'<Effective_Duration/delta;
        
        % check for success
        jjj = abs(C_p(ii))>5;
        if sum(jjj)~=0
            continue
        end
        
        % Puts hold on to add simulation to the figure
        pp = pp+1;
        if pp>sb1*sb2; break; end
        subplot(sb1,sb2,pp)
        hold all
        plot([0,Time(end)],[0,0],':k')
        H1 = plot(Time(ii),C_v(ii),'b','linewidth',2);
        H2 = plot(Time(ii),C_p(ii),':r','linewidth',2);
        xlabel('Time (s)')
        ylabel('Cursor Vel. (cm/s)')
        if tr==1
            legend([H1,H2],'Pos','Vel','location','best')
        end
        ylim([-10 10])
        text(.2,-7,sprintf('Cv_{mn} = %.2f',mean(C_v(ii))))
        text(3,-7,sprintf('Cv_{rms} = %.2f',rms(C_v(ii))))
        xlim([0,Effective_Duration]);
        set(gca,'fontsize',9)
        
    end
    
end








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L_out, Limits , Label , PerfLevel] = normLambda(L_in, L_c , Sigma_c , Case)

% Convert lambda values to normalized values
%                  Inputs
% -------------------------------------------
% L_in: raw lambda data
% L_c: critical lambda
% Sigma_c: critical slope
% Case: normalization case
%                  Outputs
% -------------------------------------------
% L_out: converted/normalized lambda
% Limits: xlim for the converted/normalized lambda
% Label: the label for new lambda
% PerfLevel: value of normalized lambda at 50% success

switch Case
    
    case 1 % no normalization
        L_out = L_in;
        Limits = [1 6];
        Label = '\lambda';
        PerfLevel = NaN;
    case 2 % normalize by critical lambda
        L_out = L_in/L_c;
        Limits = [0,1.5];
        Label = '\lambda / \lambda_c';
        PerfLevel = 1;
    case 3 % normalize by critical lambda and critical slope
        L_out = (L_in-L_c)/Sigma_c;
        Limits = [-10,5];
        Label = '\lambda_N';
        PerfLevel = 0;
    case 4 % convert to success rate
        L_N = (L_in-L_c)/Sigma_c;
        L_out = 50*(1-erf(L_N/sqrt(2)));
        Limits = [0,100];
        Label = '% Success';
        PerfLevel = 50;
end

end
















%% 
% 
% sb1 = 3;
% sb2 = 4;
% y1  = D.CursRMS_P; y1 = y1(:,lmbd_idx,:); % mean over trials
% y2  = D.CursRMS_V; y2 = y2(:,lmbd_idx,:); % mean over trials
% 
% figure(16)
% clf
% for kk=1:2%length(StrategyList)
%     i = StrategyList(kk);
%     for lm = 1:min(sb1*sb2,size(y1,2))  % loop over lambda
%         subplot(sb1,sb2,lm)
%         hold all
%         Slope = y1(:,lm,i)./y2(:,lm,i);
%         h1(kk) = plot(Slope,'.','color',CL(i,:));
%         plot([0 sum(~isnan(Slope))],nanmean(Slope)*[1,1],'-' , 'color',CL(i,:)/2)
%         %plot(y1(:,lm,i) ,Lambda_List(lm)*y1(:,lm,i),'.k' )
%         
%         ylim([0,2])
%         %xlim([0 .04])
%         set(gca,'fontsize',12)
%         %xlabel('RMS_p (m)')
%         %ylabel('RMS_v (m/s)')
%         jjj = ~isnan(y1(:,lm,i));
%         if sum(jjj)~=0
%         crl = corr(y1(jjj,lm,i),y2(jjj,lm,i));
%         %text(.03,.05,sprintf('\\lambda = %.2f',Lambda_List(lm)));
%         %text(-.038,.012-kk*.003,sprintf('R = %.2f',crl));
%         end
%     end
%     
% end
% legend(h1,'C_p','C_v')
% 
% 
% 





%%  __________________________________/ Control Gain: Cursor 
% figure(4)
% clf
% sb1=3;
% sb2=4;
% pp=0;
% 
%
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+1)
%     hold all
%     y2 = squeeze( nanmean(D.ControlGain(:,j,:,:),1) );  y2 = y2(lmbd_idx,:); % mean over trials
%     plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
%     for i=1:size(y2,2)
%     X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     %xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('L:  Pos.');
%         title('Cursor')
%     elseif j==2
%          ylabel('L:  Vel.');
%     else
%          ylabel('L:  Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% % Control Gain: Hand 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+2)
%     hold all
%     y2 = squeeze( nanmean(D.ControlGain(:,3+j,:,:),1) );  y2 = y2(lmbd_idx,:);  % mean over trials
%     plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
%     for i=1:size(y2,2)
%     X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     %xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('L: Pos.');
%         title('Hand')
%     elseif j==2
%          ylabel('L: Vel.');
%     else
%          ylabel('L: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% 
% 
% 
% 
%% _________________________________________/ Kalman Gain: Cursor 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+3)
%     hold all
%     y2 = squeeze( nanmean(D.KalmanGain(:,j,:,:),1) );  y2 = y2(lmbd_idx,:);  % mean over trials
%     plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
%     for i=1:size(y2,2)
%     X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     %xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('K: Pos.');
%         title('Cursor')
%     elseif j==2
%          ylabel('K: Vel.');
%     else
%          ylabel('K: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% % Kalman Gain: Hand 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+4)
%     hold all
%     y2 = squeeze( nanmean(D.KalmanGain(:,3+j,:,:),1) );  y2 = y2(lmbd_idx,:); % mean over trials
%     plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
%     for i=1:size(y2,2)
%     X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     %xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('K: Pos.');
%         title('Hand')
%     elseif j==2
%          ylabel('K: Vel.');
%     else
%          ylabel('K: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% legend(Conditions,'location','best')




%% Control Gain initial and end
% figure(3)
% clf
% sb1=3;
% sb2=4;
% pp=0;
% 
% 
% % Control Gain: Cursor 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+1)
%     hold all
%     X = Lambda_List; X(1)=[];
%     Z = D.ControlGain_I(:,j,:,:)-D.ControlGain(:,j,:,:);
%     y2 = squeeze( nanmean(Z,1) );  y2(1,:)=[]; % mean over trials
%     for i=1:size(y2,2)
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('L:  Pos.');
%         title('Cursor')
%     elseif j==2
%          ylabel('L:  Vel.');
%     else
%          ylabel('L:  Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% % Control Gain: Hand 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+2)
%     hold all
%     Z = D.ControlGain_I(:,j+3,:,:)-D.ControlGain(:,j+3,:,:);
%     y2 = squeeze( nanmean(Z,1) );  y2(1,:)=[]; % mean over trials
%     for i=1:size(y2,2)
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('L: Pos.');
%         title('Hand')
%     elseif j==2
%          ylabel('L: Vel.');
%     else
%          ylabel('L: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% 
% 
% 
% 
% % Control Gain: Cursor 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+3)
%     hold all
%     Z = D.ControlGain_I(:,j,:,:)-D.ControlGain_F(:,j,:,:);
%     y2 = squeeze( nanmean(Z,1) );  y2(1,:)=[]; % mean over trials
%     for i=1:size(y2,2)
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('K: Pos.');
%         title('Cursor')
%     elseif j==2
%          ylabel('K: Vel.');
%     else
%          ylabel('K: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% 
% 
% % Control Gain: Cursor 
% for j=1:3
%     subplot(sb1,sb2,sb2*(j-1)+4)
%     hold all
%     Z = D.ControlGain_I(:,j+3,:,:)-D.ControlGain_F(:,j+3,:,:);
%     y2 = squeeze( nanmean(Z,1) );  y2(1,:)=[]; % mean over trials
%     for i=1:size(y2,2)
%         plot(X,y2(:,i),'.-','color',CL(i,:),'linewidth',1.5,'markersize',10)
%     end
%     xlim([0,Lambda_List(end)+1]) 
%     if j==1
%         ylabel('K: Pos.');
%         title('Hand')
%     elseif j==2
%          ylabel('K: Vel.');
%     else
%          ylabel('K: Acc.');   
%     end
%            
%     xlabel('\lambda')
% end
% legend(Conditions,'location','best')
% 
% 

