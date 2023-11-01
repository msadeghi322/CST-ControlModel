



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




%% Analysis

FileName = 'sim';


Case = 1;

switch Case
    
    case 1
        
%         FilePath = 'Reduced data/CST_R1.2_Sensitivity_Delay01';
%         StrategyList = [1,2,3]; 
%         DelayList      = [20, 50, 100]/1000;
        
        FilePath = 'Reduced data/CST_R1.2_Sensitivity_Delay02';
        StrategyList = [1,2,3]; 
        DelayList      = [30, 50, 70]/1000;
        

        Obs_List = [1 1 1 1 1];
        Q_List = [1e5 0 0 0 0 
                  0   1e10 0 0 0];

        fpth = sprintf('%s/R1',FilePath);
        Lambda_List = [.5:.4:6.5]';
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
        for j=1:size(DelayList,2)
            Cond = Cond+1;
            
            for k=1:length(Lambda_List)
                
                flp = sprintf('%s/Q%d_H%d',FilePath,i,j);
                fln = sprintf('%s/%s_L%d.mat',flp,FileName,k);
                fprintf('%s \n',fln)
                D = load(fln);
                Sim = D.Sim;
                clear D;
                FN = fieldnames(Sim);
                FN(13:end)=[];
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
                    
                    
                    
                    % calculate correlation
                    %Corr_P(n,k,Cond) = corr(C_p(ii),H_p(ii));
                    %Corr_V(n,k,Cond) = corr(C_v(ii),H_v(ii));
                    %Corr_A(n,k,Cond) = corr(C_a(ii),H_a(ii));
                    rp = xcorr(C_p(ii)-mean(C_p(ii)),H_p(ii)-mean(H_p(ii)),'coeff');
                    Corr_P(n,k,Cond) = min(rp);
                    rv = xcorr(C_v(ii)-mean(C_v(ii)),H_v(ii)-mean(H_v(ii)),'coeff');
                    Corr_V(n,k,Cond) = min(rv);
                    
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



%% Initialize

% Lambda case: plotting agaist:
%_____ Lambda:                     Case=1
%_____ Lambda/Lambda_c:            Case=2
%_____ (Lambda-Lambda_c)/Sigma_c:  Case=3;
%_____ Success rate:               Case=4

LambdaCase = 1;
LCLabels = {'L','LC','LN','Success'};

SaveFigures = 0;
% GroupColor = [.5,.4,0
%       0,.6,.6
%       .8,.5,.5
%       .6,0,0
%       .5,.8,.5;
%       0,.4,0
%       .5,.5,.8
%       0,0,.7
%       .6,0,.9
%       ];
GroupColor = flip(jet(18));
GroupColor = GroupColor(1:3:end, :);


FontSize = 11;  
LineWidth = 1.5;
MarkerSize = 12;

%% 

figure(101)
clf
sb1=4;
sb2=9;
pp=0;

ii = Lambda_List<5.5 & Lambda_List>=1.7;

lmbd_idx = ii;%1:length(Lambda_List);
Lambda_c = NaN(9,1);
Sigma_c = NaN(9,1);

% Success rate ______________________________________
pp=pp+1;
subplot(sb1,sb2,pp)
hold all
X  = Lambda_List(lmbd_idx); 
y  = squeeze( nanmean(D.Success,1) );  y = y(lmbd_idx,:); % mean over trials
ft=fittype('50*( 1-erf((x - mu)/sqrt(2)/sigma) )', 'independent', 'x', 'dependent', 'y' );
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    B0 = [-2,4,4];
    FF = fit(X,100*y(:,i),ft, 'StartPoint',[3.5,0.2],'Robust','Bisquare');
    l = 0:.01:Lambda_List(end);
    yf = FF(l);
    %[~,ind] = min(abs(yf-.5)); % find critical lambda
    Lambda_c(i) = FF.mu;
    Sigma_c(i)  = FF.sigma;
    hh(i) = plot(l,yf,'color',GroupColor(i,:),'linewidth',LineWidth);
    plot(X,100*y(:,i),'.','color',GroupColor(i,:),'linewidth',LineWidth,'markersize',MarkerSize)
    a=1;
end
plot([0,Lambda_List(end)],[50,50],':k')
ylim([-5 105])
%xlim([0,Lambda_List(end)+1])
xlim([1,6])
ylabel('% Success')
xlabel('\lambda')
set(gca,'fontsize',FontSize)
%legend(hh,Conditions,'location','best')


% Lag Position ______________________________________________
pp=pp+sb2;
subplot(sb1,sb2,pp)
hold all
indL = lmbd_idx(3:end-10);
y = squeeze( nanmean(D.Lag_P,1) );  y = y(lmbd_idx,:); % mean over trials
s = squeeze( nanstd(D.Lag_P,0,1) ); s = s(lmbd_idx,:);
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',GroupColor(i,:),'linewidth',LineWidth,'markersize',MarkerSize)
end
plot([PerfLevel PerfLevel],[0,1000],'color',[.7,.7,.7],'linewidth',LineWidth)
xlabel(Label)
xlim(Limits)
ylim([0 400])
ylabel('Lag (ms)');
set(gca,'fontsize',FontSize)


% Correlation Pos
pp=pp+sb2;
subplot(sb1,sb2,pp)
hold all
y = squeeze( nanmean(D.Corr_P,1) );  y = y(lmbd_idx,:); % mean over trials
clear HH
for kk=1:length(StrategyList)
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    HH(kk) = plot(X,y(:,i),'.-','color',GroupColor(i,:),'linewidth',LineWidth,'markersize',MarkerSize);
end
plot([PerfLevel PerfLevel],[-10,10],'color',[.7,.7,.7],'linewidth',LineWidth)
xlabel(Label)
xlim(Limits)
ylim([-1 -.7])
ylabel('Correlation');
set(gca, 'YDir','reverse')
set(gca,'fontsize',FontSize)



pp=pp+sb2;
subplot(sb1,sb2,pp)
hold all
G = D.HandRMS_P ./ D.CursRMS_P;
y  = squeeze( nanmean(G,1) );  y = y(lmbd_idx,:); % mean over trials
for kk=1:length(StrategyList)
    
    i = StrategyList(kk);
    [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
    plot(X,y(:,i),'.-','color',GroupColor(i,:),'linewidth',LineWidth,'markersize',MarkerSize);
    
end
plot([PerfLevel PerfLevel],[-10,10],'color',[.7,.7,.7],'linewidth',LineWidth)
xlabel(Label)
xlim(Limits)
ylim([.8 1.6])
ylabel('RMS Ratio');
set(gca,'fontsize',FontSize)
set(gca,'ytick',[.8:.4:1.6])


% Critical lambda based on % success
% pp=pp+1;
% subplot(sb1,sb2,sb2+1)
% hold all
% for kk=1:length(StrategyList)
%     i = StrategyList(kk);
%     bb = bar(2*kk-1,Lambda_c(i),.7);
%     bb.FaceColor = CL(i,:);
% end
% set(gca,'xtick',2*[1:9]-1,'xticklabel',Conditions(StrategyList))
% xtickangle(70)
% ylim([0,Lambda_List(end)+1]);
% ylabel('\lambda_C')
% %legend(Conditions(StrategyList),'location','best','orientation','horizontal')
% grid
% set(gca,'fontsize',12)
%%%%%%%%%%%%%%%%%%%%%%%%%% Corr and 



if SaveFigures
    FLN = sprintf('ModelLagCorMI_%s',LCLabels{LambdaCase});
    print(FLN,'-dpng')
end





%%

idx = 1:length(Lambda_List);

y1 = D.CursRMS_P; y1 = y1(:,idx,:); % mean over trials
y2 = D.CursRMS_V; y2 = y2(:,idx,:); % mean over trials

figure(15)
clf
sb1 = 4;
sb2 = 9;
DiffLvl = {'Easy','Moderate','Hard'};
SkipData = 2; % skip every other
%____________________________________/ Regression without offset
PhaseSlope = [];
for kk=length(StrategyList):-1:1%
    
    i = StrategyList(kk);
    % Difficulty levels
    Thrsh = .9; % threshold for difficulty
    Lambda_n1 = Lambda_List(idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<1;
    DL(3).ii = Lambda_n1>1;
    
    for lm = 1:1%length(DL)  % loop over lambda
        subplot(sb1,sb2,kk+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:SkipData:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:SkipData:end);
        B = regress( Y , X );
        PhaseSlope(kk,lm) = B(1);
        Y_h =X*B; 
        h1(kk) = plot(X,Y,'o','MarkerSize',MarkerSize-11);
        h1(kk).MarkerFaceColor = GroupColor(i,:);
        h1(kk).MarkerEdgeColor = GroupColor(i,:);
        
        [XX,sri]=sort(X);
        plot(XX,Y_h(sri),'-','color',[.2 .2 .2],'linewidth',LineWidth);
        ylim([0,10])
        xlim([0 5])
        text(.5,9,sprintf('a = %.2f',B(1)),'fontsize',FontSize);
        if kk==1
        %xlabel('Cursor pos. RMS (cm)')
        ylabel(sprintf('%s \n RMS Cursor Vel. (cm/s)',DiffLvl{lm}))
        end
        set(gca,'fontsize',FontSize)
        
        
        subplot(sb1,sb2,5+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:SkipData:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:SkipData:end);
        B = regress( Y , X );
        PhaseSlope(kk,lm) = B(1);
        Y_h =X*B; 
        h1(kk) = plot(X,Y,'o','MarkerSize',MarkerSize-11);
        h1(kk).MarkerFaceColor = GroupColor(i,:);
        h1(kk).MarkerEdgeColor = GroupColor(i,:);
        
        [XX,sri]=sort(X);
        plot(XX,Y_h(sri),'-','color',[.2 .2 .2],'linewidth',LineWidth);
        ylim([0,10])
        xlim([0 5])
        set(gca,'fontsize',FontSize)
        
        
    end    
end
%legend(h1,'Pos','Vel','location','southeast')



 %%

figure(14)
clf
sb1 = 4;
sb2 = 9;
y1 = D.CursorMeanPos; y1 = y1(:,lmbd_idx,:); % mean over trials
y2 = D.CursorMeanVel; y2 = y2(:,lmbd_idx,:); % mean over trials

[ RcorrModel , PhaseEigRatio ] = deal(NaN(3,3)); % model x difficulty 

txtLbl = {'R_p','R_v','R_p_v'};
for kk=length(StrategyList):-1:1%
    
    i = StrategyList(kk);
    % Difficulty levels
    Lambda_n1 = Lambda_List(lmbd_idx)./Lambda_c(i);
    DL(1).ii = Lambda_n1<=Thrsh;
    DL(2).ii = Lambda_n1>Thrsh & Lambda_n1<Thrsh;
    DL(3).ii = Lambda_n1>Thrsh;
    
    
    for lm = 1:1%length(DL) % loop over lambda
        
        subplot(sb1,sb2,kk+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:SkipData:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:SkipData:end);
        h1(kk) = plot(X,Y,'.','color',GroupColor(i,:));
        plot([-5,5],[0,0],':k');
        plot([0,0],[-2 2],':k');
        ylim([-2,2])
        xlim([-5 5])
        set(gca,'fontsize',FontSize)
        jjj = ~isnan(X);
        if sum(jjj)~=0
            crl = corr(X(jjj),Y(jjj));
            RcorrModel(kk,lm) = crl;
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca([X(jjj),Y(jjj)]);
            PhaseEigRatio(kk,lm) = max(LATENT)/min(LATENT);
            %text(-4,1.5,sprintf('%s = %.2f',txtLbl{kk},crl),'fontsize',FontSize);
            text(-4,-1.5,sprintf('R = %.2f',crl),'fontsize',FontSize);
        end
        if kk==1
        %xlabel('Mean cursor pos. (cm)')
        ylabel(sprintf('%s \n Mean Cursor Vel. (cm/s)',DiffLvl{lm}))
        end
    end
    
end
%legend(h1,'Pos','Vel','location','southeast')




for kk=2:-1:1%length(StrategyList)
    
    i = StrategyList(kk);
    % Difficulty levels
    for lm = 1:length(DL) % loop over lambda
        subplot(sb1,sb2,5+(lm-1)*sb2)
        hold all
        X = y1(:,DL(lm).ii,i)*100; X=X(:); X=X(1:SkipData:end);
        Y = y2(:,DL(lm).ii,i)*100; Y=Y(:); Y=Y(1:SkipData:end);
        plot(X,Y,'.','color',GroupColor(i,:));
        plot([-5,5],[0,0],':k');
        plot([0,0],[-2 2],':k');
        ylim([-2,2])
        xlim([-5 5])
        set(gca,'fontsize',FontSize)
        jjj = ~isnan(X);
        if sum(jjj)~=0
        crl = corr(X(jjj),Y(jjj));
        %text(-4,sign(1.5-kk)*1.5,sprintf('%s = %.2f',txtLbl{kk},crl),'fontsize',FontSize);
        end
        if kk==1
        %xlabel('Mean cursor pos. (cm)')
        %ylabel(sprintf('%s \n Mean cursor vel. (cm/s)',DiffLvl{lm}))
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
        B.EdgeColor = GroupColor(i,:);
        B.FaceColor = GroupColor(i,:);
        
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
        B.EdgeColor = GroupColor(i,:);
        B.FaceColor = GroupColor(i,:);
    end
    if dl==1;  ylabel('R');   end
    set(gca,'xtick',[1,2],'xticklabel',{'Pos','Vel'})
    xlim([.5 2.5])
    ylim([0 1])
end













%%

for qq =1:2
    figure(1000+qq)
    clf
    sb1 = 4;
    sb2 = 5;
    fb = 1; % feedback type, both pos and vel
    k = 10; % lambda index
    
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



%%  __________________________________/ Control Gain: Cursor 
figure(13)
clf
sb1=3;
sb2=4;
pp=0;


for j=1:2
    subplot(sb1,sb2,j)
    hold all
    y2 = squeeze( nanmean(D.ControlGain(:,j,:,:),1) );  y2 = y2(lmbd_idx,:); % mean over trials
    %plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
    for kk=1:length(StrategyList)
        i = StrategyList(kk);
        [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
        %X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
        plot(X,y2(:,i),'.-','color',GroupColor(i,:),'linewidth',1.5,'markersize',10)
    end
    %xlim([0,Lambda_List(end)+1]) 
    if j==1
        ylabel('Position Gain');
        title('Cursor')
    elseif j==2
         ylabel('Velocity Gain');
    else
         ylabel('Acceleration Gain');   
    end
    xlim(Limits)   
    ylim([0,60])
    xlabel('\lambda')
    set(gca,'fontsize',FontSize)
end


% Control Gain: Hand 
for j=1:2
    subplot(sb1,sb2,2*sb2+j)
    hold all
    y2 = squeeze( nanmean(D.ControlGain(:,2+j,:,:),1) );  y2 = y2(lmbd_idx,:);  % mean over trials
    plot([1,1],[0 max(y2(:))],'color',[.7,.7,.7],'linewidth',2)
    for kk=1:length(StrategyList)
        i = StrategyList(kk);
        [X, Limits , Label , PerfLevel] = normLambda(Lambda_List(lmbd_idx), Lambda_c(i) , Sigma_c(i) , LambdaCase);
        %X = Lambda_List(lmbd_idx)/Lambda_nrm(i);
        plot(X,y2(:,i),'.-','color',GroupColor(i,:),'linewidth',1.5,'markersize',10)
    end
    %xlim([0,Lambda_List(end)+1]) 
    if j==1
        ylabel('Position Gain');
        title('Hand')
    elseif j==2
         ylabel('Velocity Gain');
    else
         ylabel('Acceleration Gain');   
    end
    xlim(Limits)   
    ylim([0,60])      
    xlabel('\lambda')
    set(gca,'fontsize',FontSize)
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












