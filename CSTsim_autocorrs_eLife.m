% Autocorrelation computations for eLife paper with Mohsen et al.
clear 
clc
close all
set(0,'DefaultFigureWindowStyle','docked') 
set(0,'defaultAxesFontSize',15)
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
set(groot,'defaultAxesBox','off')
set(0, 'DefaultFigureRenderer', 'painters');


load('SimData.mat')
% Data structure info (from Mohsen's email):
% SimData( i , j ) -> control objective i (i=1 for pos and i=2 for vel), and lambda index j. The structures contains:
% Lambda: scalar 
% C_p (Cursor position): 801-by-50(columns are different trials)
% C_v (Cursor velocity)
% H_p (Hand position)
% H_v (Hand velocity)

% Let's focus on 1 Lambda value and a few trials each of position &
% velocity to start with

fs=100; % sampling freq [Hz]
Nsamples=801; % length of simulated cursor and hand vectors per trial
t=[0:Nsamples-1]/fs;
idxP=1; idxV=2; idxL=4; Ntrials=4;
Lamb=SimData(idxP,idxL).Lambda; 
CpP=SimData(idxP,idxL).C_p;
HpP=SimData(idxP,idxL).H_p;
CpV=SimData(idxV,idxL).C_p;
HpV=SimData(idxV,idxL).H_p;


acf_idx=[1:2*Nsamples-1]-Nsamples; % acf sample index, with peak at 0 
figure('Position',[25 10 1400 800]);nfig=1;
for nn=1:Ntrials
    subplot(Ntrials,4,nfig)
    trl = randi(50);
    plot(t,CpP(:,trl),'b',t,HpP(:,trl),'r','linewidth',2),grid on;
    if nn==1
        title('Position Control'); 
        ylabel('Position (m)')
    end
    if nn==Ntrials
        xlabel('Time (s)'); 
        legend('Cursor','Hand','Location','southwest'); 
    end
    xlim([0,8])
    ylim([-.05,.05])
    
    
    
    
    subplot(Ntrials,4,nfig+1)
    acf_CpP=xcorr(HpP(:,trl),'normalized');
    rect_width=fix(sum(abs(acf_CpP)));
    equiv_rect=stepfun(acf_idx,-floor(rect_width/2))-stepfun(acf_idx,ceil(rect_width/2));
    clip_acf=diff(acf_CpP.*(acf_CpP>0.5));
    [~,mainlobewidth]=min(clip_acf);
%     std=sqrt(sum((acf_idx).^2.*abs(acf_HpP').^2/sum(abs(acf_HpP').^2)));
%     std=sqrt(sum(acf_idx.^2.*acf_HpP'.^2)/sum(abs(acf_HpP').^2));
    plot(acf_idx,acf_CpP,'k','linewidth',2),
    grid on;xlim([-(Nsamples-1) Nsamples])
    total_intensity=sum(abs(acf_CpP(Nsamples:end)));
    cumul_intensity=cumsum(abs(acf_CpP(Nsamples:end)))/total_intensity;
    hold,plot(acf_idx,equiv_rect,'k:','linewidth',2);
%     hold,plot(acf_idx(Nsamples:end),cumul_intensity,'k:','linewidth',2)
%     cumul_intensity=cumul_intensity.*(cumul_intensity<0.5);
    [~,width]=max(cumul_intensity.*(cumul_intensity<0.5));
    txt = {['0.5 width=',num2str(2*(mainlobewidth-Nsamples),3)],...
        ['rect width=',num2str(rect_width,3)]};
    text(-700,-.7,txt)
    ylim([-1,1.5])
    
    
    
    
    if nn==1,title('acf Cpos'),end
    subplot(Ntrials,4,nfig+2)
    plot(t,CpV(:,trl),'b',t,HpV(:,trl),'r','linewidth',2),grid on;
    if nn==1,title('Velocity Control'); ylabel('Position (m)'); ;end
    xlim([0,8])
    ylim([-.05,.05])
    if nn==Ntrials
        xlabel('Time'); 
        legend('Cursor','Hand','Location','southwest'); 
    end
    
    
    
    subplot(Ntrials,4,nfig+3)
    acf_CpV=xcorr(HpV(:,trl),'normalized');
    rect_width=fix(sum(abs(acf_CpV))/2);
    equiv_rect=stepfun(acf_idx,-floor(rect_width/2))-stepfun(acf_idx,ceil(rect_width/2));
    clip_acf=diff(acf_CpV.*(acf_CpV>0.5));
    [~,mainlobewidth]=min(clip_acf);
%     std=sqrt(sum(acf_idx.^2.*abs(acf_HpV').^2)/sum(abs(acf_HpV').^2));
    plot(acf_idx,acf_CpV','k','linewidth',2)
    grid on;
    xlim([-(Nsamples-1) Nsamples]);
    total_intensity=sum(abs(acf_CpV(Nsamples:end)));
    cumul_intensity=cumsum(abs(acf_CpV(Nsamples:end)))/total_intensity;
    hold,plot(acf_idx,equiv_rect,'k:','linewidth',2);
%     cumul_intensity=cumul_intensity.*(cumul_intensity<0.5);
%     hold,plot(acf_idx(Nsamples:end),cumul_intensity,'k:','linewidth',2)
    [~,width]=max(cumul_intensity.*(cumul_intensity<0.5));
    txt = {['0.5 width=',num2str(2*(mainlobewidth-Nsamples),3)],...
        ['rect width=',num2str(rect_width,3)]};
    if nn==1,title('acf Cpos'),end
    nfig=nfig+4;
    text(-700,-.7,txt)
    ylim([-1,1.5])
end

%% %%% Let's take a crack at looping through all data  %%%%%%%
fig1=figure('Position',[25 10 1400 800]);hold on;
Nstrategy=2; Nlambda=4; Ntrials=50;
rect_width=zeros(Nstrategy,Nlambda,Ntrials);
mainlobe_width=zeros(Nstrategy,Nlambda,Ntrials);
Lambda=zeros(1,Nlambda);
for ii=1:Nstrategy % position=1, velocity=2
    for jj=1:Nlambda 
        Lambda(jj)=SimData(1,jj).Lambda;
        for nn=1:Ntrials 
            Cp=SimData(ii,jj).C_p(:,nn);
            Hp=SimData(ii,jj).H_p(:,nn);
            acf_Cp=xcorr(Cp,'normalized');
            clip_acf=-diff(acf_Cp.*(acf_Cp>0.5)); % create spikes where acf crosses 0.5
            clip_acf=clip_acf.*(clip_acf>0.25); % keep only pos spikes > 0.25
            [~,locs]=findpeaks(clip_acf); % find indices of spikes
            pospeaklocs=nonzeros(locs.*(locs>Nsamples)); % spike indices > Nsamples
%             [~,idx]=min(diff(acf_Cp.*(acf_Cp>0.5)));
%             if idx<Nsamples  % debug -- find acf's with negative width
%                 ii,jj,nn
%             end
            mainlobe_width(ii,jj,nn)=2*(pospeaklocs(1)-Nsamples)+1;
            rect_width(ii,jj,nn)=sum(abs(acf_Cp));
            subplot(2,1,ii),hold on; plot(acf_Cp)
        end
    end
    if ii==1,title('Pos Control'); else, title('Vel Control'); end
    grid on;
end
HistoBinWidth=50;
figure('Position',[35 10 1400 800]);title('ACF width');

for i=1:length(Lambda)
subplot(2,2,i)
histogram(rect_width(1,i,:),'BinWidth',HistoBinWidth);
hold on; histogram(rect_width(2,i,:),'BinWidth',HistoBinWidth);
legend('Pos','Vel');
xlabel('Width [samples]'),ylabel('Count')
title(sprintf('acf Cp rect width for Pos vs. Vel Control \n Lamba=%.1f',Lambda(i)))
end
return




figure('Position',[45 5 1400 800]);title('ACF width');
subplot(2,2,1)
histogram(mainlobe_width(1,1,:),'BinWidth',HistoBinWidth);
hold on; histogram(mainlobe_width(2,1,:),'BinWidth',HistoBinWidth);
legend('Pos','Vel');
xlabel('Width [samples]'),ylabel('Count')
title(['acf Cp main-lobe (0.5) width for Pos vs. Vel Control, Lamba=',num2str(Lambda(1),2)])

subplot(2,2,2)
histogram(mainlobe_width(1,2,:),'BinWidth',HistoBinWidth);
hold on; histogram(mainlobe_width(2,2,:),'BinWidth',HistoBinWidth);
legend('Pos','Vel');
xlabel('Width [samples]'),ylabel('Count')
title(['acf Cp main-lobe (0.5) width for Pos vs. Vel Control, Lamba=',num2str(Lambda(2),2)])

subplot(2,2,3)
histogram(mainlobe_width(1,3,:),'BinWidth',HistoBinWidth);
hold on; histogram(mainlobe_width(2,3,:),'BinWidth',HistoBinWidth);
legend('Pos','Vel');
xlabel('Width [samples]'),ylabel('Count')
title(['acf Cp main-lobe (0.5) width for Pos vs. Vel Control, Lamba=',num2str(Lambda(3),2)])

subplot(2,2,4)
histogram(mainlobe_width(1,4,:),'BinWidth',HistoBinWidth);
hold on; histogram(mainlobe_width(2,4,:),'BinWidth',HistoBinWidth);
legend('Pos','Vel');
xlabel('Width [samples]'),ylabel('Count')
title(['acf Cp main-lobe (0.5) width for Pos vs. Vel Control, Lamba=',num2str(Lambda(4),2)])

return;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
