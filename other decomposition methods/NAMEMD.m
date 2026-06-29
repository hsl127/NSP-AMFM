function  IMF = NAMEMD(no_noi,sig,fs,wmM,support,numdir,factor)

% NAMEMD calculated the intrinsic mode functions (IMF) of a signal (sig)
% using noise-assisted MEMD.   使用噪声辅助MEMD
% Inputs:
% no_noi:   No. noise channels (usually 1-12) 噪声通道一般为1-12道
% sig:      Signal to analyze (channels x samples) 要分析的信号
% fs:       Sampling frequency   采样率
% support:  Whether to use channel 2:end as support for the MEMD only (0)
% or to use these channels for IMF extractions as well (1)
% numdir:   No. directions to use in the MEMD (default 64)   在MEMD中使用的方向(默认64)
% factor:   Factor difference between signal and noise channel power (default 1)
%信号与噪声信道功率之差因数(默认为1)
% Outputs:
% IMFavg:   Mean of IMFs across realizations for each band  各波段实现中imf的平均值
% IMFmed:   Median of IMFs across realizations for each band   各波段实现中imf的中位数
% nans:     No. times an IMF in any band is missing  
% IMFreps:  All maximum powered IMFs for all realization and bands
% 所有频带中能量最大的IMF
%
% Created by Sofie Therese Hansen 2016, edited jan. 2019.
% Ref: "Unmixing oscillatory brain activity by EEG source localization and
% empirical mode decomposition", by ST Hansen et al.
%
% Function MEMD applies the "Multivariate Empirical Mode Decomposition" algorithm (Rehman and Mandic, Proc. Roy. Soc A, 2010)
% to multivariate inputs. Code verified by the authors in simulations for signals containing 3-16 channels.


if nargin<7
    factor=1;
end

if nargin<6
    numdir=64;
end
if nargin<5
    support=0; % Use extra channels for support only
end
if nargin<4
    wmM=1;
end

[no_chan,samps]=size(sig);sampsO=samps;
if support==1
    no_chan=1;
end
st=75;samps=samps-2*st; % remove edge effects
calc_max=1; % find the imf with maximum power within band
reps=30; % No. realizations to average over
IMFreps=NaN(3,reps,sampsO);
% frequency bands of interest:
alpha = [7 14];
sbeta = [14,22];
fbeta = [22,30];
nans=0;
% for rep=1:reps
%     rng(rep)
    neutral_pow = mean(sig(:).^2); % estimate average signal power
    noise_chan=randn(no_noi,size(sig,2))*0.25;
   
    IMF=memd([sig;noise_chan],numdir); % permform NA-(M)EMD
    
%     no_imfs=size(IMF,2);
%     instAmp = NaN(no_chan,no_imfs,sampsO);instFreq = NaN(no_chan,no_imfs,sampsO);PHI = NaN(no_chan,no_imfs,sampsO);
%     for chan=1:no_chan
%         [instAmp(chan,:,:),instFreq(chan,:,:),PHI(chan,:,:)] = INST_FREQ_local(squeeze(IMF(chan,:,:)));% calculating the instantaneous amplitude and the instantaneous frequency
%     end
%     
%     instAmp=instAmp(:,:,st+1:end-st);instFreq=instFreq(:,:,st+1:end-st); % reduce edge effects in frequency estimation
%     if wmM==1
%         wm = mean(instFreq.*fs,3);
%     else
%         wm = sum(instFreq.*(instAmp./repmat(sum(instAmp,3),[1,1,samps]))*fs,3);
%     end
%     
%     ampSum=sum((instAmp(:,:,:)).^2,3);
%     if calc_max==1
%         for fr=1:3
%             switch fr
%                 case 1
%                     freq=alpha;
%                 case 2
%                     freq=sbeta;
%                 case 3
%                     freq=fbeta;
%             end
%             al= find(wm>=freq(1)& wm<freq(2)); % find imf in band
%             if isempty(al)
%                 nans=1+nans;
%             else
%                 [~,idx]=max(ampSum(al)); % find imf with max power in band
%                 [ch, idx]=ind2sub(size(wm),al(idx));
%                 imf=squeeze(IMF(ch,idx,:));
%                 IMFreps(fr,rep,:)=imf;
%             end
%         end
%     end
% end
% IMFavg=squeeze(nanmean(IMFreps,2));
% IMFmed=squeeze(nanmedian(IMFreps,2));
% 
