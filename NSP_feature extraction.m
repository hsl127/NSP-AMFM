function BCI_IV_1_CSP(subject_id,BAND,varargin)
% clear
if nargin < 1
    clear;clc;
    subject_id = 4;
end

tdnn_opts = [];
tdnn_opts.trainFcn  = @tdnn_train;
tdnn_opts.testFcn   = @tdnn_test;
tdnn_opts.v_range   = 0;
tdnn_opts.dn  = 5;
tdnn_opts.nn  = 4;
tdnn_opts.cpl_flg   = 4;

if nargin < 2
    BAND = [8 30];%CWT的频段,与滤波频段无关
elseif nargin > 2
    tdnn_opts.dn  = varargin{1};
    tdnn_opts.nn  = varargin{2};
end

subject_name_group = {'a','b','f','g'};
subject    = subject_name_group{subject_id};
fprintf('subject:%s(%02d-%02d)\n', subject,BAND(1),BAND(2));

train_erp_time      = [ 0.0, 4.0];
extract_train_time  = [-2.0, 8.0];




eegfilepath = sprintf('E:\\EEG_PRO\\HsirCSP_LDA\\BCI4_1_csp1\\BCICIV_calib_ds1%s', subject);
load(eegfilepath);
cnt= 0.1*double(cnt);%To convert it to uV values
fs = nfo.fs;


sel_channels = [11 15 29 27 31 26 32 25 33 36 39];%Channel selection, corresponding to {FC3 ', FC4', Cz ', C3', C4 ', C5', C6 ', T7', T8 ', CCP3', CCP4 '}

cnt = cnt(:,sel_channels);

time_start = extract_train_time(1);
time_end   = extract_train_time(2);
trials = length(mrk.pos);%not include rest state
base_line_start = mrk.pos+time_start*fs;
time = time_end - time_start;
Nosamples  = time*fs;
samplesidx = bsxfun(@plus,base_line_start(:),0:Nosamples-1)';
EEG = cnt(samplesidx(:),:); % samples - trials x channels
EEG = reshape(EEG',size(EEG,2),[],trials);
LABELS = mrk.y';
[LABELS,idx] = sort(LABELS);
EEG = EEG(:,:,idx);
% EEG = eeg_filt(EEG ,fs ,BAND);%EEG = eeglab_eegfiltnew(EEG ,fs ,BAND(1) ,BAND(2),[],0,[],0);
trim_start = fix( ( train_erp_time(1) - extract_train_time(1) ) * fs + 1 );
trim_end   = fix( ( train_erp_time(2) - extract_train_time(1) ) * fs);
EEG       = EEG(:,201:600,:);
NSP1206_g1=zeros(size(EEG,1),size(EEG,2),size(EEG,3));
NSP1206_g2=zeros(size(EEG,1),size(EEG,2),size(EEG,3))
NSP1206_g3=zeros(size(EEG,1),size(EEG,2),size(EEG,3))
NSP1206_g4=zeros(size(EEG,1),size(EEG,2),size(EEG,3))
for i=1:size(EEG,3)
    for j=1:size(EEG,1)
        sig=squeeze(EEG(j,:,i))
        res=AMFMNSP_4(sig',fs);
        NSP1206_g1(j,:,i)=res(1,:);
        NSP1206_g2(j,:,i)=res(2,:);
        NSP1206_g3(j,:,i)=res(3,:);
        NSP1206_g4(j,:,i)=res(4,:);
        
    end
end

save NSP1206_g1
save NSP1206_g2
save NSP1206_g3
save NSP1206_g4
end

















