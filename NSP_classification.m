clear ;
close all;
load BCICIV_g_fea.mat;
load BCICIV_calib_ds1g.mat;
LABELS = mrk.y';

EEG=BCICIV_g_fea;
svm_acc=zeros(1,100);
each_acc=zeros(1,length(0.001:0.002:5));
 for i=1:100
    idx_r=randperm(100);   %Generate random numbers
    A=idx_r(1:70);
    B=A+100;
    C=idx_r(71:100);
    D=C+100;
    trainidx=[A B]; %Take 70 samples from each class as the training set and record the sample IDs in the training set
    testidx = [C D];%Take 30 samples from each category as the test set and record the sample IDs in the test set

    gndtrain = LABELS(trainidx); gndtest = LABELS(testidx);
    gndtrain(gndtrain==-1) = 2; gndtest(gndtest==-1) = 2;
    [features,W] = CSPfeature(EEG,trainidx,gndtrain,5);
    featrain  = features(trainidx,:); %Training set features
    featest = features(testidx,:);  %Test set features
    ii=1;
    for parameter_set=0.0001:0.0005:5
        model=libsvmtrain(gndtrain ,featrain,'-c 2 -t 2 -g parameter_set');  
        [predictlabel,ac,decv]=libsvmpredict(gndtest,featest,model ,'libsvm_options' ); 
        each_acc(ii)=ac(1);
        ii=ii+1;
    end
    [max_acc,idex]=max(each_acc);
    svm_acc(i)=max_acc;
end
mean_svm=mean(svm_acc);
fprintf(['ACC = ',num2str(mean_svm)]);