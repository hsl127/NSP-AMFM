function [W, D, invW, nClass, labset] = CSP(EEG,LABELS,Nobases,multi_c_method)
% Seek spatial feature subspace for EEG ( channel x samples x trial)
% Phan Anh Huy, 2010

% EEG = double(EEG);
EEG = bsxfun(@rdivide,EEG,sqrt(sum(sum(EEG.^2,2)))); %normalize EEG's power
In = size(EEG);

% EEG1 = reshape(EEG(:,:,j == 1),In(1),[]);
% TF = isnan(EEG);
% EEG(TF(:,1),:) = [];
% LABELS1 = LABELS(TF(:,1)==0);

[labset,~,j] = unique(LABELS);

nClass = length(labset);%Number of classes
ClsIdxs= cell(nClass,1);
Ns    =  zeros(nClass,1);
for i=1:nClass
    ClsIdxs{i}=find(LABELS==labset(i));
    Ns(i)=length(ClsIdxs{i});
end

switch nClass
    case 2
% %         Ytrain = LABELS;
% %         Ytrain(Ytrain==2) = -1;
% %         S = diag(Ytrain);
% %         C1 = zeros(size(EEG,1));
% %         C2 = zeros(size(EEG,1));
% %         for i = 1:length(Ytrain)
% %             Y_i = EEG(:,:,i);
% %             for j = 1:length(Ytrain)
% %                 if (S(i,j) == 0)
% %                     continue;
% %                 end
% %                 Y_j = EEG(:,:,j);
% %                 if i==j
% %                    if Ytrain(i)==1
% %                       C1 = C1 + S(i,j)*Y_i*Y_j';%P_i*(P_j');
% %                    else
% %                       C2 = C2 + S(i,j)*Y_i*Y_j';%*P_i*(P_j'); 
% %                    end
% %                 else
% %                    if Ytrain(i)==labels(1)
% %                       S1 = S1 + S(i,j)*Y_i*(Y_j');%P_i*(P_j');
% %                    else
% %                       S2 = S2 + S(i,j)*Y_i*(Y_j');%P_i*(P_j'); 
% %                    end
% % %                    S2 = S2 + S(i,j)*P_i*(P_j'); 
% %                 end
% %             end  
% %         end
% %          C1 = C1/Ns(1);C2 = -C2/Ns(2);
        
        
%         COVtrain = covariances_h(EEG);
%         C1 = mean(COVtrain(:,:,j == labset(1)),3);%C1 = max(C1,C1');
%         C2 = mean(COVtrain(:,:,j == labset(2)),3);%C2 = max(C2,C2');
        
        EEG1 = reshape(EEG(:,:,j == labset(1)),In(1),[]);%EEG1 = reshape(EEG(:,:,j == 1),In(1),[]);     %
        C1 = EEG1*EEG1'/Ns(1);


        C1 = max(C1,C1');         %C1 = EEG1*EEG1'/sum(j == 1); C1 = max(C1,C1');%

        EEG2 = reshape(EEG(:,:,j == labset(2)),In(1),[]);%EEG2 = reshape(EEG(:,:,j == 2),In(1),[]);
        C2 = EEG2*EEG2'/Ns(2); C2 = max(C2,C2');         %C2 = EEG2*EEG2'/sum(j == 1); C2 = max(C2,C2');
        Ct = C1+C2;
        
        % Find bases W
        [W1,D] = eig(C1,Ct); 
        [D,index]=sort(diag(D),'descend');
        
        %after 20201216
        W2 = W1(:,index);
        W2 = bsxfun(@rdivide,W2,sqrt(sum(W2.^2)));
        invW = inv(W2');
        invW = invW(:,[1:Nobases/2 end-Nobases/2+1:1:end]);
        W2 = W2(:,[1:Nobases/2 end-Nobases/2+1:end]);
        W  = W2';
%         W1 = W1(:,index([end:-1:end-Nobases/2+1 1:Nobases/2]));%20200924
% % %         W1 = W1(:,index([1:Nobases/2 end-Nobases/2+1:1:end]));%before 20201216
% % %         W1 = bsxfun(@rdivide,W1,sqrt(sum(W1.^2)));
% % %         W  = W1';
        
%         obj = zeros(Nobases,1);
% %         fprintf('obj:');
%         for ki=1:Nobases
%           obj(ki)  = (W(ki,:)*C1*W(ki,:)')/(W(ki,:)*Ct*W(ki,:)');
%           fprintf('%.3f  ',obj(ki));
%         end
%         fprintf('\n');
     case {3, 4}
        if isempty(multi_c_method)
            W = cell(nClass,1);
            for c=1:nClass
                EEG1 = reshape(EEG(:,:,j == c),In(1),[]);%chn x times
                EEG2 = reshape(EEG(:,:,j ~= c),In(1),[]);%chn x times
                ivec_n    = 1:nClass;
                ivec_n(c) = [];
                C1 = EEG1*EEG1'/Ns(c); C1 = max(C1,C1');
                C2 = EEG2*EEG2'/sum(Ns(ivec_n)); C2 = max(C2,C2');

                [W0,D] = eig(C1, C1+C2); 
                [D,index]=sort(diag(D));
                W0 = W0(:,index([end:-1:end-Nobases/2+1 1:Nobases/2]));
                W0 = bsxfun(@rdivide,W0,sqrt(sum(W0.^2)));
                W{c,1} = W0';
            end
        else %仅适用用于BCI IV_1 其它应用需要调整
            switch multi_c_method
                 case 0%DC
                    W = cell(nClass-1,1);

                    EEG1 = reshape(EEG(:,:,j == 2),In(1),[]);%chn x times
                    EEG2 = reshape(EEG(:,:,j ~= 2),In(1),[]);%chn x times
                    ivec_n    = 1:nClass;
                    ivec_n(2) = [];
                    C1 = EEG1*EEG1'/Ns(2); C1 = max(C1,C1');
                    C2 = EEG2*EEG2'/sum(Ns(ivec_n)); C2 = max(C2,C2');
                    [W0,D] = eig(C1, C1+C2); 
                    [D,index]=sort(diag(D));
                    W0 = W0(:,index([end:-1:end-Nobases/2+1 1:Nobases/2]));
                    W0 = bsxfun(@rdivide,W0,sqrt(sum(W0.^2)));
                    W{1,1} = W0';

                    EEG1 = reshape(EEG(:,:,j == 1),In(1),[]);%chn x times
                    EEG2 = reshape(EEG(:,:,j == 3),In(1),[]);%chn x times
                    C1 = EEG1*EEG1'/Ns(1); C1 = max(C1,C1');
                    C2 = EEG2*EEG2'/Ns(3); C2 = max(C2,C2');
                    [W0,D] = eig(C1, C1+C2); 
                    [D,index]=sort(diag(D));
                    W0 = W0(:,index([end:-1:end-Nobases/2+1 1:Nobases/2]));
                    W0 = bsxfun(@rdivide,W0,sqrt(sum(W0.^2)));
                    W{2,1} = W0';
                case 1%PW
                    pw_grp = [2 1;2 3; 1 3];
                    W = cell(size(pw_grp,1),1);
                    for c=1:size(pw_grp,1)
                        EEG1 = reshape(EEG(:,:,j == pw_grp(c,1)),In(1),[]);%chn x times
                        EEG2 = reshape(EEG(:,:,j == pw_grp(c,2)),In(1),[]);%chn x times
                        C1 = EEG1*EEG1'/Ns(pw_grp(c,1)); C1 = max(C1,C1');
                        C2 = EEG2*EEG2'/Ns(pw_grp(c,2)); C2 = max(C2,C2');
                        [W0,D] = eig(C1, C1+C2); 
                        [D,index]=sort(diag(D));
                        W0 = W0(:,index([end:-1:end-Nobases/2+1 1:Nobases/2]));
                        W0 = bsxfun(@rdivide,W0,sqrt(sum(W0.^2)));
                        W{c,1} = W0';
                    end
                    
                    
            end
        end
            
        
        
        
    otherwise
        return
end
% if nClass==2 %2分类
%     
% else
%    for i=1:nClass
%        EEG1 = reshape(EEG(:,:,j == labset(i)),In(1),[]);%chn x times
%        cov  = covm(EEG1','E');%/Ns(i);
%        if i==1
%           n = size(cov,1); 
%           ECM = zeros(nClass,n,n); 
%        end
%        ECM(i,:,:) = cov;
%    end
%    	pat = 1;
% 	sz = size(ECM); 
% 
% 	for k=1:sz(1)
% 		[mu,sd,COV(k,:,:),xc,N,R2]=decovm(squeeze(ECM(k,:,:)));
%     end 
% 
% 
% 	V = repmat(NaN,sz(2)-1,2*pat*sz(1));
% 	d = V(1,:);
%     R = permute(COV,[2,3,1]);
%     for k = 1:sz(1) 
%         [W,D] = eig(R(:,:,k),sum(R,3));
%         V(:,2*pat*k+[1-2*pat:0]) = W(:,[1:pat,end-pat+1:end]);
%         d(1,2*pat*k+[1-2*pat:0]) = diag(D([1:pat,end-pat+1:end],[1:pat,end-pat+1:end]));
%     end
%     W = V';
% end




return
end


function COV = covariances_h(X)


[Ne , ~, Nt] = size(X);

COV = zeros(Ne,Ne,Nt);
for i=1:Nt
    COV(:,:,i) = X(:,:,i)*X(:,:,i)';%/trace(X(:,:,i)*X(:,:,i)');
end

end