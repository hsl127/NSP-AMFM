function [features,W,D, invW,nClass, labset] = CSPfeature(X,trainidx,gndtrain,Nofeat,varargin)
% Copyright of this implementation by Anh Huy Phan, 2010.

if nargin==5
   all_opts = varargin{1};
   multi_c_method = all_opts.mc;
   if multi_c_method==0
      multi_c_method = []; 
   end
   complex_flg = all_opts.complex_flg;
else
   multi_c_method = [];
   complex_flg = 0;
end
if isempty(trainidx)
    [W, D,invW,nClass, labset] =CSP(X,gndtrain,2*Nofeat,multi_c_method);
else
    [W, D,invW,nClass, labset] =CSP(X(:,:,trainidx),gndtrain,2*Nofeat,multi_c_method);
%     [W, nClass, labset] =CSP2(X(:,:,trainidx),gndtrain,2*Nofeat);
%     [W, nClass, labset] =STRCSP(X(:,:,trainidx),gndtrain,2*Nofeat,1,0);
end


if nClass>2
   if multi_c_method~=0%isempty(multi_c_method)
      kn = nClass;   %ovr or pw
   else
      kn = nClass-1; %dc
   end
   features = cell(kn,1);
   for i=1:kn
       features{i,1} = feaCSP(X,W{i,1});
   end
else
%2 class
   if complex_flg==0
      features = feaCSP(X,W);
   else
      features = feaComplexCSP(X,W);  
   end

end
