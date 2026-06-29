function q = memd(x, varargin)
%
%
% function MEMD applies the "Multivariate Empirical Mode Decomposition" algorithm (Rehman and Mandic, Proc. Roy. Soc A, 2010)
% to multivariate inputs. We have verified this code by simulations for signals containing 3-16 channels.
%函数MEMD应用“多元经验模态分解”到多变量的输入。我们对包含3-16信道的信号进行了验证仿真
% Syntax:
%
% imf = MEMD(X)
%   returns a 3D matrix 'imf(N,M,L)' containing M multivariate IMFs, one IMF per column, computed by applying
%   the multivariate EMD algorithm on the N-variate signal (time-series) X of length L.
% 返回一个3D矩阵'imf(N,M,L)'，包含M个多元imf，每列一个imf，通过应用计算
% 长度为L的n变量信号(时间序列)X上的多变量EMD算法。
%    - For instance, imf_k = IMF(k,:,:) returns the k-th component (1 <= k <= N) for all of the N-variate IMFs.
%imf_k = IMF(k，:，:)返回第K个分量的所有IMF(1 <= k <= N)
%   For example,  for hexavariate inputs (N=6), we obtain a 3D matrix IMF(6, M, L)
%   where M is the number of IMFs extracted, and L is the data length.
%对于6变量的输入，M代表提取IMF的数量，L是向量的长度

% imf = MEMD(X,num_directions)
%   where integer variable num_directions (>= 1) specifies the total number of projections of the signal
%     - As a rule of thumb, the minimum value of num_directions should be twice the number of data channels,
%     - for instance, num_directions = 6  for a 3-variate signal and num_directions= 16 for an 8-variate signal
%   The default number of directions is chosen to be 64 - to extract meaningful IMFs, the number of directions
%   should be considerably greater than the dimensionality of the signals
%
% imf = MEMD(X,num_directions,'stopping criteria')
%   uses the optional parameter 'stopping criteria' to control the sifting process.
%    The available options are
%      -  'stop' which uses the standard stopping criterion specified in [2]
%      -  'fix_h' which uses the modified version of the stopping criteria specified in [3]
%    The default value for the 'stopping criteria' is 'stop'.
%
%  The settings  num_directions=64 and 'stopping criteria' = 'stop' are defaults.
%     Thus imf = MEMD(X) = MEMD(X,64) = MEMD(X,64,'stop') = MEMD(X,[],'stop'),
%
% imf = MEMD(X, num_directions, 'stop', stop_vec)
%   computes the IMFs based on the standard stopping criterion whose parameters are given in the 'stop_vec'%根据标准停止准则计算imf，其参数在'stop_vec'中给出
%     - stop_vec has three elements specifying the threshold and tolerance values used, see [2].
%     - the default value for the stopping vector is   step_vec = [0.075 0.75 0.075].
%     - the option 'stop_vec' is only valid if the parameter 'stopping criteria' is set to 'stop'.
%
% imf = MEMD(X, num_directions, 'fix_h', n_iter)
%   computes the IMFs with n_iter (integer variable) specifying the number of consecutive iterations when
%   the number of extrema and the number of zero crossings differ at most by one [3].
%     - the default value for the parameter n_iter is set to  n_iter = 2.
%     - the option n_iter is only valid if the parameter  'stopping criteria' = 'fix_h'
%
%
% This code allows to process multivaraite signals having 3-16 channels, using the multivariate EMD algorithm [1].
%   - to perform EMD on more than 16 channels, modify the variable 'Max_channels' on line 510 in the code accordingly.
%   - to process 1- and 2-dimensional (univariate and bivariate) data using EMD, we recommend the toolbox from
%                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%
% Acknowledgment: Part of this code is based on the bivariate EMD code, publicly available from
%                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html. We would also like to thank 
%                 Anh Huy Phan from RIKEN for helping us in optimizing the code and making it computationally efficient. 
%
%
% Copyright: Naveed ur Rehman and Danilo P. Mandic, Oct-2009
%
%
% [1]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010
% [2]  G. Rilling, P. Flandrin and P. Goncalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
%      Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003
% [3]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
%      Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Usage %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 1:

% inp = randn(1000,3);
% imf = memd(inp);
% imf_x = reshape(imf(1,:,:),size(imf,2),length(inp)); % imfs corresponding to 1st component
% imf_y = reshape(imf(2,:,:),size(imf,2),length(inp)); % imfs corresponding to 2nd component
% imf_z = reshape(imf(3,:,:),size(imf,2),length(inp)); % imfs corresponding to 3rd component


% Case 2:

% load syn_hex_inp.mat
% imf = memd(s6,256,'stop',[0.05 0.5 0.05])

global N N_dim;  %设置全局变量
[x, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stop_crit, stp_cnt] = set_value(x, nargin, varargin{:});

r=x;  %把X赋值给r
n_imf=1;  %这个n_imf等于1
q = zeros(N_dim,1,N);% q此刻变成了3维
while ~stop_emd(r, seq, ndir)
    % current mode  这是当前模式，因为第一次赋值是r等于x
    m = r;  
     
    % computation of mean and stopping criterion  均值和停止计算的标准
    if(strcmp(stop_crit,'stop'))   %tf = strcmp(s1,s2) 比较 s1 和 s2，如果二者相同，则返回 1 (true)，否则返回 0 (false)。
        [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
    else
        counter=0;
        [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
    end
    
    % In case the current mode is so small that machine precision can cause
    % spurious extrema to appear
    if (max(abs(m))) < (1e-10)*(max(abs(x)))
        if ~stop_sift
            warning('emd:warning','forced stop of EMD : too small amplitude')
        else
            disp('forced stop of EMD : too small amplitude')
        end
        break
    end
    
    % sifting loop
    while ~stop_sift && nbit<MAXITERATIONS
        %sifting
        m = m - env_mean;
        % computation of mean and stopping criterion
        if(strcmp(stop_crit,'stop'))
            [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
        else
            [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
        end
    
        nbit=nbit+1;
        
        if(nbit==(MAXITERATIONS-1) &&  nbit > 100)
            warning('emd:warning','forced stop of sifting : too many iterations');
        end
    end
    
    q(:,n_imf,:)=m';   %这里就是计算imf的地方，最后的返回值也是返回q
    
    n_imf = n_imf+1;
    r = r - m;  %计算残差
    nbit = 0;
end
% Stores the residue
q(:,n_imf,:)=r';  %最后的残差分量在这

%sprintf('Elapsed time: %f\n',toc);
end

%---------------------------------------------------------------------------------------------------
function stp = stop_emd(r, seq, ndir)
global N_dim;
ner = zeros(ndir,1);
for it=1:ndir%64个方向就会有64个角度，然后就是在这计算停止的投影
    if (N_dim~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        %hammersley序列的多元信号(对于N_dim ~=3)在-1.00 -1.00范围内对hammersley序列进行线性归一化


        b=2*seq(1:end,it)-1;
        
        % Find angles corresponding to the normalised sequence找到与归一化序列对应的角度
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];
        dir_vec(1:N_dim-1) =  cos(tht) .*dir_vec(1:N_dim-1);

    else % Trivariate signal with hammersley sequence含hammersley序列的三元信号
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1))=1;
        tt((tt<-1))=-1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    % Projection of input signal on nth (out of total ndir) direction
    % vectors
    y = r * dir_vec';   %输入信号在第n个(非全ndir)方向上的投影
    % Calculates the extrema of the projected signal
    [indmin, indmax] = local_peaks(y);
    
    ner(it) = length(indmin) + length(indmax);
end

% Stops if the all projected signals have less than 3 extrema
stp = all(ner < 3)%在这里仅仅是计算64个方向上的投影，每个方向上的曲线极值点的情况以及相关的反馈
end

%---------------------------------------------------------------------------------------------------
% computes the mean of the envelopes and the mode amplitude estimate 计算包络线的平均值和模幅估计
function [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir) %new
global N N_dim;
NBSYM = 2;
count=0;

env_mean=zeros(length(t),N_dim);   %先定义一个数组
amp = zeros(length(t),1);      %用来装赋值的叭？
nem = zeros(ndir,1);
nzm = zeros(ndir,1);
for it=1:ndir
    if (N_dim ~=3) % Multivariate signal (for N_dim ~=3) with hammersley sequence  多变量信号(对于N_dim ~=3)与hammersley序列
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00  hammersley序列在-1.00 -1.00范围内的线性归一化
        b=2*seq(1:end,it)-1;  %hammersley序列在-1.00 -1.00范围内的线性归一化
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere  求n球上单位方向向量的坐标
        dir_vec(1:N_dim) = [1 cumprod(sin(tht))];  %累积乘积  n维球面
        dir_vec(1:N_dim-1) =  cos(tht) .*dir_vec(1:N_dim-1);  %n-1维的方向向量
    else % Trivariate signal with hammersley sequence  含hammersley序列的三元信号
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0  hammersley序列在-1.0 ~ 1.0范围内的线性归一化
        tt = 2*seq(1,it)-1;
        tt((tt>1))=1;
        tt((tt<-1))=-1;
        
        % Normalize angle from 0 - 2*pi   从0 - 2归一化角度
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    
    % Projection of input signal on nth (out of total ndir) direction vectors  输入信号在第n个(总ndir)方向向量上的投影
    y(1:N)  = dir_vec * m(1:N,:)';  %第一次是原始信号在这个方向向量上的投影，每个点都会对应一个投影
    %为什么会有64个方向是因为在前几行代码中每个方向向量都定义了一个角度
    %每一次it循环都会有不同的相切的角度，相当于把原始信号投影到不同的维度
    
    % Calculates the extrema of the projected signal  计算投影信号的极值
    [indmin, indmax] = local_peaks(y);  %计算峰值
    
    
    nem(it) = length(indmin) + length(indmax);
    
    indzer = zero_crossings(y);
    nzm(it) = length(indzer);
    
    [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,y,m,NBSYM);
    
    % Calculate multidimensional envelopes using spline interpolation
    % Only done if number of extrema of the projected signal exceed 3 
    %使用样条插值计算多维包络
    %仅当投影信号的极值数超过3时执行
    if(mode)
        env_min = spline(tmin,zmin.',t).';
        env_max = spline(tmax,zmax.',t).';
        
        amp = amp + sqrt(sum((env_max-env_min).^2,2))/2;
        env_mean = env_mean + (env_max+env_min)/2;%整体信号在64个方向各自的平均相加就是均值
    else % if the projected signal has inadequate extrema
        count=count+1;
    end
end
if(ndir>count)
    env_mean = env_mean/(ndir-count);
    amp = amp/(ndir-count);
else
    env_mean = zeros(N,N_dim);
    amp = zeros(N,1);
    nem = zeros(1,ndir);
end
end

%-------------------------------------------------------------------------------
% Stopping criterion
function [stp,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir)
global N N_dim;
try
    [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir);
    sx = sqrt(sum(env_mean.^2,2));
    if(amp) % something is wrong here
        sx = sx./amp;
    end
    stp = ~((mean(sx > sd) > tol | any(sx > sd2)) & any(nem > 2));
catch
    env_mean = zeros(N,N_dim);
    stp = 1;
end
end

function [stp,env_mean,counter]= stop_sifting_fix(m,t,seq,ndir,stp_count,counter)
global N N_dim;%定义全局变量
try
    [env_mean,nem,nzm] = envelope_mean(m,t,seq,ndir);
    if (all(abs(nzm-nem)>1))
        stp = 0;
        counter = 0;
    else
        counter = counter+1;
        stp = (counter >= stp_count);
    end
catch
    env_mean = zeros(N,N_dim);
    stp = 1;
end
end

%---------------------------------------------------------------------------------------
% defines new extrema points to extend the interpolations at the edges of the
% signal (mainly mirror symmetry)定义新的极值点来扩展信号边缘的插值(主要是镜像对称)
function [tmin,tmax,zmin,zmax,mode] = boundary_conditions(indmin,indmax,t,x,z,nbsym)
%这是求上下包络的
lx = length(x);
if (length(indmin) + length(indmax) < 3)
    mode = 0;
    tmin=NaN;tmax=NaN;zmin=NaN;zmax=NaN;
    return
else
    mode=1; %the projected signal has inadequate extrema 投影信号极值不足
end
% boundary conditions for interpolations :插值边界条件:
if indmax(1) < indmin(1)
    if x(1) > x(indmin(1))
        lmax = fliplr(indmax(2:min(end,nbsym+1)));
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = indmax(1);
    else
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
        lsym = 1;
    end
else
    
    if x(1) < x(indmax(1))
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = fliplr(indmin(2:min(end,nbsym+1)));
        lsym = indmin(1);
    else
        lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = 1;
    end
end

if indmax(end) < indmin(end)
    if x(end) < x(indmax(end))
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = fliplr(indmin(max(end-nbsym,1):end-1));
        rsym = indmin(end);
    else
        rmax = [lx,fliplr(indmax(max(end-nbsym+2,1):end))];
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = lx;
    end
else
    if x(end) > x(indmin(end))
        rmax = fliplr(indmax(max(end-nbsym,1):end-1));
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = indmax(end);
    else
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-nbsym+2,1):end))];
        rsym = lx;
    end
end
tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough 如果对称的部分没有足够的延伸
if tlmin(1) > t(1) || tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,nbsym)));
    else
        lmin = fliplr(indmin(1:min(end,nbsym)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(lx) || trmax(end) < t(lx)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
    else
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
    end
    if rsym == lx
        error('bug')
    end
    rsym = lx;
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
end
zlmax =z(lmax,:);
zlmin =z(lmin,:);
zrmax =z(rmax,:);
zrmin =z(rmin,:);

tmin = [tlmin t(indmin) trmin];
tmax = [tlmax t(indmax) trmax];
zmin = [zlmin; z(indmin,:); zrmin];
zmax = [zlmax; z(indmax,:); zrmax];
end

function [indmin, indmax] = local_peaks(x)
if(all(x < 1e-5))
    x=zeros(1,length(x));
end
m = length(x);
% Calculates the extrema of the projected signal 计算投影信号的极值
% Difference between subsequent elements: 后续元素的区别:
dy = diff(x); a = find(dy~=0);
lm = find(diff(a)~=1) + 1;
d = a(lm) - a(lm-1);
a(lm) = a(lm) - floor(d/2);
a(end+1) = m;
ya  = x(a);

if(length(ya) > 1)
    % Maxima
    [pks_max,loc_max]=peaks(ya);
    % Minima
    [pks_min,loc_min]=peaks(-ya);
    
    if(~isempty(pks_min))
        indmin = a(loc_min);
    else
        indmin = NaN;
    end
    if(~isempty(pks_max))
        indmax = a(loc_max);
    else
        indmax = NaN;
    end
else
    indmin=NaN;
    indmax=NaN;
end
end

function [pks_max,locs_max] =peaks(X)
dX = sign(diff(X));
locs_max = find((dX(1:end-1) >0) &  (dX(2:end) <0)) + 1;
pks_max = X(locs_max);
end


function indzer = zero_crossings(x)
indzer = find(x(1:end-1).*x(2:end)<0);
if any(x == 0)
    iz = find( x==0 );
    if any(diff(iz)==1)
        zer = x == 0;
        dz = diff([0 zer 0]);
        debz = find(dz == 1);
        finz = find(dz == -1)-1;
        indz = round((debz+finz)/2);
    else
        indz = iz;
    end
    indzer = sort([indzer indz]);
end
end

function seq = hamm(n,base)
seq = zeros(1,n);  %定义一个数组
if ( 1 < base )
    seed = 1:1:n;
    base_inv = inv(base);  %矩阵求逆
    while ( any ( seed ~= 0 ) )
        digit = mod (seed(1:n), base);
        seq = seq + digit * base_inv;
        base_inv = base_inv / base;
        seed = floor (seed / base );
    end
else
    temp = 1:1:n;
    seq = (mod(temp,(-base + 1 ))+0.5)/(-base);  %为什么要这样
end
end

function [q, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stp_crit, stp_cnt] = set_value(q, narg, varargin)  %设置初始值的函数

error(nargchk(1,4,narg));  %nargchk函数用来显示输入的参数是否足够，error用来显示nargchk报了什么错误
ndir = [];   %赋初始值
stp_crit = [];%赋初始值
stp_vec = [];%赋初始值
stp_cnt  = [];%赋初始值
MAXITERATIONS  = [];%赋初始值
sd=[];%赋初始值
sd2=[];%赋初始值
tol=[];%赋初始值
prm= [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149];%这一长串是在干嘛？？？

% Changes the input vector to double vector   将输入向量改为二重向量
q = double(q);  %把输入改为双精度型

% Specifies maximum number of channels that can be processed by the code
% 指定代码可处理的最大通道数，这里是指定了16个通道
% Its maximum possible value is 32.
Max_channels = 16;

if(narg==2)
    ndir=varargin{1};  %64指默认的64个方向  这是指N维空间里的好几个方向向量
end

if(narg==3)
    if(~isempty(varargin{1}))
        ndir=varargin{1};
    else
        ndir=64;
    end
    stp_crit=varargin{2};
end

if(narg==4 && strcmp(varargin{2},'fix_h'))
    if(isempty(varargin{1}))
        ndir=64;
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    end
elseif (narg==4 && strcmp(varargin{2},'stop'))
    if(isempty(varargin{1}))
        ndir=64;
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    end
elseif (narg==4 && ~xor(strcmp(varargin{2},'fix_h'),strcmp(varargin{2},'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

%%%%%%%%%%%%%% Rescale input signal if required   如果需要的话，改变输入信号的尺度
if (any(size(q)) == 0)   %确定任何数组元素是否为非零
    datamsgid = generatemsgid('emptyDataSet');
    error(datamsgid,'Data set cannot be empty.');  %检查有没有出错
end
if size(q,1) < size(q,2)  %需要输入的维度其实是采样点*通道数
    q=q';
end

%%%%%%%%%%%% Dimension of input signal
N_dim = size(q,2);  %确定输入的维度
if(N_dim < 3 || N_dim > Max_channels)  %如果通道数小于3或者通道数大于刚才定义的最大通道数16，就会报错
    error('Function only processes the signal having 3 and 16 channels.');
end

%%%%%%%%%%%% Length of input signal  输入信号的长度
N = size(q,1);   

%%%%%%%%%%%%% Check validity of Input parameters  检查输入信号的有效性
if ~isempty(ndir) && (~isnumeric(ndir) || ~isscalar(ndir) || any(rem(ndir,1)) || (ndir < 6))
    Nmsgid = generatemsgid('invalid num_dir');
    error(Nmsgid,'num_dir should be an integer greater than or equal to 6.');%报出方向必须大于等于6的整数
end

if ~isempty(stp_crit) && (~ischar(stp_crit) || ~xor(strcmp(stp_crit,'fix_h'),strcmp(stp_crit,'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');   %停止标准应该是固定h或停止
end

if ~isempty(stp_vec) && (~isnumeric(stp_vec) || length(stp_vec)~=3 || ~strcmp(stp_crit,'stop'))
    Nmsgid = generatemsgid('invalid stop_vector');
    error(Nmsgid,'stop_vector should be an array with three elements e.g. default is [0.075 0.75 0.075] ');  %这个应该输入3个数
end

if ~isempty(stp_cnt) && (~isnumeric(stp_cnt) || ~isscalar(stp_cnt) || any(rem(stp_cnt,1)) || (stp_cnt < 0) || ~strcmp(stp_crit,'fix_h'))
    Nmsgid = generatemsgid('invalid stop_count');
    error(Nmsgid,'stop_count should be a nonnegative integer');%停止计数应该是非负整数
end

if (isempty(ndir))
    ndir=64; % default  %默认值是输入64个方向
end

if (isempty(stp_crit))
    stp_crit='stop'; % default
end

if (isempty(stp_vec))
    stp_vec=[0.075,0.75,0.075]; % default
end

if (isempty(stp_cnt))
    stp_cnt=2; % default
end

if(strcmp(stp_crit,'stop'))
    sd = stp_vec(1);
    sd2 = stp_vec(2);
    tol = stp_vec(3);
end

%%%%%%%%%%%%% Initializations for Hammersley function       Hammersley函数的初始化
base(1) = -ndir;  %base(1)等于负的64是几个意思

%%%%%%%%%%%%%% Find the pointset for the given input signal  %找到输入信号的点集
if(N_dim==3)
    base(2) = 2;
    for it=1:N_dim-1
        seq(it,:) = hamm(ndir,base(it));
    end
  
else
    for iter = 2 : N_dim
        base(iter) = prm(iter-1);
    end
    
    for it=1:N_dim
        seq(it,:) = hamm(ndir,base(it));
    end
end

%%%%%%%%%%%% Define t
t=1:N;

% Counter
nbit=0;
MAXITERATIONS=1000; % default

% tic
end