function [U, R] = ObassAMFMSps( S, lambda, epsi, winsize, max_iter, freq_thresh, Ureal, Phireal, display )
%% Operator-based signal separation using AM-FM signal model
%% Written by Xiyuan, Hu, IACAS

%% initial local parameters
if ~exist( 'display', 'var' ), display = false; end
if ~exist('Phireal', 'var') || isempty( Phireal ), 
    Phireal = cell(1,1);
    Phireal{1} = zeros( size( S ) );
    realfreq = false;
else
    realfreq = true;
    freq_max = 0;
    freq_min = 100;
    for i = 1:length( Phireal ),
        freq_max = max( freq_max, max( Phireal{i}(:) ) );
        freq_min = min( freq_min, min( Phireal{i}(:) ) );
    end
    freq_max = freq_max + 0.1;
    freq_min = freq_min - 0.1;
end
if ~exist( 'Ureal', 'var' ) || isempty( Ureal ), 
    Ureal = zeros( size( S ) );
    realsig = false;
else
    realsig = true;
end
if ~exist('freq_thresh', 'var'), freq_thresh = 0; end
if ~exist('max_iter', 'var'), max_iter = 200; end
if ~exist('winsize', 'var'), winsize = 51; end
if ~exist('epsi', 'var'), epsi = 1e-3; end

filt1 = [-0.5 0 0.5]; 
filt2 = [1 -2 1];
gamma = 1;
psnr = 0;

%% initial debug variables
ls = length(S);
ErrPrev = ones(1,max_iter)*1e20;
R_curr = zeros(size(S));
R_prev = zeros(size(S));

% show the figures or not
if display,
    figure('Name','Iteration');
    noi=S-Ureal;
    axismax = max( 10*log10((Ureal*Ureal')/(noi*noi')),0 ) + 30;
    subplot(2,2,4), plot(1,2,'r+'), axis([1  max_iter  0  axismax]), hold on;
end

ind = 1;
while ind < max_iter,
    % estimate parameters IF and ILA
    [Q_curr, P_curr] = UpdateParameterPQAMFMSps( S, R_prev, filt1, filt2, gamma, winsize, 1e-7 );
    
    % whether display the result
    if display,
        subplot(2,2,1), plot( sqrt(abs(Q_curr-(P_curr.^2)/2)), 'r' ); axis([1 length(S) -0.5 0.5]);
        hold on;
        if realfreq,
            for ip = 1:length( Phireal ), 
                subplot(2,2,1), plot( Phireal{ip}, 'k' ); axis([1 length(S) -abs(freq_min) freq_max]);
            end
        end
        subplot(2,2,1), plot( P_curr, 'b' ); title('Estimated P and Q');
        hold off;
        drawnow;
    end

    % compute the new lambda
    [lambda] = UpdateParmsLamdaSps( S, filt1, filt2, P_curr, Q_curr, lambda, gamma, winsize );
    
    % compute the residual signal
    [R_curr] = UpdateResidualAMFMSps( S, filt1, filt2, P_curr, Q_curr, lambda, 1, freq_thresh );

    % compute the parameter gamma
    [gamma] = UpdateParmsGamma( S, R_curr, winsize );
    
    % compute the new extract component
    U = ( 1+gamma ) * ( S - R_curr );
    ErrCurr = dot(R_curr-R_prev, R_curr-R_prev)/ls;
    ErrPrev(ind) = ErrCurr;
    if realsig, 
        psnr = 10*log10( dot(Ureal,Ureal)/dot(Ureal-U,Ureal-U) );
    end
    
    % whether show the result
    if display,
        subplot(2,2,2), plot( U ), title('First Component');
        drawnow;
        if realsig,
            subplot(2,2,3), plot( U-Ureal ), title('Error Signal'), axis([1 length(S) -1 1]);
        else
            subplot(2,2,3), plot( S-U ), title('Residual Component');
        end
        drawnow;
        subplot(2,2,4), plot(ind, psnr, 'r+'); title(sprintf('PSNR: %.3f', psnr));
        drawnow;pause(.1);
    end

    % stop criteria 
    if ind > 5,
        smerr = sqrt( mean( ErrPrev(ind-5:ind-1) ) );
        if  smerr < epsi && abs(sqrt(ErrCurr)-smerr) < 0.05*smerr,
            break;
        elseif ind == max_iter,
            break;
        end
    end

    % update residual signal
    R_prev = R_curr;
    ind = ind + 1;
    
end

U = ( 1+gamma ) * ( S-R_curr );
R = S - U;

%% compute evaluation result
if realsig,
    noi=S-Ureal;
    psnr_prev = 10*log10((Ureal*Ureal')/(noi*noi'));
    noi=U-Ureal;
    psnr_curr = 10*log10((Ureal*Ureal')/(noi*noi'));
    fprintf( 'PSNR improved from %.4f to %.4f\n', psnr_prev, psnr_curr );
end

if display,
    figure;
    subplot(2,1,1), plot( abs( fftshift( fft( S ) ) ) );
    hold on;
    plot( abs( fftshift( fft( U ) ) ), 'r' );
    hold off;
    subplot(2,1,2), plot( abs( fftshift( fft( R ) ) ) );
end