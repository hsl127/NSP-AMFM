function res=AMFMNSP_4(s,fs)
s=squeeze(s)';
y =s;

t2=0:1/fs:4-1/fs;

res = zeros(4,length(y));

% parameters definition
lambda = 0.01;         % lagrange parameters
epsi = 0.0001;           % stop criteria
winsize = 51;       % window size of estimate parameters P and Q

fprintf( 'extracting first subcomponent...\n' );
[res(1,:), r] = ObassAMFMSps( y, lambda, epsi, winsize );

lambda = 0.01;
epsi = 0.0001;
winsize = 31;
fprintf( 'extracting second subcomponent...\n' );
[res(2,:), r] = ObassAMFMSps( r, lambda, epsi, winsize );

fprintf( 'extracting third subcomponent...\n' );
[res(3,:), r] = ObassAMFMSps( r, lambda, epsi, winsize );

fprintf( 'extracting fourth subcomponent...\n' );
[res(4,:), r] = ObassAMFMSps( r, lambda, epsi, winsize );

end