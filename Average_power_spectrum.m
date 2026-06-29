clear all;close all
% Average power spectrum for the first 4 components

load NSP_a1.mat
load NSP_a2.mat
load NSP_a3.mat 
load NSP_a4.mat 



aa=NSP_a1
bb=NSP_a2
cc=NSP_a3
dd=NSP_a4
aa1(1,1,1:1800)=0;
aa2(1,1,1:1800)=0;
aa3(1,1,1:1800)=0;
aa4(1,1,1:1800)=0;
for i=1:size(aa,1)
    for ii=1:size(aa,2)
    whitenois1=aa(i,ii,1:1800)
    n=length(whitenois1);
    fft_white1=fft(whitenois1);
    fft_white_shift1=fftshift(fft_white1);
    realy1 =abs(fft_white_shift1).^2/n;
    realf1 = (-n/2:n/2-1)*(450/n);
    aa1=aa1+realy1
   
    end
end
aa1=aa1/2200;

for i=1:size(bb,1)
    for ii=1:size(bb,2)
    whitenois2=bb(i,ii,:)
    n=length(whitenois2);
    fft_white2=fft(whitenois2);
    fft_white_shift2=fftshift(fft_white2);
    realy2 =abs(fft_white_shift2).^2/n;
    realf2 = (-n/2:n/2-1)*(450/n);
    aa2=aa2+realy2
   
    end
end
aa2=aa2/2200;

for i=1:size(cc,1)
    for ii=1:size(cc,2)
    whitenois3=cc(i,ii,1:1800)
    n=length(whitenois3);
    fft_white3=fft(whitenois3);
    fft_white_shift3=fftshift(fft_white3);
    realy3 =abs(fft_white_shift3).^2/n;
    realf3 = (-n/2:n/2-1)*(450/n);
    aa3=aa3+realy3
   
    end
end
aa3=aa3/2200;


for i=1:size(dd,1)
    for ii=1:size(dd,2)
    whitenois4=dd(i,ii,1:1800)
    n=length(whitenois4);
    fft_white4=fft(whitenois4);
    fft_white_shift4=fftshift(fft_white4);
    realy4 =abs(fft_white_shift4).^2/n;
    realf4 = (-n/2:n/2-1)*(450/n);
    aa4=aa4+realy4;
   
    end
end
aa4=aa4/2200;




aa1=squeeze(aa1);
aa2=squeeze(aa2);
aa3=squeeze(aa3);
aa4=squeeze(aa4);
 figure
 plot(realf1,aa1,'r');
 xlim([0 50])
 hold on;
 plot(realf2,aa2,'b');
 xlim([0 50])
 hold on;
 plot(realf3,aa3,'k');
 xlim([0 50])
 hold on;
 plot(realf4,aa4,'g');
 xlim([0 50])
 hold on;

 legend('c_1(t)','c_2(t)','c_3(t)','C_4(t)');
 
 
 