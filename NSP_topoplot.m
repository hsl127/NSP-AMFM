
clc;
clear;
close all;
load('BCICIV_g_fea.mat');

imf=BCICIV_g_fea;%Processed data, the first 100 samples are left-handed motor imagery, and the last 100 are right-handed motor imagery

channel =        { 'AF3','AF4','F5','F3','F1',...
                                 'Fz','F2','F4','F6','FC5',...
                                 'FC3','FC1','FCz','FC2','FC4',...
                                 'FC6','CFC7','CFC5','CFC3','CFC1',...
                                 'CFC2','CFC4','CFC6','CFC8','T7',...
                                 'C5','C3','C1','Cz','C2',...
                                 'C4','C6','T8','CCP7','CCP5',...
                                 'CCP3','CCP1','CCP2','CCP4','CCP6',...
                                 'CCP8','CP5','CP3','CP1','CPz',...
                                 'CP2','CP4','CP6','P5','P3',...
                                 'P1','Pz','P2','P4','P6',...
                                 'PO1','PO2','O1','O2'
                                }; 
channel_names =        {  
                        'FC3','FC4','Cz','C3','C4','C5','C6','T7','T8','CCP3','CCP4'      
                        }; 
                    
location = chanlocsseek(channel);
pos_id = [11 15 29 27 31 26 32 25 33 36 39];  
energy = zeros(1,11 );
for ee = 1:11
    energy(ee) = sum(squeeze(imf(ee,:,1)).^2);
end
t_csp_w =sum(energy) ;
t_csp_en=energy/t_csp_w;
figure

t_csp_w1=zeros(59,1);
t_csp_w1(pos_id,1) = t_csp_en*(10e-1) ;
contourIncrement = 10; 
contourLevels = 1; 
scales= [min(t_csp_w1),max(t_csp_w1)];
%please installe EEGLAB toolbox for using topoplot function
topoplot(t_csp_w1,location,'maplimits',scales,'electrodes','off','plotrad',0.6,'headrad',0.6,'gridscale',300,'style','both');

pause(1);
title('left','FontSize',20);
