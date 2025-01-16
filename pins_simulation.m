clear all
addpath(genpath('./_src'))
%========================================================================== 
%% settings
%========================================================================== 
Nspin       = 4000;       % number of positions allong the z direction
dt          = 10*10^(-6); % 10us timestep
d           = 1*10^-3;    % slice thickness for PINS (saturation thickness)
D           = 4.5*10^-3;  % slice gap + slice thickness
RFduration  = 7*10^(-3);  % [s]
Nsub        = 14;         % number of PINS sub peaks
FA          = pi/2;       % flip angle [rad]
gamma       = 2*pi*42.577*10^6;    % [rad Hz/T]
prms.gMax   = 40*10^-3;   % Tesla/m, max gradient amplitude
prms.gSlew  = 180;        % T/m/s, max gradient slew
prms.B1max  = 15*10^-6;   % maximum B1 scanner can do[Tesla]
%========================================================================== 
%% preapare RF and gradient blips
%========================================================================== 
[rf,gzblips] = generatePINSpulse(d,D,dt,RFduration,Nsub,FA,prms);
%========================================================================== 
%% preapare slice refocuse gradient
%========================================================================== 
% we dont need refocusing gradient in actual experiment.
Aref   = sum(gzblips)*dt/2; % Half of total Gblips
[gref] = makeGref(Aref,dt,prms);

%========================================================================== 
%% show RF and Gz
%========================================================================== 
gz = [gzblips(1,1:size(rf,2)),gref];
rf = [rf,zeros(1,size(gref,2))];
showRFandGz(dt,rf,gz);
%========================================================================== 
%% prep for Bloch sim
%========================================================================== 
posZ        = zeros(1,Nspin); %variable to hold Nspin positions allong the z direction
mtFinal     = zeros(1,Nspin); %variable to hold the final magnetization calculated for each position
mzFinal     = zeros(1,Nspin); %variable to hold the final magnetization calculated for each position
T1         = 5; % [sec]
T2         = 5; % [sec]
for ii = 1:size(posZ,2)
    posZ(ii)  = (ii-2000)*10^-5; %Distance from iso center in meters
end  
%========================================================================== 
%% off-respnance
%========================================================================== 
dfmax       = 0; % offresonance [Hz]
dB_off      = 2*pi/gamma*dfmax.*linspace(-1,1,size(posZ,2));
%========================================================================== 
%% Bloch sim for PINS saturation profile
%========================================================================== 


for jj = 1:size(posZ,2) %loop over different positions allong the z direction

    %calculate delta b0 at position j and the 1st time step
    dB0        =  dot([0,0,gz(1)],[0,0,posZ(jj)]); 
    [mT,mZ]    =  bloch(dt, dB_off(jj)+dB0,rf(1),T1,T2,0,1);   % start from fully relaxed spin state  

    for tt = 2:size(rf,2)
        %update dB0AtPosJ for the ith time step
        dB0       =  dot([0,0,gz(tt)],[0,0,posZ(jj)]); 
        [mT,mZ]   =  bloch(dt, dB_off(jj)+dB0,rf(tt),T1,T2,mT,mZ);   % start from fully relaxed spin state
    end

    mtFinal(jj)  = mT; %store final m for each j
    mZFinal(jj) = mZ; %store final m for each j
end

%========================================================================== 
%% show
%========================================================================== 

xmax = 20;
figure
subplot(311)
hold on
plot(posZ.*1000,abs(mtFinal),'linewidth',2)
% plot(posZ.*100,abs(mFinal2))
title('Slice profile')
xlabel('position [mm]') 
xlim([ -xmax xmax])
grid on
ylabel('signal intensity [a.u.]')

subplot(312)
hold on
plot(posZ.*1000,abs(mZFinal),'linewidth',2)
title('')
xlabel('position [mm]') 
ylabel('mz')
ylim([-0 1 ])
xlim([ -xmax xmax])
grid on
subplot(313)
hold on
plot(posZ.*1000,angle(mtFinal),'linewidth',2)
title('')
xlabel('position [mm]') 
ylabel('phase [rad]')
ylim([-pi pi])
xlim([ -xmax xmax])
grid on
