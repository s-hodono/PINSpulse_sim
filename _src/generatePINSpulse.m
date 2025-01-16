function [rf,gzblips] = generatePINSpulse(d,D,dt,RFduration,Nsub,FA,prms)
% ===== outputs ===========================================================
% rf          = rf waveform with amplitudes 
% gzblips     = slice gradients. blips.
% ===== inputs ============================================================
% d           = slice thicknes [m]
% D           = slice gap + slice thickness[m]
% dt          = timestep [sec]
% RFduration  = RF duration [sec]
% Nsub        = number of subpeaks in RF
% FA          = flip angle [rad]
% prms        = other stuff
% prms.gMax   = gradient amplitude [T/m]
% prms.gSlew  = slew rate
% prms.B1max  = maximum B1 amplitude Scanner allows us to use [T]

%% settings
gamma       = 2*pi*42.577*10^6;    % [rad Hz/T]

gMax        = prms.gMax;
gSlew       = prms.gSlew;

Ablips      = 2*pi/(gamma*D); % gradient moment [T s / m]


rampidxs    = ceil(gMax/gSlew/dt);
AtriMax     = rampidxs*dt*gMax; % maximum gradient moment of triangle


BWTP        = d/D * Nsub;

rf          = zeros(1,round(RFduration/dt));
gzblips     = zeros(1,round(RFduration/dt));
%% get sinc pulse "shape" first
time           = dt:dt:RFduration;
t0             = RFduration/BWTP;
rfSinc         = zeros(size(time));
for ii = 1:round(RFduration/dt)
    rfSinc(ii) = sinc(time(ii)/t0 -(BWTP/2));   
end
hannFilter         = hann(RFduration/dt)';
rfSinc             = hannFilter.*rfSinc;

%% get gradients
% one gradient blip
if AtriMax > Ablips
  t1                = sqrt(Ablips/gSlew); % time for ramp UP
  rampidxs          = ceil(t1/dt);
  gShape            = [(0:rampidxs)./rampidxs fliplr((0:rampidxs)./rampidxs)];
  gShape(rampidxs+1)=[];
  
else
  nflat     = ceil((Ablips-AtriMax)/gMax/dt);
  gShape    = [(0:rampidxs)./rampidxs ones([1 nflat]) (rampidxs:-1:0)./rampidxs];
end

gzblip    = gShape*(Ablips/(sum(gShape)*dt));
oneBlips  = size(gzblip,2);
totBlips  = size(gzblip,2)*(Nsub-1);
totPINS   = size(rfSinc,2)- totBlips;
onePINS   = floor(totPINS/Nsub);

delta     = (oneBlips + onePINS)*dt;
BW        = d/(D*delta);

if onePINS < 0
    disp('===============================================================')
    disp('program stopped.')
    disp('Nsub is too many to RF duration. make RF longer or Nsub smaller')
    disp('===============================================================')
    return;
end


if onePINS*Nsub+totBlips > RFduration/dt
    disp('===============================================================')
    disp('program stopped.')
    disp(['expected RF duration is ',num2str((onePINS*Nsub+totBlips)*dt*1000),'ms'])
    disp('===============================================================')
    return;
end
gztmp        = zeros(onePINS+oneBlips,Nsub-1);

for ii = 1:Nsub-1
    gztmp(:,ii) = [zeros(1,onePINS),gzblip]';
end
idx = (Nsub-1)*(oneBlips + onePINS);
gzblips(1,1:idx) = gztmp(:)';


%% PINS RF
PINScenter_idx = onePINS/2:onePINS+oneBlips:size(rfSinc,2);
% make each sub peak
Subpkeaks      = zeros(onePINS+oneBlips,Nsub);
for ii = 1:Nsub
    Subpkeaks(:,ii)    = rfSinc(ceil(PINScenter_idx(ii))).*[ones(1,onePINS),zeros(1,oneBlips)]';

end
PINS = Subpkeaks(:)';
PINS = PINS(1,1:size(rfSinc,2));


Amp                = FA/(sum(PINS)*dt*gamma);
PINS               = Amp.*PINS;
rf                 = PINS;

%% scaling factor
Amp_sq                = (pi/2)/(1*10^-3*gamma);
Amp_pins              = max(PINS);
scaling_gactor = 0.5/(pi/2) * (pi/2)/(totPINS*dt*1000);
disp(['scaling factor = ', num2str(scaling_gactor)])

%% 
disp(['BWTP = ', num2str(BWTP),', D = ', num2str(D*1000),'mm, BW = ', num2str(BW),', B1 max = ', num2str(max(PINS)*10^6),'uT'])
if Amp > prms.B1max
    disp('===============================================================')
    disp('B1 maybe too large.')
    disp('===============================================================')
    
    return;
end
    
%%

peak_amps = Subpkeaks(1,:);
peak_amps = peak_amps./max(peak_amps);

peak_phs  = rad2deg(angle(peak_amps));




% figure(1)
% subplot(311)
% hold on
% plot(time(1:size(rf,2))*1000,rf*10^6)
% title('RF amplitude')
% xlabel('time [ms]') 
% ylabel('B1^+ [uT]')
% xlim([0 size(rf,2)*dt*1000])
% ylim([-0.2*prms.B1max*10^6 prms.B1max*10^6])
% grid on
% subplot(312)
% hold on
% plot(time(1:size(rf,2))*1000,angle(rf))
% title('RF phase')
% xlabel('time [ms]') 
% ylabel('RF phase [rad]')
% ylim([-pi pi])
% xlim([0 size(rf,2)*dt*1000])
% grid on
% subplot(313)
% hold on
% plot(time(1:size(rf,2))*1000,gzblips(1,1:size(rf,2))*1000)
% title('gradient magnetic field amplitude')
% xlabel('time [ms]') 
% ylabel('G_z[mT/m]')
% xlim([0 size(rf,2)*dt*1000])
% grid on

end