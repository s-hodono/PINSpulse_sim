function [] = showRFandGz(dt,rf,gz)
time = 0:dt:size(rf,2)*dt-dt;

figure
subplot(311)
hold on
plot(time*1000,abs(rf)*10^6,'linewidth',2)
title('RF amplitude')
xlabel('time [ms]') 
ylabel('B1^+ [uT]')
xlim([0 size(rf,2)*dt*1000])
grid on
subplot(312)
hold on
plot(time*1000,angle(rf),'linewidth',2)
title('RF phase')
xlabel('time [ms]') 
ylabel('RF phase [rad]')
ylim([-pi pi])
xlim([0 size(rf,2)*dt*1000])
grid on
subplot(313)
hold on
plot(time*1000,gz*1000,'linewidth',2)
title('gradient')
xlabel('time [ms]') 
ylabel('G_z[mT/m]')
xlim([0 size(rf,2)*dt*1000])
grid on
