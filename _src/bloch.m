function [ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec

   gamma = 2*pi*42.577*10^6; %Hz per Tesla 
   %calulate mtOut and mzOut:

    % mT
    mtOut = 1i * exp(1i*angle(B1)).*sin(gamma*dt*abs(B1)) .* mz; % from chiayin's input
%     mtOut = 1i*sin(gamma*dt*B1)*mz;
    mtOut = mtOut + cos(gamma*dt*imag(B1))*real(mt) + 1i*cos(gamma*dt*real(B1))*imag(mt);

    
    % mZ
    mzOut = cos(gamma*dt*B1)*mz;
    mzOut = mzOut + sin(gamma*dt*imag(B1))*real(mt) -    sin(gamma*dt*real(B1))*imag(mt);

    
    %apply B0:
    mtOut = mtOut*exp(1i*gamma*dt*dB0);

   %apply relaxation:
    mtOut = mtOut*exp(-dt/T2);
    mzOut = mzOut*exp(-dt/T1) + 1*(1-exp(-dt/T1));
end