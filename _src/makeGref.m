function [gref] = makeGref(A,dt,prms)
% A           = gradient moment we want to refocus. half of Gslice.
% dt          = timestep [sec]
% prms        = other stuff
% prms.gMax   = gradient amplitude [T/m]
% prms.gSlew  = slew rate
% prms.B1max  = maximum B1 amplitude Scanner allows us to use [T]
% 

gMax        = prms.gMax;
gSlew       = prms.gSlew;


rampidxs   = ceil(gMax/gSlew/dt);
AtriMax    = rampidxs*dt*gMax; % maximum gradient moment of triangle


if AtriMax > A
    t1                = sqrt(A/gSlew); % time for ramp UP
    rampidxs          = ceil(t1/dt);
    gShape            = [(0:rampidxs)./rampidxs fliplr((0:rampidxs)./rampidxs)];
    gShape(rampidxs+1)=[];
else
  nflat       = ceil((A-AtriMax)/gMax/dt);
  gShape      = [(0:rampidxs)./rampidxs ones([1 nflat]) (rampidxs:-1:0)./rampidxs];
end

gref      = -gShape*(A/(sum(gShape)*dt));

end