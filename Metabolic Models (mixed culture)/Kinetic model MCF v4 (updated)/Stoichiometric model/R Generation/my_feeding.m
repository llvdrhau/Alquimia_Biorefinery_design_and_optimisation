% My_feeding.- Calculation of the feeding stream for any simulation time

function R = my_feeding(t, idR, R)

% This module provides the current conditions of feeding stream
mFd  = R(idR).mFd;  
Now = mFd.Now;
TM = mFd.Time;

if Now == length(TM) || t < TM(Now+1)
    return
else
    % From the Excel feeding program mFd, depending on the current time value
    % appropriate feed composition is assigned
    InfV = mFd.FdV(Now+1,:); 
    mFd.Now = Now+1;
    Inf_.Qgas = InfV(1);
    Inf_.Q = InfV(2);
    Inf_.InfV = InfV';

    % Structure with names updated
    for i=1:length(InfV),
        Inf_.(char(mFd.FdNames(i))) = InfV(i);
    end

    % Output to the global R
    R(idR).Inf_ = Inf_;
    R(idR).mFd = mFd;
end