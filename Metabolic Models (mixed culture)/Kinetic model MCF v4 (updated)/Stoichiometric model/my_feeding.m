% My_feeding.- Calculation of the feeding stream for any simulation time

function R = my_feeding(t, R)

% This module provides the current conditions of feeding stream
mFd  = R.mFd;          % Field mFd feeding
%For a matrix of feeding stream values, first column time and rest of
%colums values of the feed rate and concentrations
Now = mFd.Now;              % Represents the row index in the feeding matrix
TM = mFd.Time;              % Represents the time at which the row index changes

if Now == length(TM) || t < TM(Now+1)   %If the row index is the last OR if time is less than the next change of row (change of value)
    return
else                                    % Otherwise update the value of the influent stream concentrations, flowrates, and field Now
    % WARNING: CHANGING FEED CONDITIONS DISCONTINUOUSLY CAN LEAD TO
    % DIFFICULTIES IN INTEGRATION. CONSIDER USING INTERPOLATION
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
    R.Inf_ = Inf_;
    R.mFd = mFd;
end