function [trout,envsm_ave,envsm] = trace_equal(trin,sampint,op_length,trip)
% Syntax
%  trout=aec(trin);
%  trout=aec(trin,sampint,op_length,trip)
%  trout=aec(trin,sampint,op_length,trip,profileflag)
%
% Description 
%   AEC performs an automatic amplitude adjustment.
%
% Method
%   1) Compute Hilbert envelope of the input trace TRIN
%   2) Convolve envelope with triangular smoother of half-length
%      OP_LENGTH
%   3) Divide input trace by smoothed envelope
%   4) Balance the output to a maximum of 1.0
%
% Inputs
%   trin      = input trace or gather of traces.
%   t or sampint   = sample interval for trin 
%               For backwards compatibility, if sampint is supplied as a time coordinate 
%               vector, the difference of the first two elements is used as
%               the sample interval.
%   **********  Default is 0.001 seconds.  (1 millisecond) ********
%   op_length = half-length of triangular smoother in seconds
%   **********  default is 1/8th of the trace length *******
%               ***** must be less than half the trace length *****
%   trip      = front end time before which the smoothed envelope is
%               set to a constant
%   **********  default is op_length/10 *******

% Outputs
%   trout     = output trace or gather of traces
%   envsm     = smoothed hilbert envelope of trace or traces
%   ma        = maximum absolut value of the corretect trace befor
%               normalization to 1.
% NOTES:
%   1) To remove trace normalization to 1, trout_unnorm=trout*ma; (single
%           trace) or trout_unnorm(:,k)=trout(:,k)*ma(a); (multitrace)
%   2) To recover the input trace: trin_recovered=trout*ma.*envsm; (single
%           trace) or trin_recovered(:,k)=trout(:,k)*ma.*envsm(:,k) (multitrace)
% NOtice: Modified from CREWS


% set defaults
if nargin < 2 || isempty(sampint)
    sampint = 0.001;   % sample interval in seconds
else
    % Backwards compatibility with time coordinate vector
    if length(sampint) > 1
        sampint = sampint(2) - sampint(1);
    end
end
if nargin < 3 || isempty(op_length)
    op_length = sampint*length(trin)/8;
end
if nargin < 4 || isempty(trip);
    trip=op_length/10.;
end

% Number of traces
[nt,ntr] = size(trin);

% average smoothed evenlope
[envsm_avg,envstack]=average_envelope(trin,sampint,op_length,trip);

% smoothed evenlopes for each trace
trout0=zeros(size(trin));
envsm=trout0;
for i=1:ntr
    [trout0(:,i),envsm(:,i)] = aec_vector(trin(:,i)',sampint,op_length,trip);
end

% change traces' magnitides according to the differences of their evnsm 
% between envsm_avg
envmul = zeros(size(envsm));
for i = 1:ntr
    envmul(:,i) = smooth(envsm(:,i)./envsm_avg);
end
trout = trin./envmul;

function [trout_,envsm] = aec_vector(trin_,sampint_,op_length_,trip_)
% the original trace is trin_=trout*ma.*envsm

% double the operator length
op2=op_length_*2;
% form new trace padded to a power of 2
trinew=padpow2(trin_,0);
% compute the envelope
env=abs(hilbm(trinew));
env=env(1:length(trin_));
% compute the smoothed envelope
nop=round(op2/sampint_)+1;
envsm=conv(env,triang(nop));
% grab the central length(trin) samples
envsm=envsm(round(nop/2):length(trin_)+round(nop/2)-1);
% stabilize the envelope at the ends
ntrip=round(trip_/sampint_)+1;

%for some reason, envsm(ntrip+1:length(envsm) is sometimes a column vector
%and sometimes a row vector. KWH, 2014
envsm_part = envsm(ntrip+1:length(envsm));
m = size(envsm(ntrip+1:length(envsm)),1);
if m > 1
    envsm_part = envsm_part';
end

envsm=[envsm(ntrip)*ones(1,ntrip) envsm_part];
envsm=[envsm(1:length(envsm)-ntrip) envsm(length(envsm)-ntrip)*ones(1,ntrip)];
% correct the trace
if(sum(abs(envsm))==0)
    trout_=zeros(size(trin_));
%     ma=0;
else
    trout_=trin_./envsm;
%     % balance the output to have a maximum of 1
%     ma=max(abs(trout_));
%     %ma=1;
%     trout_=trout_/ma;
end


% balance the output to have the same mean power as input
%trout=balans(trout,trin);
end

function [envsm,envstack]=average_envelope(trin_,sampint_,op_length_,trip_)
%get the average envelope and smooth it
[nt_,ntr_]=size(trin_);
envstack=zeros(nt_,1);
for k=1:ntr_
%     tr=padpow2(trin_(:,k),0);%make sure its a power of 2 for hilbert
    envstack=envstack+env(trin_(:,k));
end
envstack=envstack(1:nt_)/ntr_;
% double the operator length
op2=op_length_*2;
% compute the smoothed envelope
nop=round(op2/sampint_)+1;
envsm=conv(envstack,triang(nop))/sum(triang(nop));
% grab the central length(trin) samples
% envsm=envsm(round(nop/2):length(trin_)+round(nop/2)-1);
envsm=envsm(round(nop/2):nt_+round(nop/2)-1); % modified by mbw at 20191226
% stabilize the envelope at the ends
ntrip=round(trip_/sampint_)+1;

%for some reason, envsm(ntrip+1:length(envsm) is sometimes a column vector
%and sometimes a row vector. KWH, 2014
envsm_part = envsm(ntrip+1:length(envsm));
m = size(envsm(ntrip+1:length(envsm)),1);
if m > 1
    envsm_part = envsm_part';
end

envsm=[envsm(ntrip)*ones(1,ntrip) envsm_part];
envsm=[envsm(1:length(envsm)-ntrip) envsm(length(envsm)-ntrip)*ones(1,ntrip)]';
end

function w=triang(n)

if(iseven(n))
    n2=n/2;
    d=1/n2;
    val=d*(1:n2)-d/2;
    %val=linspace(0,1,n2+2);
    w=[val val(n2:-1:1)];
    w=w';    
else
    n2=ceil(n/2);
    d=1/n2;
    val=d*(1:n2-1);
    w=[val 1 val(n2-1:-1:1)];
    w=w';
end
end

end