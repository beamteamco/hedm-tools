function Rp = ErrorRp(yData, yCalc, varargin)
% ErrorRp - calculates residual of the pattern.  for more details
% look at Rietvald Method (Young).
%
%   USAGE:
%
%   Rp  = ErrorRp(yData, yCalc)
%
%   INPUT:
%   yData
%       experimental intensity data
%
%   yCalc
%       calculated intensity
%
%   OUTPUT:
%
%   Rp
%       residual

% default options
optcell = {...
    'Threshold', 0.0001, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

i1  = yData > opts.Threshold;
i2  = yData < -opts.Threshold ;
i   = i1 | i2;

if sum(i) == 0
    Rp  = 0;
else
    Rp  = sum((yData(i) - yCalc(i)).^2)/sum(yData(i).^2);
    Rp  = Rp.^.5;
end