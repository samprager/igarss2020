function [x_nl,ft] = nonlinearfm(x,fs,B,varargin)
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Samuel Prager
% University of Southern California
% email: sprager@usc.edu
% Created: 2017/11/08 14:18:00; Last Revised: 2017/11/08 14:18:00

%------------- BEGIN CODE --------------

mode = 'time';
if (nargin>3)
    mode = varargin{1};
end

N = numel(x);
T=N/fs;
df = fs/numel(x);
f = [-fs/2:df:fs/2-df];
df = fs/numel(x);
dt = 1/fs;
t = -T/2:dt:(T/2-dt);

if (strcmp(mode,'spectrum'))
    xfft = x;
else
    xfft = fftshift(fft(x));
end

xfft = xfft(abs(f)<=(B/2));
Tg = cumsum(abs(xfft).^2);
% c1*Tg(1)+c2=-T/2
% c2 = -T/2-c1*Tg(1)
% c1*Tg(end)+c2=T/2
% c1*Tg(end) -T/2-c1*Tg(1)= T/2)
% c1(Tg(end)-Tg(1))=T
c1 = T/(Tg(end)-Tg(1));
c2 = -T/2-c1*Tg(1);
Tg = c1*Tg+c2;
frng = f(abs(f)<=(B/2));
if(length(Tg) == length(unique(Tg)))
    ft = interp1(Tg, f(abs(f)<=(B/2)), t);
else
    fprintf('Warning: Tg values not unique')
    [C,ia,ic] = unique(Tg);
    ft = interp1(Tg, f(abs(f)<=(B/2)), t);
end
phi = 2*pi*cumsum(ft)/fs;
phi = phi+pi/2-phi(1);
x_nl = exp(1i*phi);

end
%------------- END OF CODE --------------

