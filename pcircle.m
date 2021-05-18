function [XYCoords] = pcircle(center, radius, numofpoints)
%--------------------------------------------------------------
% [XYCoords]=PCIRCLE(CENTER,RADIUS,NOP,STYLE)
% This routine draws a circle with center defined as
% a vector CENTER, radius as a scaler RADIUS. NUMOFPOINTS is 
% the number of points on the circle. 
%
%   Usage Examples
%
%   pcircle([1, 3], 3, 360); 
%   pcircle([2, 4], 2, 360);
%
%--------------------------------------------------------------

if (nargin <3)
    error('Please see help for INPUT DATA.');
end

theta   = linspace(0, 2*pi, numofpoints);
rho     = ones(1, numofpoints)*radius;
[x, y]  = pol2cart(theta, rho);
x       = (x + center(1));
y       = (y + center(2));
XYCoords= [x; y];