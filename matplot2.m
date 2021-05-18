% This function opens an image in joint time-freq. domain
% Inputs====> x: matrix in time-freq domain
%             r: dynamic range of the screen
%             F: Frequency vector
%             T: Time Vector
function [pp, p]=matplot2(T,F,x,r)
b=max(max(abs(x)));
ra=b/(10^(r/20));
p=x.*(abs(x)>=ra)+ra*ones(size(x)).*(abs(x)<ra);
pp=20*log10(abs(p)/b);
colormap(jet(256)),imagesc(T,F,pp);
