function [f] = gaussmf(x,sig,c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
f = exp(-(x-c).^2/(2*sig.^2));
end

