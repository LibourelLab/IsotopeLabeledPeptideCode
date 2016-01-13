function [ R ] = rescaleMatrix(R,r)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

R=R.*r+(1-r);
end

