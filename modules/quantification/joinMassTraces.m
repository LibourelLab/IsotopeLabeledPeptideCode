function [ mout ] = joinMassTraces(mt1,mt2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


if min(mt1.scanNum) > max(mt2.scanNum)
   m1 = mt2;
   m2 = mt1;
else
    m1 = mt1;
    m2 = mt2;
end
    mout = m1;
    mout.scanNum = [m1.scanNum,m2.scanNum];
    mout.rt = [m1.rt,m2.rt];
    mout.idx = [m1.idx,m2.idx];
    mout.mz = [m1.mz,m2.mz];
    mout.intensity = [m1.intensity,m2.intensity];
    mout.noise = [m1.noise,m2.noise];
    mout.fid = [m1.fid,m2.fid];
end

