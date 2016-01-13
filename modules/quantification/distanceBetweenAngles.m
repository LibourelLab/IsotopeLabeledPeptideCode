function [ distance] = distanceBetweenAngles( a1,a2 )
%distanceBetweenAngles: finds the distance between two angles: used to
%determine the difference between moduli of m/z values
%notes: a1...should be the reference peak (i.e. obserced mass for the
%peptide), and a2 should be the mass you are comparing.


angle = sort([a1,a2]);

distance = angle(2) - angle(1);

if (abs(distance) >180)
    distance = angle(2) -(angle(1)+360);
end



end

