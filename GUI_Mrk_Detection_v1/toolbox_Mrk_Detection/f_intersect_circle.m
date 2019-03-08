function [i] = f_intersect_circle(o1,r1,o2,r2)
% Xavier Gasparutto - 2019 - Willy Taillard Kinesiology Laboratory, HUG, UNIGE, Geneva

% Find points of intersection of circle 1 and 2
% Circle 1: origin o1, radius r1
% Circle 2: origin o2, radius r2
% Use triangle o1,o2,i w/ i = intersection point
r3 = norm(o2 - o1);                                                        % Distance between circle centres
u  = (o2-o1)/r3;                                                           % Unit vector between circle centres
% law of cosine - angle between O1O2 and O1I
a = acos((r1*r1 + r3*r3 - r2*r2)/(2*r1*r3));                               % angle can be >0 or <0
R1 = [cos(a) -sin(a);sin(a) cos(a)];                                       % With Positive angle
R2 = R1';                                                                  % With Negative angle
% Intersection point - 2 solutions
i(1,:) = o1 + R1*u*r1;                                                     % Solution 1
i(2,:) = o1 + R2*u*r1;                                                     % Solution 2