function [subIMG,h,w] = f_subIMG(IMG)
% Select sub IMG 
% conditions on rectangle to match dimension of input image in case
% rectangle is out of bounds
[HB,WB] = size(IMG);

rect = getrect;

h = round(rect(2):rect(2)+rect(4)); % Height of rectangle
w = round(rect(1):rect(1)+rect(3)); % Widht of rectangle
if w(1)   < 1;  w = 1:w(end); end 
if w(end) > WB; w = w(1):WB;  end
if h(1)   < 1;  h = 1:h(end); end
if h(end) > HB; h = h(1):HB;  end

subIMG = IMG(h,w);
h = h(1);
w = w(1);