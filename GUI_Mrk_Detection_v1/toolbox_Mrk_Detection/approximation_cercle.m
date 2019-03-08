function [xc,yc,rc]=approximation_cercle(x,y)

% calcul du cercle de régression par les moindres carrés
abc=[x y ones(length(x),1)]\-(x.^2+y.^2);
a=abc(1); b=abc(2); c=abc(3);
xc=-a/2;
yc=-b/2;
rc=sqrt((xc^2+yc^2)-c);