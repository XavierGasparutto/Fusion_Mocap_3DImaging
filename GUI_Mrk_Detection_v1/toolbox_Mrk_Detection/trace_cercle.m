function[]=trace_cercle(x,y,r,color,line,width)

% tracé des centres de cercles 
b = plot(x,y,'x','color',color);

%définition des équations des cercles
Theta=-pi:0.01:pi;
xc=x+r*cos(Theta);
yc=y+r*sin(Theta);

switch nargin
    case 5
    %tracé des cercles
    a = plot(xc,yc,'color',color,'linestyle',line);
    case 6
    %tracé des cercles
    a = plot(xc,yc,'color',color,'linestyle',line,'linewidth',width);
end

% remove from legend
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b.Annotation.LegendInformation.IconDisplayStyle = 'off';