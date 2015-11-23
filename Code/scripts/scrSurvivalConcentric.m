% calculate hazard and survival for signal transfered between concentric rings 
R     = 3.5;
dr    = 0.01;
r     = 2:dr:R;
f     = @(r)r;
y     = @(r) (0.15*r.^2) -0.0087.*r -0.64;
intFR = @(r) (0.15/4).*r.^4 -(0.0087/3).*r.^3 -(0.64/2).*r.^2;
h     = @(r)  intFR(r)./(f(r)-intFR(r))+1;% the hazard function 
s     = zeros(1,numel(r));
tr    = s;
s(1)  = 1;

for rIdx = 2:numel(r)
    s(rIdx)  = exp(-trapz(r(1:rIdx),h(r(1:rIdx))));
    tr(rIdx) = trapz(r(1:rIdx),r(1:rIdx).*s(1:rIdx));
end

H   = -log(s);
F   = 1-s; % cdf
f   = (1-F).*h(r);
fig = figure;
ax  = axes('Parent',fig);
line('XData',r,'YData',s,'Color','b','Parent',ax,'DisplayName','survival');
line('XData',r,'Ydata',h(r),'Color','r','Parent',ax,'DisplayName','hazard');
line('XData',r,'YData',H,'Color','g','Parent',ax,'DisplayName','cumulative hazard');
line('XData',r,'YData',F,'Color','k','Parent',ax,'DisplayName','cdf');
line('XData',r,'YData',f,'Color','m','Parent',ax,'DisplayName','pdf');
% line('XData',r,'YData',r,'Color','y','Parent',ax,'DisplayName','line');
daspect(ax,[1 1 1]);
legend(ax,get(ax,'Children'))
