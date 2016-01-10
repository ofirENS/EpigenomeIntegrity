% calculate hazard and survival for signal transfered between concentric
% rings aroung a uv laser focal point. The loss function is described at
% steady state 15 minutes post UV-C. Here we are interested in the material
% passing between layers and not necesserily the end amount. 

R     = 3.5; % [mu/m]
dr    = 0.005;
startRadius = 1e-10;
r     = startRadius:dr:R;
f     = @(r)r;
y     = @(r) (0.15*r.^2) -0.0087.*r -0.64;
intFR = @(r) (0.15/4).*r.^4 -(0.0087/3).*r.^3 -(0.64/2).*r.^2;
h     = @(r)  f(r)./(f(r)-intFR(r));% the hazard function 
s     = zeros(1,numel(r));
tr    = s;
s(1)  = 1;

for rIdx = 2:numel(r)
    s(rIdx)  = exp(-trapz(r(1:rIdx),h(r(1:rIdx))));
    tr(rIdx) = trapz(r(1:rIdx),r(1:rIdx).*s(1:rIdx));
end
tr = tr./trapz(r,tr);
H   = -log(s);
F   = 1-s; % cdf
f   = (1-F).*h(r);
fig = figure('Units','norm');
ax  = axes('Parent',fig);
line('XData',r-r(1),'YData',s,'Color','b','Parent',ax,'DisplayName','survival');
line('XData',r-r(1),'Ydata',h(r),'Color','r','Parent',ax,'DisplayName','hazard');
line('XData',r-r(1),'YData',H,'Color','g','Parent',ax,'DisplayName','cumulative hazard');
line('XData',r-r(1),'YData',F,'Color','k','Parent',ax,'DisplayName','cdf');
line('XData',r-r(1),'YData',f,'Color','m','Parent',ax,'DisplayName','pdf');
% line('XData',r,'YData',r,'Color','y','Parent',ax,'DisplayName','line');
daspect(ax,[1 1 1]);
xlabel(ax,'radius [\mum]')
legend(ax,get(ax,'Children'))
set(ax,'FontSize',25)
set(get(ax,'Children'),'LineWidth',3)
drawnow
