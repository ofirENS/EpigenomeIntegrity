% Functions and figures for the model fitting 
% based on the analytical model developed for epigenomeIntegrity project, 
% the functions are here plotted using parameters fitted from the
% experimental data

close all 
% plot properties 

fontSize   = 40;
markerSize = 12;
lineWidth  = 4;

% Experimental measurements 
% uvc dose (the point u=100 is excluded for now due to irregular
% measurement)
uData = [ 5	10	15	20	25	30	35	40	45	50	55	60	65	70	75];%	100];
% histone loss data
hData = [ 10.7143   10.8680   22.1973   24.1895   27.9165   23.1343   36.7809 ...
    42.0486   38.1288   45.2075   43.6863   44.8139   46.0792   47.6219   48.8158]./100;% 39.5242]./100;
% DNa loss data
dData   =[ 1.5704212005	1.689934217	6.2022046651	11.6868181521	12.8785877917...
    13.4063551786	18.2744867455	18.0307962375	23.0606990052	23.7519692784...
    19.9985465308	14.2129016791	18.6451205628	19.8890159764	25.5485722258]./100;%	23.34475654]./100;
      
% Analytical solutions of the model for histones and DNA loss vs UV dose
h  = @(a1,a2,u) 1-exp(-a1*u)./(1+a2.*(1-exp(-a1*u)));
d  = @(a1,a2,u) a2.*(1-exp(-a1.*u))./(1+a2.*(1-exp(-a1.*u)));
R = @(a1,a2,u) 1+a2.*(1-exp(-a1.*u));
% parameter values
c1 = 0.007;% 0.0070;
c2 = 0.7;

%____ plot histone loss, h
fig1 = figure();
ax1  = axes('Parent',fig1,'NextPlot','Add');
line('Xdata',uData,'YData',h(c1,c2,uData),'Color','r','LineWidth',lineWidth,...
    'Parent',ax1,'DisplayName','histone loss model');
line('XData',uData,'YData',hData,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','r','MarkerEdgeColor','k',...
    'LineStyle','none','Parent',ax1,'DisplayName','histone loss exp. data'), 

%____ plot DNA loss, d
line('XData',uData,'YData',d(c1,c2,uData),'Color','g','LineWidth',lineWidth,...
    'DisplayName','DNA loss fraction, model','Parent',ax1,'LineStyle','--');
line('XData',uData,'YData',dData,'Marker','^','Color','g','MarkerSize',markerSize,'MarkerFaceColor','g',...
    'MarkerEdgeColor','k','LineStyle','none',...
    'Parent',ax1,'DisplayName','DNA loss fraction, exp. data'), 

legend(get(ax1,'Children'),'Location','NW')
title(' Model fit to experimental data','Parent',ax1,'FontSize',fontSize), 
xlabel('U.V dose','Parent',ax1,'FontSize',fontSize)
ylabel('Loss fraction','Parent',ax1,'FontSize',fontSize)
set(ax1,'FontSize',fontSize)

%_____ plot sliding fraction, h-d
fig2 = figure;
ax2  = axes('Parent',fig2);
line('XData',uData,'YData',(hData-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax2,'DisplayName','histone sliding loss, exp.data','MarkerSize',markerSize,'LineStyle','none');
line('XData',uData,'YData',h(c1,c2,uData)-d(c1,c2,uData),'LineWidth',lineWidth,'Parent',ax2,...
    'DisplayName','histone sliding loss, model','Color','k')
xlabel('U.V dose','Parent',ax2,'FontSize',fontSize);
ylabel('h(U)-d(U)','Parent',ax2,'FontSize',fontSize);
title('Fraction of histones lost by sliding','Parent',ax2,'FontSize',fontSize);
legend(get(ax2,'Children'),'Location','NW')
set(ax2,'FontSize',fontSize,'LineWidth',lineWidth)


%_____ fraction of sliding out of the damage region (h-d)/(1-d)
fig3 = figure; 
ax3  = axes('Parent',fig3,'NextPlot','add');
line('XData',uData,'Ydata',(hData-dData)./(1-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax3,'DisplayName','histone sliding out of DR, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uData,'YData',(h(c1,c2,uData)-d(c1,c2,uData))./(1-d(c1,c2,uData)),'Parent',ax3,...
    'LineWidth',lineWidth,'DisplayName','histone sliding out of DR, model')
% add linear fit for comparison
linearFitModel = fittype('a*x');
[fitValues, fitScore] = fit(uData',((hData-dData)./(1-dData))',linearFitModel,'Robust','LAR');
line('XData',uData,'YData',fitValues(uData),'Parent',ax3,'DisplayName','linear fit','Color','r',...
    'LineStyle','--','LineWidth',lineWidth)
title('Fraction of histone sliding out of the DR','Parent',ax3,'FontSize',fontSize)
xlabel('U.V dose','Parent',ax3,'FontSize',fontSize);
ylabel('(h(u)-d(u))/(1-d(u)', 'FontSize', fontSize);
legend(get(ax3,'Children'),'Location','NW')
set(ax3,'FontSize',fontSize,'LineWidth',lineWidth)

%____ relative sliding contribution (h-d)/h
fig4 = figure; 
ax4  = axes('Parent',fig4,'NextPlot','add','FontSize',fontSize);
line('XData',uData,'Ydata',(hData-dData)./(hData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax4,'DisplayName','relative histone sliding contribution, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uData,'YData',(h(c1,c2,uData)-d(c1,c2,uData))./(h(c1,c2,uData)),'Parent',ax4,...
    'LineWidth',lineWidth,'DisplayName','relative histone sliding contribution, model')
xlabel('U.V dose','Parent',ax4,'FontSize',fontSize);
ylabel('(h(U)-d(U))/h(U)','FontSize',fontSize,'Parent',ax4)
legend(get(ax4,'Children'),'Location','NW');
title('Contibition of sliding to the total histone loss','Parent',ax4,'FontSize',fontSize)
set(ax4,'FontSize',fontSize,'LineWidth',lineWidth)

%___ Expansion factor 
fig5 = figure;
ax5  = axes('Parent',fig5,'NextPlot','add','FontSize',fontSize,'LineWidth',lineWidth);
line('XData',uData,'YData',R(c1,c2,uData),'Parent',ax5)
xlabel('U.V dose','Parent',ax5,'FontSize',fontSize)
ylabel('\alpha(U)','Parent',ax5,'FontSize',fontSize);
set(ax5,'FontSize',fontSize,'LineWidth',lineWidth)

%___histone loss in time 
fig6  = figure;
kappa = linspace(0,1,100); 
ax6  = axes('Parent',fig6,'NextPlot','add','FontSize',fontSize,'LineWidth',lineWidth);
for uIdx = 1:numel(uData)
 line('XData',kappa,'YData',h(c1,c2,kappa.*uData(uIdx)),'Parent',ax6,'displayName',['histone loss U=' num2str(uData(uIdx))],'Color','r')
 line('XData',kappa,'YData',d(c1,c2,kappa.*uData(uIdx)),'Parent',ax6,'displayName',['DNA loss U=' num2str(uData(uIdx))],'Color','b')
end
xlabel('time','Parent',ax6,'FontSize',fontSize)
ylabel('h(t;U)','Parent',ax6,'FontSize',fontSize);
% legend(get(ax6,'Children'))
set(ax6,'LineWidth',lineWidth)

%___calculate the time for which there is gamma percent loss
gamma = 0.6; % the constant value in (h-d)/h
uValues = linspace(0,uData(end),100);
S     = zeros(1,numel(uValues));
tStar = zeros(1,numel(uValues));

for uIdx = 1:numel(uValues)
%     gamma = (hData(uIdx)-dData(uIdx))./hData(uIdx);
    % the fraction of the time for saturation for which we have gamma loss
    % of histones
tStar(uIdx) = -(1./(c1*uValues(uIdx))).*log(((c2+1).*(1-gamma*h(c1,c2,uValues(uIdx))))./(1+c2-gamma*h(c1,c2,uValues(uIdx))));
% plug this time into the equation for the expansion factor
% and calculate the fraction of the expansion attributed to sliding: R(tStar,U)/R(t_s,U)
S(uIdx) = R(c1,c2,tStar(uIdx).*uValues(uIdx))./R(c1,c2,uValues(uIdx));
end

S1 = 1-S;% relative contirbution of chromatin opening to the expansion 

%___ expansion attributed to sliding/ opening 
fig7 = figure('Units','norm');
ax7  = axes('Parent',fig7,'NextPlot','Add');
% line('XData',uData,'YData',S,'LineWidth',lineWidth,'LineStyle','--')
bar(uValues, [S;S1]','stacked','BarWidth',2,'LineStyle','none','Parent',ax7);   
plot(ax7,uValues,S,'k','LineWidth',lineWidth)
set(ax7,'YLim',[0 1],'XLim',[0,uData(end)])
xlabel('U.V dose','Parent',ax7,'FontSize',fontSize);
ylabel('fraction of DR expansion','Parent', ax7,'FontSize',fontSize)
title('Relative expansion attibuted to sliding and chromatin opening','Parent',ax7,'FontSize',fontSize);
annotation(fig7,'textbox',...
    [0.58 0.8 0.2 0.1],...
    'String',{'Chromatin opening'},...
    'FitBoxToText','on','LineStyle','none','Color','w','FontSize',fontSize);
annotation(fig7,'textbox',...
    [0.38 0.45 0.2 0.1],...
    'String',{'Histone sliding'},...
    'FitBoxToText','on','LineStyle','none','Color','w','FontSize',fontSize);
set(ax7,'FontSize',fontSize,'LineWidth',lineWidth)
% add a patch 

%___Calculate SSE 
sseHistone   = sum((h(c1,c2,uData(1:end-1))-hData(1:end-1)).^2);
sseDNA       = sum((d(c1,c2,uData)-dData).^2);
sseHminusD   = sum((h(c1,c2,uData)-d(c1,c2,uData)-(hData-dData)).^2);
sseRelativeH = sum(((h(c1,c2,uData)-d(c1,c2,uData))./h(c1,c2,uData)-(hData-dData)./hData).^2);


