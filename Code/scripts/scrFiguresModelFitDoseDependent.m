% Functions and figures for the model fitting 
% based on the analytical model developed for epigenomeIntegrity project, 
% the functions are here plotted using parameters fitted from the
% experimental data

close all 
% plot properties 

fontSize   = 12;
markerSize = 12;
lineWidth  = 4;

% Plot figues
showHAndDFit        = true;
showSlidingFraction = false;
showSlidingOutOfDR  = true;
showRelativeSliding = true;
showExpansionFactor = false;
showLossInTime      = false;
showRelativeOpening = false;

% Experimental measurements 
% uvc dose (the point u=100 is excluded for now due to irregular
% measurement)
uData = [5 10	15	20	25	30	35	40	45	50	55	60	65	70	75 100];
% histone loss data
% previous 
% hData = [ 10.7143   10.8680   22.1973   24.1895   27.9165   23.1343   36.7809 ...
%     42.0486   38.1288   45.2075   43.6863   44.8139   46.0792   47.6219   48.8158]./100;% 39.5242]./100;
% current 
hData = [10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872...
         21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765...
         42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168]./100;%	40.8353794651
hData = [10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872	21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765	42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168	40.8353794651
]./100;

% hData = dataH./100;
% hData = [ 0.1159    0.1312    0.1963    0.2286    0.2212    0.2137    0.3059    0.3727    0.3735    0.4251    0.4291    0.4285    0.4387    0.4258    0.4564 0.4084];

%___DNA loss data____
% % previous
% dData   =[ 1.5704212005	1.689934217	6.2022046651	11.6868181521	12.8785877917...
%     13.4063551786	18.2744867455	18.0307962375	23.0606990052	23.7519692784...
%     19.9985465308	14.2129016791	18.6451205628	19.8890159764	25.5485722258]./100;%	23.34475654]./100;
% current 
dData = [1.5704212005	1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153...
         12.6333239455	20.0360763966	22.3129622161 22.5107680397	22.7887958612	20.4006799168...
         21.1679155925	22.757261652	26.9902966182]./100;%	26.4974599239
% dData = [ 0.1184    0.1613    0.1148    0.1701    0.1427    0.1354    0.1856    0.2263    0.2318    0.2251    0.2410    0.2366    0.2455    0.2408    0.2699 0.2891];
dData = [1.5704212005	1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153	12.6333239455	20.0360763966	22.3129622161	22.5107680397	22.7887958612	20.4006799168	21.1679155925	22.757261652	26.9902966182	26.4974599239
]./100;

% Analytical solutions of the model for histones and DNA loss vs UV dose
% dData = dataD./100;
T = @(a1,u) 1-exp(-a1.*u);
N = @(a1,a2,u) exp(-a2.*T(a1,u));
R = @(a1,a2,a3,a4,u) 1+a3.*(1-N(a1,a2,u))+a4.*T(a1,u);
% N  = @(a1,a2,u) 1-a2.*(1-exp(-a1.*(u)));% N(U)/N(0)
% R  = @(a1,a2,a3,u) 1+a3.*(1-a2).*(1-exp(-a1.*u))+0.5*a2.*a3.*(1-exp(-2*a1.*u));%   (1-exp(-a1.*u)).*(a3+a2.*a3.*exp(-a1.*u));%R(U)/R(0)
d  = @(a1,a2,a3,a4,u) (R(a1,a2,a3,a4,u)-1)./(R(a1,a2,a3,a4,u));
h  = @(a1,a2,a3,a4,u) 1-(N(a1,a2,u))./R(a1,a2,a3,a4,u) ;%./R(a1,a2,a3,u);
% h  = @(a1,a2,u) 1-exp(-a1*u)./(1+a2.*(1-exp(-a1*u)));
% d  = @(a1,a2,u) a2.*(1-exp(-a1.*u))./(1+a2.*(1-exp(-a1.*u)));
fo = fitoptions('Methods','NonlinearLeastSquares','StartPoint',[0.07 0.05 0.07 0.05],'Lower',[0 0 0],...
    'Robust','off','TolX',1e-18,'TolFun',1e-18,'MaxIter',30000,'MaxFunEval',10000,'Normalize','off','DiffMinChange',0.0001,...
    'DiffMaxChange',0.1);
% ftH = fittype('(a3-(1+a3).*(1-a2+exp(-a1.*x)))./(a3.*a2.*(1-exp(-a1.*x)))','options',fo);
   
ftH = fittype('1-exp(-a2.*(1-exp(-a1.*u)))./(1+a3.*(1-exp(-a2.*(1-exp(-a1.*u))))+a4.*(1-exp(-a1.*u)))',...
              'independent','u','options',fo);
% ftD = fittype('a3.*a2.*exp(-a1.*x)./(1+a3.*a2.*exp(-a1.*x))','options',fo);
%  1-(1-a2.*(1-exp(-a1.*x)))./(1+a3.*a2.*exp(-a1.*x))
% [fitParamsH,gof,output] = fit(uData',hData',ftH,fo);
% % % [fitParamsD]=fit(uData',dData',ftD);
% % % % parameter values
% c1 = fitParamsH.a1;%0.0069;
% c2 = fitParamsH.a2;%0.78;
% c3 = fitParamsH.a3;
% c4 = fitParamsH.a4;


c1 = 0.02;% curve h
c2 = 0.387; % lift h
c3 = 0.415; % lift d
c4 = 0.245; % lift h and d
%-- plot ---
if showHAndDFit
%____ plot histone loss, h
fig1 = figure('Name','histone and DNA fit');
ax1  = axes('Parent',fig1,'NextPlot','Add');
line('Xdata',uData,'YData',h(c1,c2,c3,c4,uData),'Color','r','LineWidth',lineWidth,...
    'Parent',ax1,'DisplayName','histone loss, model');
line('XData',uData,'YData',hData,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','r','MarkerEdgeColor','k',...
    'LineStyle','none','Parent',ax1,'DisplayName','histone loss, exp. data'), 

%____ plot DNA loss, d
line('XData',uData,'YData',d(c1,c2,c3,c4,uData),'Color','g','LineWidth',lineWidth,...
    'DisplayName','DNA loss fraction, model','Parent',ax1,'LineStyle','--');
line('XData',uData,'YData',dData,'Marker','^','Color','g','MarkerSize',markerSize,'MarkerFaceColor','g',...
    'MarkerEdgeColor','k','LineStyle','none',...
    'Parent',ax1,'DisplayName','DNA loss fraction, exp. data'), 

legend(get(ax1,'Children'),'Location','NW')
title(' Model fit to experimental data','Parent',ax1,'FontSize',fontSize), 
xlabel('U.V dose','Parent',ax1,'FontSize',fontSize)
ylabel('Loss fraction','Parent',ax1,'FontSize',fontSize)
set(ax1,'FontSize',fontSize)
end

if showSlidingFraction 
%_____ plot sliding fraction, h-d
fig2 = figure;
ax2  = axes('Parent',fig2);
line('XData',uData,'YData',(hData-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax2,'DisplayName','histone sliding loss, exp.data','MarkerSize',markerSize,'LineStyle','none');
line('XData',uData,'YData',h(c1,c2,c3,uData)-d(c1,c2,c3,uData),'LineWidth',lineWidth,'Parent',ax2,...
    'DisplayName','histone sliding loss, model','Color','k')
xlabel('U.V dose','Parent',ax2,'FontSize',fontSize);
ylabel('h(U)-d(U)','Parent',ax2,'FontSize',fontSize);
title('Fraction of histones lost by sliding','Parent',ax2,'FontSize',fontSize);
legend(get(ax2,'Children'),'Location','NW')
set(ax2,'FontSize',fontSize,'LineWidth',lineWidth)
end

if showSlidingOutOfDR  
%_____ fraction of sliding out of the damage region (h-d)/(1-d)
fig3 = figure; 
ax3  = axes('Parent',fig3,'NextPlot','add');
line('XData',uData,'Ydata',(hData-dData)./(1-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax3,'DisplayName','histone sliding out of DR, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uData,'YData',(h(c1,c2,c3,c4,uData)-d(c1,c2,c3,c4,uData))./(1-d(c1,c2,c3,c4,uData)),'Parent',ax3,...
    'LineWidth',lineWidth,'DisplayName','histone sliding out of DR, model')
% add linear fit for comparison
linearFitModel = fittype('a*x');
[fitValues, fitScore] = fit(uData',((hData-dData)./(1-dData))',linearFitModel,'Robust','LAR','StartPoint',1);
line('XData',uData,'YData',fitValues(uData),'Parent',ax3,'DisplayName','linear fit','Color','r',...
    'LineStyle','--','LineWidth',lineWidth)
title('(H-D)/(1-D), Fraction of histone sliding out of the DR','Parent',ax3,'FontSize',fontSize)
xlabel('U.V dose','Parent',ax3,'FontSize',fontSize);
ylabel('(h(u)-d(u))/(1-d(u)', 'FontSize', fontSize);
legend(get(ax3,'Children'),'Location','NW')
set(ax3,'FontSize',fontSize,'LineWidth',lineWidth)

end

if showRelativeSliding 
%____ relative sliding contribution (h-d)/h
fig4 = figure; 
ax4  = axes('Parent',fig4,'NextPlot','add','FontSize',fontSize);
line('XData',uData,'Ydata',(hData-dData)./hData,'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax4,'DisplayName','relative histone sliding contribution, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uData,'YData',(h(c1,c2,c3,c4,uData)-d(c1,c2,c3,c4,uData))./(h(c1,c2,c3,c4,uData)),'Parent',ax4,...
    'LineWidth',lineWidth,'DisplayName','relative histone sliding contribution, model')
xlabel('U.V dose','Parent',ax4,'FontSize',fontSize);
ylabel('(h(U)-d(U))/h(U)','FontSize',fontSize,'Parent',ax4)
legend(get(ax4,'Children'),'Location','NW');
title('1-D/H, Contibition of sliding to the total histone loss','Parent',ax4,'FontSize',fontSize)
set(ax4,'FontSize',fontSize,'LineWidth',lineWidth)
end

if showExpansionFactor
%___ Expansion factor 
fig5 = figure;
uValues = 0:1:100;
ax5  = axes('Parent',fig5,'NextPlot','add','FontSize',fontSize,'LineWidth',lineWidth);
line('XData',uValues,'YData',R(c1,c2,c3,c4,uValues),'Parent',ax5)
xlabel('U.V dose','Parent',ax5,'FontSize',fontSize)
ylabel('\alpha(U)','Parent',ax5,'FontSize',fontSize);
set(ax5,'FontSize',fontSize,'LineWidth',lineWidth)
end


if showLossInTime
%___histone/DNA loss in time 
fig6  = figure;
kappa = linspace(0,1,100); 
ax6   = axes('Parent',fig6,'NextPlot','add','FontSize',fontSize,'LineWidth',lineWidth);
for uIdx = 1:numel(uData)
 line('XData',kappa,'YData',h(c1,c2,kappa.*uData(uIdx)),'Parent',ax6,'displayName',['histone loss U=' num2str(uData(uIdx))],'Color','r')
 line('XData',kappa,'YData',d(c1,c2,kappa.*uData(uIdx)),'Parent',ax6,'displayName',['DNA loss U=' num2str(uData(uIdx))],'Color','b')
end
xlabel('time','Parent',ax6,'FontSize',fontSize)
ylabel('h(t;U)','Parent',ax6,'FontSize',fontSize);
title('histone and DNA loss in time','Parent', ax6,'FontSize',fontSize);
% legend(get(ax6,'Children'))
set(ax6,'LineWidth',lineWidth)
end


%___calculate the time for which there is gamma percent loss
gamma            = 0.4; % the constant value in (h-d)/h

uValues          = linspace(0,uData(end),100);
expansionSliding = zeros(1,numel(uValues));
expansionOpening = zeros(1,numel(uValues));
tStar            = zeros(1,numel(uValues)); %the fraction of ts for which there is gamma loss of histones

% for uIdx = 1:numel(uValues)
% %     gamma = (hData(uIdx)-dData(uIdx))./hData(uIdx);
%     % the fraction of the time for saturation for which we have gamma loss
%     % of histones
% ts          = uValues(uIdx);% when the fraction equals 1. (at saturation)
% tStar(uIdx) = -(1./(c1*ts)).*log(((c2+1).*(1-gamma*h(c1,c2,c3,ts)))./(1+c2*(1-gamma*h(c1,c2,c3,ts))));
% 
% % plug this time into the equation for the expansion factor
% % and calculate the fraction of the expansion attributed to sliding: (R(tStar,U)-R_0)/R(t_s,U)
% expansionSliding(uIdx) = (R(c1,c2,c3,tStar(uIdx).*ts)-1)./(R(c1,c2,c3,ts)-1);
% % relative contribution of chromatin opening to the expansion 
% expansionOpening(uIdx) =  (R(c1,c2,c3,ts)-R(c1,c2,c3,tStar(uIdx).*ts))./(R(c1,c2,c3,ts)-1);
% end

if showRelativeOpening
%___ Expansion attributed to sliding/ opening 
fig7 = figure('Units','norm');
ax7  = axes('Parent',fig7,'NextPlot','Add');
% line('XData',uData,'YData',S,'LineWidth',lineWidth,'LineStyle','--')
bar(uValues, [expansionSliding;expansionOpening]','stacked','BarWidth',2,'LineStyle','none','Parent',ax7);   
plot(ax7,uValues,expansionSliding,'k','LineWidth',lineWidth)
set(ax7,'YLim',[0 1],'XLim',[0,uData(end)])
xlabel('U.V dose','Parent',ax7,'FontSize',fontSize);
ylabel('fraction of DR expansion','Parent', ax7,'FontSize',fontSize)
title('Relative expansion attibuted to sliding and chromatin opening','Parent',ax7,'FontSize',fontSize);
annotation(fig7,'textbox',...
    [0.38 0.6 0.2 0.1],...
    'String',{'Chromatin opening'},...
    'FitBoxToText','on','LineStyle','none','Color','k','FontSize',fontSize);
annotation(fig7,'textbox',...
    [0.38 0.25 0.2 0.1],...
    'String',{'Histone sliding'},...
    'FitBoxToText','on','LineStyle','none','Color','w','FontSize',fontSize);
set(ax7,'FontSize',fontSize,'LineWidth',lineWidth)
% add a patch 
end

%___Calculate SSEs for fittings 
sseHistone   = sum((h(c1,c2,c3,c4,uData(1:end-1))-hData(1:end-1)).^2);
sseDNA       = sum((d(c1,c2,c3,c4,uData)-dData).^2);
% sseHminusD   = sum((h(c1,c2,c3,uData)-d(c1,c2,c3,uData)-(hData-dData)).^2);
% sseLinearFit = fitScore.sse;
% sseRelativeH = sum(((h(c1,c2,c3,uData)-d(c1,c2,c3,uData))./h(c1,c2,c3,uData)-(hData-dData)./hData).^2);
 sseRelativeH = sum(((h(c1,c2,c3,c4,uData)-d(c1,c2,c3,c4,uData))./(1-d(c1,c2,c3,c4,uData))-(hData-dData)./(1-dData)).^2);
 [sseHistone sseDNA sseRelativeH]
% figure, plot(uData,h(c1,c2,uData)./d(c1,c2,uData))
