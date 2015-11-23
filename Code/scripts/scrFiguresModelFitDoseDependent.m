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
showSlidingFraction = true;
showSlidingOutOfDR  = true;
showRelativeSliding = true;
showExpansionFactor = false;
showRelativeOpeningDNA  = true;
showRelativeOpeningHistone  = true;

% Experimental measurements 
% uvc dose (the point u=100 is excluded for now due to irregular
% measurement)
uData = [0 5 10	15	20	25	30	35	40	45	50	55	60	65	70	75];% 100];
% histone loss data
% previous 
% hData = [ 10.7143   10.8680   22.1973   24.1895   27.9165   23.1343   36.7809 ...
%     42.0486   38.1288   45.2075   43.6863   44.8139   46.0792   47.6219   48.8158]./100;% 39.5242]./100;
% current 
hData = [10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872...
         21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765...
         42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168]./100;%	40.8353794651
hData = [0 10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872	21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765	42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168 %	40.8353794651
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
dData = [0 1.5704212005	1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153	12.6333239455	20.0360763966	22.3129622161	22.5107680397	22.7887958612	20.4006799168	21.1679155925	22.757261652	26.9902966182	%26.4974599239
]./100;

% Analytical solutions of the model for histones and DNA loss vs UV dose
% dData = dataD./100;
% %--- full (damages) system
T = @(a1,u)  1-exp(-a1.*u);% T(u)/T_max

c1 = 0.021;%0.021;% curve h
c2 = 0.4;  %0.393;% lift h
c3 = 0.39; %0.44; % lift d
c4 = 0.23; %0.23; % lift h and d

% % -- linear (damages) system
% T = @(a1,u) a1.*u; %T/T_max
% 
% c1 = 0.01;%0.023;% curve h
% c2 = 0.54;  %0.35;% lift h
% c3 = 0.2; %0.45; % lift d
% c4 = 0.378; %0.22; % lift h and d

% %-- quadratic (damages) system 
% T  = @(a1,u) a1.*u.^2;
% 
% c1 = 0.0001;%0.023;% curve h
% c2 = 1.0;  %0.35;% lift h
% c3 = 0.2; %0.45; % lift d
% c4 = 0.2; %0.22; % lift h and d

%---
N = @(a1,a2,u) exp(-(a2).*T(a1,u)); % N(u)/N_0
R = @(a1,a2,a3,a4,u) 1+a3.*(1-N(a1,a2,u))+a4.*(T(a1,u)); %R(u)/R_0
d = @(a1,a2,a3,a4,u) (R(a1,a2,a3,a4,u)-1)./(R(a1,a2,a3,a4,u));
h = @(a1,a2,a3,a4,u) 1-(N(a1,a2,u))./R(a1,a2,a3,a4,u) ;%./R(a1,a2,a3,u);

fo = fitoptions('Methods','NonlinearLeastSquares','StartPoint',0.2*[1 1 1 1 ],'Lower',[0 0 0],...
    'Robust','off','TolX',1e-15,'TolFun',1e-15,'MaxIter',30000,'MaxFunEval',10000,'Normalize','off','DiffMinChange',0.000001,...
    'DiffMaxChange',1);
% ftH = fittype('(a3-(1+a3).*(1-a2+exp(-a1.*x)))./(a3.*a2.*(1-exp(-a1.*x)))','options',fo);
   
ftH = fittype('1-exp(-(a2/a1).*(1-exp(-a1.*u)))./(1+a3.*(1-exp(-(a2/a1).*(1-exp(-a1.*u))))+a4.*(1-exp(-a1.*u)))',...
              'independent','u','options',fo);
% ftD = fittype('a3.*a2.*exp(-a1.*x)./(1+a3.*a2.*exp(-a1.*x))','options',fo);
% [fitParamsH,gof,output] = fit(uData',hData',ftH,fo);

% % % [fitParamsD]=fit(uData',dData',ftD);
% % % % parameter values

% ---autofit
% c1 = fitParamsH.a1;%0.0069;
% c2 = fitParamsH.a2;%0.78;
% c3 = fitParamsH.a3;
% c4 = fitParamsH.a4;





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
line('XData',uData,'YData',h(c1,c2,c3,c4,uData)-d(c1,c2,c3,c4,uData),'LineWidth',lineWidth,'Parent',ax2,...
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
% linearFitModel = fittype('a*x');
% [fitValues, fitScore] = fit(uData',((hData-dData)./(1-dData))',linearFitModel,'Robust','LAR','StartPoint',1);
% line('XData',uData,'YData',fitValues(uData),'Parent',ax3,'DisplayName','linear fit','Color','r',...
%     'LineStyle','--','LineWidth',lineWidth)
title('Fraction of histone sliding out of the IDR','Parent',ax3,'FontSize',fontSize)
xlabel('U.V dose','Parent',ax3,'FontSize',fontSize);
ylabel('(H(u)-D(u))/(1-D(u)', 'FontSize', fontSize);
legend(get(ax3,'Children'),'Location','NW')
set(ax3,'FontSize',fontSize,'LineWidth',lineWidth)

end

if showRelativeSliding 
%____ relative sliding contribution (h-d)/h
fig4 = figure; 
ax4  = axes('Parent',fig4,'NextPlot','add','FontSize',fontSize);
line('XData',uData,'Ydata',1 -(dData)./hData,'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax4,'DisplayName','relative histone sliding contribution, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uData,'YData',1-(d(c1,c2,c3,c4,uData))./(h(c1,c2,c3,c4,uData)),'Parent',ax4,...
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


if showRelativeOpeningDNA
fig6 = figure;
% this is the DNA loss due to opening 
% show A_openning/A(u)
uVals = 0:3:max(uData);
ax6 = axes('Parent', fig6,'NextPlot','Add','FontSize',fontSize,'LineWidth',lineWidth);
title('Relative contribution to DNA loss','FontSize',fontSize,'Parent',ax6)
xlabel(ax6,'UV dose','FontSize',fontSize)
ylabel(ax6,'Fraction of loss','FontSize',fontSize)
relativeOpening = [(c4.*T(c1,uVals))./(c3.*(1-exp(-c2.*T(c1,uVals)))+c4.*T(c1,uVals));...
    1-(c4.*T(c1,uVals))./(c3.*(1-exp(-c2.*T(c1,uVals)))+c4.*T(c1,uVals))]';
bar(uVals,relativeOpening,1,'Stacked')
c = get(ax6,'Children');
set(c(1),'LineStyle','none','FaceColor','y');
set(c(2),'LineStyle','none')
 legend('opening','sliding');
 annotation(fig6,'textbox',...
    [0.404834260977118 0.714880332986472 0.226581941867656 0.100936524453697],...
    'String',{'Sliding contribution'},...
    'FontSize',25,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'LineWidth',2.5);
annotation(fig6,'textbox',...
    [0.402730531520393 0.369529983792544 0.270322620519159 0.11345218800648],...
    'String',{'Opening contribution'},...
    'FontSize',25,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'LineWidth',2.5,...
    'Color',[1 1 1]);
end

if showRelativeOpeningHistone  
    % histone loss due to opening
  fig7 = figure;
  uVals = 0:0.5:max(uData);
  ax7 = axes('Parent',fig7,'NextPlot','Add');
  relativeOpeningH = [(1-(1./(1+c4.*T(c1,uVals))))./(1-N(c1,c2,uVals)./R(c1,c2,c3,c4,uVals));...
      1-(1-(1./(1+c4.*T(c1,uVals))))./(1-N(c1,c2,uVals)./R(c1,c2,c3,c4,uVals))];
  bar(uVals,relativeOpeningH',1.5,'Stacked')
  title('Relative contribution to histone loss','FontSize',fontSize)
  xlabel(ax7,'UV dose','FontSize',fontSize); 
  ylabel(ax7,'Fraction of loss','FontSize',fontSize)
  set(ax7,'FontSize',fontSize,'LineWidth',lineWidth)
  c = get(ax7,'Children');
  set(c(1),'LineStyle','none','FaceColor','y');
  set(c(2),'LineStyle','none')
  annotation(fig7,'textbox',...
    [0.398404202719404 0.170872887242418 0.264760197775033 0.147995889003083],...
    'String',{'Opening contribution'},...
    'FontSize',25,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'LineWidth',2.5,...
    'Color',[1 1 1]);
annotation(fig7,'textbox',...
    [0.412619283065513 0.621788283658787 0.238802224969098 0.112024665981501],...
    'String',{'Sliding contribution'},...
    'FontSize',25,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'LineWidth',2.5);
%   legend('opening','sliding');
end

%___Calculate SSEs and R2 for fittings 
sseHistone   = sum((h(c1,c2,c3,c4,uData(1:end-1))-hData(1:end-1)).^2);
rsHistone    = 1-sseHistone./sum((hData-mean(hData)).^2);
sseDNA       = sum((d(c1,c2,c3,c4,uData)-dData).^2);
rsDNA        = 1- sseDNA./sum((dData-mean(dData)).^2);
sseSlidingOutOfIDR = sum(((h(c1,c2,c3,c4,uData)-d(c1,c2,c3,c4,uData))./(1-d(c1,c2,c3,c4,uData)) -(hData-dData)./(1-dData)).^2);
rsSlidingOutOfIDR  = 1- sseSlidingOutOfIDR./sum(((hData-dData)./(1-dData) - mean((hData-dData)./(1-dData))).^2);
[rsDNA rsHistone,rsSlidingOutOfIDR]

