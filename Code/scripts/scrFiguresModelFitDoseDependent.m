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
showExpansionFactor = true;
showRelativeOpeningDNA      = true;
showRelativeOpeningHistone  = true;

% Experimental measurements 
% uvc dose (the point u=100 is excluded for now due to irregular
% measurement)
% uData = [0 5 10	15	20	25	30	35	40	45	50	55	60	65	70	75 100];
uData = [0  10	15	20	25	30	35	40	45	50	55	60	65	70	75 100];
%__ histone loss data___
% hData = [0 10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872	21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765	42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168 	40.8353794651
% ]./100;
% excluding measurement at 5 msec
hData = [0 	10.8220788165	14.4014983755	20.8225447327	21.2024074872	21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765	42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168 	40.8353794651
]./100;

%___DNA loss data____
% dData = [0 1.5704212005	1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153	12.6333239455	20.0360763966	22.3129622161	22.5107680397	22.7887958612	20.4006799168	21.1679155925	22.757261652	26.9902966182	26.4974599239
% ]./100;
% excluding measurement at 5 msec
dData = [0 	1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153	12.6333239455	20.0360763966	22.3129622161	22.5107680397	22.7887958612	20.4006799168	21.1679155925	22.757261652	26.9902966182	26.4974599239
]./100;

% Analytical solutions of the model for histones and DNA loss vs UV dose
T       = @(a1,u)  (1-exp(-a1.*u.^2));% T(u)/T_max
N_slide = @(a1,a2,u) 1-exp(-a2.*T(a1,u));
R       = @(a1,a2,a3,a4,u)  1+a3.*N_slide(a1,a2,u)+a4.*T(a1,u);%1+a4.*N_total(a1,a2,a3,u);%+(a4).*(T(a1,u))); %R(u)/R_0
d       = @(a1,a2,a3,a4,u) 1-1./R(a1,a2,a3,a4,u);%(((R(a1,a2,a3,a4,u))-1+N_open(a1,a2,a3,u))./(R(a1,a2,a3,a4,u)));%+(T(a1,u))./R(a1,a2,a3,a4,u);
h       = @(a1,a2,a3,a4,u) ((d(a1,a2,a3,a4,u) +(N_slide(a1,a2,u))./R(a1,a2,a3,a4,u)));%d(a1,a2,a3,a4,u)+(T(a1,u)-N(a1,a2,u))./R(a1,a2,a3,a4,u) ;%./R(a1,a2,a3,u);


% % % ---autofit
% opt       = optimset('TolFun',1e-15,'TolX',1e-15,'MaxIter',1e7,'MaxFunEvals',1e7,'TolCon',1e-15,'Hessian','bfgs',...
%     'Diagnostics','off');
% % run several tests
% numTests  = 1;
% fitParams = zeros(numTests,4); 
% fval      = zeros(numTests,1);
% for tIdx = 1:numTests
% [fitParams(tIdx,:),fval(tIdx),exitFlag,output]=...
%     fmincon(@FitDandH,.1*rand(1,4),-1*eye(4),zeros(4,1),[],[],zeros(4,1),10*ones(4,1),[],opt);
% end
% [~,pl] = min(fval);
%  fitParams = fitParams(pl,:);
% c1 = fitParams(1);
% c2 = fitParams(2);
% c3 = fitParams(3);
% c4 = fitParams(4);

% % for square u in T
c1 = 0.0008;%0.036;% curve h
c2 = 0.31;  %0.35;% lift h
c3 = 0.29; %0.24; % contribution of sliding 
c4 = 0.24; %0.48; % contribution of pushing

% % % for linear u in T
% c1 = 0.024;% rate of damage accumulation
% c2 = 0.36; % lift h
% c3 = 0.26; % contribution of sliding 
% c4 = 0.24; % contribution of pushing

uVals = 0.001:0.5:max(uData);


% Plot
if showHAndDFit
%____ Plot histone loss, h
fig1 = figure('Name','nucleosome and DNA fit');
ax1  = axes('Parent',fig1,'NextPlot','Add');
line('Xdata',uVals,'YData',h(c1,c2,c3,c4,uVals),'Color','r','LineWidth',lineWidth,...
    'Parent',ax1,'DisplayName','nucleosome loss, model');
line('XData',uData,'YData',hData,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','r','MarkerEdgeColor','k',...
    'LineStyle','none','Parent',ax1,'DisplayName','nucleosome loss, exp. data'), 

%____ Plot DNA loss, d
line('XData',uVals,'YData',d(c1,c2,c3,c4,uVals),'Color','g','LineWidth',lineWidth,...
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
%_____ Plot sliding fraction, h-d
fig2 = figure;
ax2  = axes('Parent',fig2);
line('XData',uData,'YData',(hData-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax2,'DisplayName','nucleosome sliding loss, exp.data','MarkerSize',markerSize,'LineStyle','none');
line('XData',uVals,'YData',h(c1,c2,c3,c4,uVals)-d(c1,c2,c3,c4,uVals),'LineWidth',lineWidth,'Parent',ax2,...
    'DisplayName','nucleosome sliding loss, model','Color','k')
xlabel('U.V dose','Parent',ax2,'FontSize',fontSize);
ylabel('h(U)-d(U)','Parent',ax2,'FontSize',fontSize);
title('Fraction of nucleosomes lost by sliding','Parent',ax2,'FontSize',fontSize);
legend(get(ax2,'Children'),'Location','NW')
set(ax2,'FontSize',fontSize,'LineWidth',lineWidth)
end

if showSlidingOutOfDR  
%_____ fraction of sliding out of the damage region (h-d)/(1-d)
fig3 = figure; 
ax3  = axes('Parent',fig3,'NextPlot','add');
line('XData',uData,'Ydata',(hData-dData)./(1-dData),'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax3,'DisplayName','nucleosome sliding out of DR, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uVals,'YData',(h(c1,c2,c3,c4,uVals)-d(c1,c2,c3,c4,uVals))./(1-d(c1,c2,c3,c4,uVals)),'Parent',ax3,...
    'LineWidth',lineWidth,'DisplayName','nucleosome sliding out of DR, model')

title('Fraction of nucleosome sliding out of the IDR','Parent',ax3,'FontSize',fontSize)
xlabel('U.V dose','Parent',ax3,'FontSize',fontSize);
ylabel('(H(u)-D(u))/(1-D(u)', 'FontSize', fontSize);
legend(get(ax3,'Children'),'Location','NW')
set(ax3,'FontSize',fontSize,'LineWidth',lineWidth)

end

if showRelativeSliding 
%____ Relative sliding contribution (h-d)/h
fig4 = figure; 
ax4  = axes('Parent',fig4,'NextPlot','add','FontSize',fontSize);
line('XData',uData,'Ydata',1 -(dData)./hData,'Marker','o','Color','k','MarkerFaceColor','c',...
    'Parent',ax4,'DisplayName','relative nucleosome sliding contribution, exp. data','MarkerSize',markerSize,...
    'LineStyle','none');
line('XData',uVals,'YData',((h(c1,c2,c3,c4,uVals))-d(c1,c2,c3,c4,uVals))./(h(c1,c2,c3,c4,uVals)),'Parent',ax4,...
    'LineWidth',lineWidth,'DisplayName','relative nucleosome sliding contribution, model')
xlabel('U.V dose','Parent',ax4,'FontSize',fontSize);
ylabel('(h(U)-d(U))/h(U)','FontSize',fontSize,'Parent',ax4)
legend(get(ax4,'Children'),'Location','NW');
title('1-D/H, Contribution of sliding to the total nucleosome loss','Parent',ax4,'FontSize',fontSize)
set(ax4,'FontSize',fontSize,'LineWidth',lineWidth)
end

if showExpansionFactor
%___ Expansion factor 
fig5 = figure;
uValues = 0.001:0.5:max(uData);
ax5  = axes('Parent',fig5,'NextPlot','add','FontSize',fontSize,'LineWidth',lineWidth);
line('XData',uValues,'YData',R(c1,c2,c3,c4,uValues),'Parent',ax5)
xlabel('U.V dose','Parent',ax5,'FontSize',fontSize)
ylabel('A(U)/A(o)','Parent',ax5,'FontSize',fontSize);
set(ax5,'FontSize',fontSize,'LineWidth',lineWidth)
end


if showRelativeOpeningDNA
fig6 = figure;
% this is the DNA loss due to opening 
% show A_openning/A(u)
uVals = 0.001:0.5:max(uData);
ax6 = axes('Parent', fig6,'NextPlot','Add','FontSize',fontSize,'LineWidth',lineWidth);
title('Relative contribution to DNA loss','FontSize',fontSize,'Parent',ax6)
xlabel(ax6,'UV dose','FontSize',fontSize)
ylabel(ax6,'Fraction of loss','FontSize',fontSize)

d_open          = ((R(c1,c2,c3,c4,uVals)-c3.*N_slide(c1,c2,uVals)-1)./(R(c1,c2,c3,c4,uVals)));%-c3.*N_slide(c1,c2,uVals)));
relativeOpening = [d_open./d(c1,c2,c3,c4,uVals);
                   1-d_open./d(c1,c2,c3,c4,uVals)];

bar(uVals,relativeOpening',1,'Stacked')
c = get(ax6,'Children');
set(c(1),'LineStyle','none','FaceColor','y');
set(c(2),'LineStyle','none')
%  legend('opening','sliding');
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
set(ax6,'XLim',[uVals(1) uVals(end)]);
end

if showRelativeOpeningHistone  
    % histone loss due to opening
  fig7  = figure;
  uVals = 0.001:0.5:max(uData);
  ax7   = axes('Parent',fig7,'NextPlot','Add');
  h_open           = ((R(c1,c2,c3,c4,uVals)-c3.*N_slide(c1,c2,uVals)-1)./(R(c1,c2,c3,c4,uVals)));%-c3.*N_slide(c1,c2,uVals)));
  relativeOpeningH = [h_open./h(c1,c2,c3,c4,uVals); 
                      1-h_open./h(c1,c2,c3,c4,uVals)];
  bar(uVals,relativeOpeningH',1.5,'Stacked')
  title('Relative contribution to nucleosome loss','FontSize',fontSize)
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
set(ax7,'Xlim',[uVals(1) uVals(end)])
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

