% scrPlotXlsHistoneAndDnaData
% this scripts reads the histone and dna signal loss from the xls and
% compares the graphs of two analysis: one using the ensemble of data and
% the second calculates all related function using the mean of H and D
% loss.

close all 
markerSize = 15; 
fontSize   = 25;
outlierSize = 13;
boxWidths   = 0.65;
outlierSymbol = 'o';
hVals = xlsread(fullfile(pwd, '..','..','Data', 'Analyse_hoechst_021115_DataCopy.xlsx'),'Bilan','B3:Q34','basic');
dVals = xlsread(fullfile(pwd, '..','..','Data', 'Analyse_hoechst_021115_DataCopy.xlsx'),'Bilan','B37:Q68','basic');
uVals = [5:5:75, 100]; % uv dose

% Calculate mean 
mh = zeros(16,1);
md = zeros(16,1);
for i=1:16; mh(i) = mean(hVals(~isnan(hVals(:,i)),i)); md(i) = mean(dVals(~isnan(dVals(:,i)),i)); end;

% --- H ---
figure,
boxplot(hVals,uVals','boxStyle','outline','outliersize',outlierSize,'symbol',outlierSymbol,'widths',boxWidths), hold on,...
plot(mh,'og-','MarkerFaceColor','r','LineWidth',3),
plot(hVals','.','MarkerSize',markerSize),
title('H'),
xlabel('UV dose')
set(gca,'FontSize',fontSize,'LineWidth',3)
% --D --
figure, 
boxplot(dVals,uVals','boxStyle','outline','outliersize',outlierSize,'symbol',outlierSymbol ,'widths',boxWidths), hold on,...
plot(md,'og-','MarkerFaceColor','r','LineWidth',2),
plot(dVals','.','MarkerSize',markerSize), 
title('D'),
xlabel('UV dose')
set(gca,'FontSize',fontSize);
% Calculate and plot various functions of the data

% --- 1 - D/H ----
m      = 1-dVals./hVals;
m1mddh = zeros(1,16);
for i=1:16; m1mddh(i) = mean(m(~isnan(m(:,i)),i)); end;
figure, 
boxplot(1-dVals./hVals,uVals','boxStyle','outline','outliersize',outlierSize,'symbol',outlierSymbol ,'widths',boxWidths), hold on,...
plot(m1mddh,'or-','MarkerFaceColor','r','LineWidth',2),
plot(1-dVals'./hVals','.','MarkerSize',markerSize), 
plot(1-md./mh,'o-g','LineWidth',2),title('1-D/H'), 
xlabel('UV dose')
set(gca,'FontSize',fontSize)
% sub axes
axes1 = axes('Parent',gcf,...
    'Position',[0.533928571428571 0.566666666666667 0.310714285714286 0.338095238095238]);
plot(axes1, uVals,1-md./mh,'g','LineWidth',4,'DisplayName','Avg.'), hold on, 
plot(axes1,uVals,m1mddh,'r','LineWidth',4,'DisplayName','Ensamble') 
xlabel(axes1,'UV dose','FontSize',fontSize)
legend(get(axes1,'Children'))

%-- (H-D)/(100-D) ----
m          = (hVals-dVals)./(100-dVals);
mhmdd100md = zeros(1,16);
for i=1:16; mhmdd100md (i) = mean(m(~isnan(m(:,i)),i)); end;
figure, 
boxplot((hVals-dVals)./(100-dVals),uVals','boxStyle','outline','outliersize',outlierSize,'symbol',outlierSymbol,'widths',boxWidths), hold on,...
plot(mhmdd100md,'or-','MarkerFaceColor','r','LineWidth',2),
plot((hVals-dVals)'./(100-dVals)','.','MarkerSize',markerSize), 
plot((mh-md)./(100-md),'o-g','LineWidth',2),
title('(H-D)/(100-D)')
xlabel('UV dose'),
set(gca,'FontSize',fontSize);
% sub axes
axes2 = axes('Parent',gcf,...
    'Position',[0.61433868974042 0.186690647482014 0.228677379480841 0.257297019527235]);
xlabel(axes2,'UV dose','FontSize',fontSize);
plot(axes2,uVals,mhmdd100md,'r','LineWidth',4,'DisplayName','Ensamble'); hold on 
plot(axes2,uVals,(mh-md)./(100-md),'g','LineWidth',4,'DisplayName','Avg.')
legend(get(axes2,'Children'));


% -- H/D ---
m    = hVals./dVals;
mhdd = zeros(1,16);
for i=1:16; mhdd(i) = mean(m(~isnan(m(:,i)),i)); end;
figure, 
boxplot((hVals)./(dVals),uVals','boxStyle','outline','outliersize',outlierSize,'symbol',outlierSymbol ,'widths',boxWidths), hold on,...
plot(mhdd,'or-','MarkerFaceColor','r','LineWidth',2),
plot((hVals)'./(dVals)','.','MarkerSize',markerSize),
plot(mh./md,'o-g','LineWidth',2),
title('H/D')
xlabel('UV dose','FontSize',fontSize), 
set(gca,'FontSize',fontSize);
%  sub axes
axes3 = axes('Parent',gcf,...
    'Position',[0.524721878862794 0.645426515930113 0.223114956736712 0.258992805755396]);
plot(axes3,uVals,mhdd,'r','LineWidth',4,'DisplayName','Ensamble'), hold on,
plot(axes3,uVals,mh./md,'g','LineWidth',4,'DisplayName','Avg.'),
legend(get(axes3,'Children'))
xlabel(axes3,'UV dose','FontSize',fontSize)

%--- D/H ----
m    = dVals./hVals;
mddh = zeros(1,16);
for i=1:16; mddh(i) = mean(m(~isnan(m(:,i)),i)); end;
figure, 
boxplot((dVals)./(hVals),uVals','boxStyle','outline','outliersize',outlierSize,...
    'symbol',outlierSymbol,'widths',boxWidths), hold on,
plot(mddh,'or-','MarkerFaceColor','r','LineWidth',2),
plot(((dVals)./(hVals))','.'),
plot(md./mh,'o-g','LineWidth',2),
title('D/H')
xlabel('UV dose','FontSize',fontSize)
set(gca,'FontSize',fontSize);
axes4 = axes('Parent',gcf,...
    'Position',[0.454882571075402 0.214799588900308 0.252163164400494 0.327852004110997]);
plot(axes4,uVals,md./mh,'g','LineWidth',4,'DisplayName','Avg.'), hold on 
plot(axes4,uVals,mddh,'r','LineWidth',4,'DisplayName','Ensamble'),
xlabel(axes4,'UV dose')
legend(get(axes4,'Children'))

%-- H-D -- 
m = hVals-dVals;
mhmd = zeros(1,16);
for i=1:16; mhmd(i) = mean(m(~isnan(m(:,i)),i)); end;
figure, boxplot(hVals-dVals,uVals','boxStyle','outline','outliersize',outlierSize,...
    'symbol',outlierSymbol,'widths',boxWidths), hold on,

plot((hVals-dVals)','.'),
plot(mh-md,'o-g','LineWidth',2),
plot(mhmd,'or','MarkerFaceColor','r','LineWidth',2)
title('H-D')
xlabel('UV dose','FontSize',fontSize)
set(gca,'FontSize',fontSize);