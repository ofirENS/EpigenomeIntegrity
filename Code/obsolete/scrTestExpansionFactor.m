% scrTestExpansionFactor
% define the lower boundary of the expansion factor as a function of: h- histone loss,
% alpha-1/pi, gamma-chromatin relaxation factor
close all 
RTmin = @(alpha,h,gamma) (1+alpha)./(1-h+alpha.*gamma);
RTmax = @(alpha,h,gamma) 1./(1-h-alpha.*(1-gamma));
slidingMin = @(alpha,h,gamma) 1-((1-h)./alpha +gamma -(1./(RTmax(alpha,h,gamma).*alpha)));
slidingMax = @(alpha,h,gamma) 1-((1-h)./alpha +gamma -(1./(RTmin(alpha,h,gamma).*alpha)));
% s =@(alpha,h,gamma) 1-(alpha./(1+alpha)).*((1-h)./alpha +gamma);
s=@(alpha,h,d,gamma) (h./(d.*alpha) -1./(alpha.*gamma) -(gamma-1)./d);
% gamma = linspace(0,1,100);
alpha = 1/pi;
h     = linspace(0.3,0.4, 5);
d     = linspace(0.23,0.23,5);
fig1  = figure; ax1 = axes('Parent',fig1);hold on
fig2  = figure; ax2 = axes('Parent',fig2); hold on
fig3  = figure; ax3 = axes('Parent',fig3); hold on
fig4 = figure; ax4 = axes('Parent',fig4); hold on
for hIdx = 1:numel(h)
    % calculate valid limits for gamma
g1 = sort(roots([1, -(1+h(hIdx)/alpha), d(hIdx)/alpha]));% valid values between roots
g2 = sort(roots([1,-(1-d(hIdx)+h(hIdx)/alpha),d(hIdx)/alpha]));% valid values smaller/bigger than roots

% truncate the range between 0 and 1
g1(1) = max([0 g1(1)]); g1(2) = min([g1(2) 1]);
% truncate the roots between 0 and 1
g2(1) = max([0 g2(1)]); g2(2) = min([g2(2),1]);
% set the valid range 
gRange = [g1(1):0.001:g2(1), g2(2):0.001:g1(2)];
gRange= gRange(1:end-1);
if max(gRange)>1
    disp('')
end
cmap = jet(numel(h));
if all(isreal(g1)) && all(isreal(g2))
% gRange = [max([0,min(g1)]):.001:max([0,min(g2)]), min([max(g1),1]):.001:min([max(g2),1])];
    % lower expansion factor bound
%     k=plot(ax1,gRange,RTmin(alpha,h(hIdx),gRange),'DisplayName',['Lower bound, histone Loss fraction: ' num2str(h(hIdx))],'Linewidth',3);
    % upper expansion factor bound 
%     plot(ax1,gRange,RTmax(alpha,h(hIdx),gRange),'DisplayName',['Upper bound, histone Loss fraction: ' num2str(h(hIdx))],...
%         'Linewidth',3,'LineStyle','-.','Color',cmap(hIdx,:))
%     k2 =plot(ax2,gRange,slidingMin(alpha,h(hIdx),gRange),'DisplayName',['Min Sliding contrib. histone loss ',num2str(h(hIdx))],'LineWidth',3);
%     plot(ax2,gRange,s(alpha,h(hIdx),gRange),'DisplayName',['Max Sliding contrib. histone loss ',num2str(h(hIdx))],'LineWidth',3,...
%         'LineStyle','-.','Color',cmap(hIdx,:));
    % histone loss as a function of expansion factor 
%     k3 = plot(ax3,s(alpha,h(hIdx),d(hIdx),gRange),RTmin(alpha,h(hIdx),gRange),'LineWidth',3,'DisplayName',['Min sliding, min Expansion, histone loss ', num2str(h(hIdx))],'Color',cmap(hIdx,:));
    k3= plot(ax3,RTmax(alpha,h(hIdx),gRange),s(alpha,h(hIdx),d(hIdx),gRange),'Color',cmap(hIdx,:),'LineWidth',3,'LineStyle','-.','DisplayName',['Max Expansion, histone loss ', num2str(h(hIdx)),', DNA loss ', num2str(d(hIdx))]);
        plot(ax3,RTmin(alpha,h(hIdx),gRange),s(alpha,h(hIdx),d(hIdx),gRange),'Color',cmap(hIdx,:),'LineWidth',3,'LineStyle','-','DisplayName',['Min Expansion, histone loss ', num2str(h(hIdx)),', DNA loss ', num2str(d(hIdx))]);
        
    k4 = plot(ax4,gRange,s(alpha,h(hIdx),d(hIdx),gRange),'-','LineWidth',3,'DisplayName', ['histone loss ', num2str(h(hIdx)),', DNA loss ', num2str(d(hIdx))],'Color',cmap(hIdx,:));
end
end
% expansion factor axes
title(ax1,'Expansion factor as a function of \gamma'), xlabel(ax1,'\gamma'), ylabel(ax1,'Expansion factor')
set(ax1,'FontSize',25,'LineWidth',3,'YLim',[ 0 1])
legend(ax1,get(ax1,'Children'))

% sliding contribution 
title(ax2,'Sliding contribution as a function of \gamma'), 
xlabel(ax2,'\gamma'), ylabel(ax2,'Sliding fraction')
set(ax2,'FontSize',25,'LineWidth',3)
legend(ax2,get(ax2,'Children'))


title(ax3,'Expansion as a function of sliding','FontSize',25), 
xlabel(ax3,'expansion factor','FontSize',25), ylabel(ax3,'Sliding fraction','FontSize',25)
set(ax3,'FontSize',25,'LineWidth',3,'YLim',[0 1])
legend(ax3,get(ax3,'Children'))

title(ax4,'Sliding contribution as a function of \gamma','FontSize',25)
    xlabel(ax4,'\gamma','FontSize',25), ylabel(ax4,'Sliding fraction','FontSize',25)
set(ax4,'FontSize',25,'LineWidth',3,'YLim',[0 0.995])
legend(ax4,get(ax4,'Children'))
