% scrTestExpansionFactor
% define the lower boundary of the expansion factor as a function of: h- histone loss,
% alpha-1/pi, gamma-chromatin relaxation factor

RTmin = @(alpha,h,gamma) (alpha+1)./(1-h+alpha.*gamma);
RTmax = @(alpha,h,gamma) 1./(1-h+alpha.*(gamma-1));
slidingMin = @(alpha,h,gamma) 1-((1-h)./alpha +gamma -(1./(RTmax(alpha,h,gamma).*alpha)));
slidingMax = @(alpha,h,gamma) 1-((1-h)./alpha +gamma -(1./(RTmin(alpha,h,gamma).*alpha)));
gamma = linspace(0,1,100);
alpha = 1/pi;
h     = linspace(0.05,0.5, 4);
fig1  = figure; ax1 = axes('Parent',fig1);hold on
fig2  = figure; ax2 = axes('Parent',fig2); hold on
fig3  = figure; ax3 = axes('Parent',fig3); hold on
for hIdx = 1:numel(h)
    % lower expansion factor bound
    k=plot(ax1,gamma,RTmin(alpha,h(hIdx),gamma),'DisplayName',['Lower bound, histone Loss fraction: ' num2str(h(hIdx))],'Linewidth',3);
    % upper expansion factor bound 
    plot(ax1,gamma,RTmax(alpha,h(hIdx),gamma),'DisplayName',['Upper bound, histone Loss fraction: ' num2str(h(hIdx))],...
        'Linewidth',3,'LineStyle','-.','Color',k.Color)
    k2 =plot(ax2,gamma,slidingMin(alpha,h(hIdx),gamma),'DisplayName',['Min Sliding contrib. histone loss fraction ',num2str(h(hIdx))],'LineWidth',3);
    plot(ax2,gamma,slidingMax(alpha,h(hIdx),gamma),'DisplayName',['Max Sliding contrib. histone loss fraction ',num2str(h(hIdx))],'LineWidth',3,...
        'LineStyle','-.','Color',k2.Color);
    % histone loss as a function of expansion factor 
    k3 = plot(ax3,RTmin(alpha,h(hIdx),gamma),slidingMax(alpha,h(hIdx),gamma),'LineWidth',3,'DisplayName',['Max sliding, min Expansion, histone loss fraction ', num2str(h(hIdx))]);
    plot(ax3,RTmax(alpha,h(hIdx),gamma),slidingMax(alpha,h(hIdx),gamma),'Color',k3.Color,'LineWidth',3,'LineStyle','-.','DisplayName',['Max sliding, max Expansion, histone loss fraction ', num2str(h(hIdx))]);
        
end
% expansio factor axes
title(ax1,'Expansion factor as a function of \gamma'), xlabel(ax1,'\gamma'), ylabel(ax1,'Expansion factor')
set(ax1,'FontSize',25,'LineWidth',3)
legend(ax1,get(ax1,'Children'))

% sliding contribution 
title(ax2,'Sliding contribution as a function of \gamma'), xlabel(ax2,'\gamma'), ylabel(ax2,'Sliding fraction')
set(ax2,'FontSize',25,'LineWidth',3)
legend(ax2,get(ax2,'Children'))
legend(ax2,get(ax2,'Children'))


title(ax3,'Sliding contribution as a function of expansion factor'), xlabel(ax3,'expansion factor'), ylabel(ax3,'Sliding fraction')
set(ax3,'FontSize',25,'LineWidth',3)
legend(ax3,get(ax3,'Children'))
legend(ax3,get(ax3,'Children'))
