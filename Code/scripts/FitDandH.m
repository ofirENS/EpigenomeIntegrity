function fVal = FitDandH(a,u)
u = [0 5 10	15	20	25	30	35	40	45	50	55	60	65	70	75];
hData = [0 10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872...
         21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765...
         42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168]./100;%	40.8353794651
     
dData = [0 1.5704212005	1.1365167475	4.545552178	8.7406190878	9.8581219326...
        10.2900341153	12.6333239455	20.0360763966	22.3129622161	22.5107680397...
        22.7887958612	20.4006799168	21.1679155925	22.757261652	26.9902966182	%26.4974599239
]./100;   
% T  = 1-exp(-a(1).*u);
T  = ((1-exp(-(a(1).^2)*u))+sqrt(pi.*u).*a(1).*(1-erf(a(1).*sqrt(u)))).^2;
N  = exp(-a(2).*T);
A  = 1+a(3)*(1-N)+a(4)*T;
D  = (A-1)./A;
H  = 1-N./A;
md = sum((dData-mean(dData)).^2);
mh = sum((hData-mean(hData)).^2);

% H = 1-(exp(-a(2).*(1-exp(-a(1).*u))))./(1+a(3).*(1-exp(-a(2).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)));
% D = (a(3).*(1-exp(-a(2).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)))./...
%     (1+a(3).*(1-exp(-a(2).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)));

% H = sqrt(sum((1-exp(-a(2).*(1-exp(-a(1).*u)))./(1+a(3).*(1-exp(-(a(2)).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)))-hData).^2)) ;
% D = sqrt(sum(((a(3).*(1-exp(-a(2).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)))./...
%     (1+a(3).*(1-exp(-a(2).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u))) -dData).^2));
fVal = (sum((H-hData).^2))./mh+(sum((D-dData).^2))./md;
    