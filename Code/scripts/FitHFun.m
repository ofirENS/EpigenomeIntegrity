function fVal = FitHFun(a,u)
u = [ 5 10	15	20	25	30	35	40	45	50	55	60	65	70	75];
hData = [10.714305725	10.8220788165	14.4014983755	20.8225447327	21.2024074872...
         21.3668579387	29.195045218	37.2706560079	37.3479226024	42.5138151765...
         42.9133041668	42.8508770934	43.8660779761	42.5763929893	44.1947934168]./100;%	40.8353794651

fVal = sum((1-exp(-a(2).*(1-exp(-a(1).*u)))./(1+a(3).*(1-exp(-(a(2)).*(1-exp(-a(1).*u))))+a(4).*(1-exp(-a(1).*u)))-hData).^2);