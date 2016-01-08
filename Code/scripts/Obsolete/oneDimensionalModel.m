% a1 = k_t
% a2 = c_dA_0
% a3 = A_0/(k_p+k_s)
% a4 = c_s/(A_0c_n)
% assuming linear relationship between the area of sliding pushing and the
% damaged DNA 
T = @(a1,a2,u) a2.*(1-exp(-a1.*u)); % damaged DNA
D = @(a1,a2,a3,u) 1-1./(1+a3.*T(a1,a2,u)); % total DNa fraction loss
H = @(a1,a2,a3,a4,u) (1+a3.*a4).*T(a1,a2,u)./(a3+T(a1,a2,u)); % total nucleosome fraction of loss

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
% % % ---autofit
opt       = optimset('TolFun',1e-15,'TolX',1e-15,'MaxIter',1e7,'MaxFunEvals',1e7,'TolCon',1e-19,'Hessian','bfgs',...
    'Diagnostics','off');
% run several tests
numTests  = 1;
fitParams = zeros(numTests,4); 
fval      = zeros(numTests,1);
for tIdx = 1:numTests
[fitParams(tIdx,:),fval(tIdx),exitFlag,output]=...
    fmincon(@OneDimensionalModelFit,0.05*rand(1,4),-1*eye(4),zeros(4,1),[],[],zeros(4,1),2*ones(4,1),[],opt);
end
a1 = fitParams(1);
a2 = fitParams(2);
a3 = fitParams(3);
a4 = fitParams(4);

figure, plot(uData,dData,'og'), hold on, plot(uData,D(a1,a2,a3,uData),'g')
 plot(uData,hData,'or'),  plot(uData,H(a1,a2,a3,a4,uData),'r')