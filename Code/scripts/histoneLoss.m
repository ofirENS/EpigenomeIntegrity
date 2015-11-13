function goalFunc = histoneLoss(coef)
% calculate actually the DNA
[~,fVals] = ode45(@(u,y)radius(u,y,coef),5:5:75,1);

dnaLoss = [1.5704212005, 1.1365167475	4.545552178	8.7406190878	9.8581219326	10.2900341153...
           12.6333239455	20.0360763966	22.3129622161 22.5107680397	22.7887958612	20.4006799168...
           21.1679155925	22.757261652	26.9902966182]./100;
     
goalFunc = sum((fVals -(1./(1-dnaLoss))').^2);
 line('XData',5:5:75,'YData',(1./(1-dnaLoss)),'Marker','o','Color','r');
line('XData',5:5:75,'YData',fVals,'Marker','o','Color','g')
function dy = radius(u,y,coef)
dy = exp(-coef(1)*u).*(coef(2)-coef(3).*u);