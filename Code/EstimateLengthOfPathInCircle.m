% function EstimateLengthOfPathInCircle
% calculate the number of path points contained within the inner circle of
% two concentric circles. 
% On the outer circle we choose numPairs pair of points and pass a Brownian
% bridge between them. the number of points falling within the inner circle
% is shown in the histogram 
close all 
dimension    = 3;
numPairs     = 150;
numPoints    = 2*numPairs;
numExperiments = 1;

dp(1)  = DomainHandlerParams('domainWidth',4,'dimension',dimension,'domainCenter',[0 0 0]);
dp(2)  = DomainHandlerParams('domainWidth',2,'dimension',dimension,'domainCenter',[0 0 0]);
domain = DomainHandler(dp);
bb     = BrownianBridge('realizations',1,'noiseSTD',sqrt(2*0.1*0.1),'numPoints',400,'dimension',dimension);
f = figure;
a = axes('Parent',f);hold on


inDomain = zeros(numExperiments,numPairs);
for expIdx =1:numExperiments
% randomly choose pair of points on the outer circle 
% choose pairs randomly
sampledPoints = domain.GetRandomBoundarySample(numPoints,1);
pairInds      = randperm(numPoints);

bbPath   = cell(numPairs,1); 

for pIdx = 1:numPairs
    % pass a brownianBridge between the pair of points 
    bb.params.startPoint = sampledPoints(pairInds(2*pIdx-1),:);
    bb.params.endPoint  = sampledPoints(pairInds(2*pIdx),:);
    bb.GetBridge;
    bbPath(pIdx) = bb.paths;
    % for each path, check the number of points in the inner circle 
    inDomain(expIdx,pIdx) = nnz(domain.InDomain(bbPath{pIdx},2));
end
end

for pIdx = 1:numPairs
line('XData',bbPath{pIdx}(:,1),'YData',bbPath{pIdx}(:,2),'ZData',bbPath{pIdx}(:,3),'Parent',a)
line('XData',bbPath{pIdx}(1,1),'YData',bbPath{pIdx}(1,2),'ZData',bbPath{pIdx}(1,3),'Parent',a,'Marker','o','MarkerFaceColor','r')
line('XData',bbPath{pIdx}(end,1),'YData',bbPath{pIdx}(end,2),'ZData',bbPath{pIdx}(end,3),'Parent',a,'Marker','o','MarkerFaceColor','g')
end

[sx,sy,sz] = sphere(20);
sx1 = sx*dp(1).domainWidth;
sy1 = sy*dp(1).domainWidth;
sz1 = sz*dp(1).domainWidth;
sx2 = sx*dp(2).domainWidth;
sy2 = sy*dp(2).domainWidth;
sz2 = sz*dp(2).domainWidth;

mesh(sx1,sy1,sz1,'FaceAlpha',0.3,'Parent',a);
mesh(sx2,sy2,sz2,'FaceAlpha',0.5,'Parent',a);
daspect(a,[1 1 1]), cameratoolbar
inDomain = inDomain(:); 
inDomain = inDomain(inDomain~=0);
figure,
hist(inDomain,20)
% end