function EstimateLengthOfPathInCircle
% calculate the number of path points contained within the inner circle of
% two concentric circles. 
% On the outer circle we choose numPairs pair of points and pass a Brownian
% bridge between them. the number of points falling within the inner circle
% is shown in the histogram 
dimension    = 3;
numPairs     = 100;
numPoints    = 2*numPairs;
numExperiments = 10;

dp(1) = DomainHandlerParams('domainWidth',5,'dimension',dimension);
dp(2) = DomainHandlerParams('domainWidth',2,'dimension',dimension);
domain = DomainHandler(dp);
bb     = BrownianBridge('realizations',1,'noiseSTD',0.5,'numPoints',350);

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

inDomain = inDomain(:); 
inDomain = inDomain(inDomain~=0);
hist(inDomain,20)
end