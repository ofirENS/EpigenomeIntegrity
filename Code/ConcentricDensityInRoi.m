function density = ConcentricDensityInRoi(particlePosition,roiRect,roiRes)
% calculate the number of particles falling within concentric band regions
% around the roi geometric center
% the density is the number of particles/(band area)
% roiRect is a 4 element row vector with the buttom left [(x,y) coordinate
% ,width, height]

rectX      = roiRect(1);
rectY      = roiRect(2);
rectWidth  = roiRect(3);
rectHeight = roiRect(4);

% Calculate the density as a function of the distance from the roi center
lx = linspace(rectX,rectX+rectWidth,2*roiRes);
ly = linspace(rectY,rectY+rectHeight,2*roiRes);
dx = diff(lx(1:2));
dy = diff(ly(1:2));

lengthIn      = zeros(1,roiRes);
rectArea      = zeros(1,roiRes);
% lengthIn(end) = dnaLengthIn;
% rectArea(end) = rectWidth*rectHeight;

for rIdx = 1:roiRes
    lengthIn(roiRes-rIdx+1) =  sum((particlePosition(:,1)<=(lx(rIdx)+rectWidth-2*(rIdx-1)*dx) & particlePosition(:,1)>=lx(rIdx) &...
    particlePosition(:,2)<=(ly(rIdx)+rectHeight-2*(rIdx-1)*dy) & particlePosition(:,2)>=ly(rIdx)));

% PolygonLengthInRoi(chainPos(:,1:2),lx(rIdx),ly(rIdx),rectWidth-2*(rIdx-1)*dx, rectHeight-2*(rIdx-1)*dy);
    rectArea(roiRes-rIdx+1)  = (rectWidth-2*(rIdx-1)*dx)*(rectHeight-2*(rIdx-1)*dy);
end
    lengthIn  = diff(lengthIn);
    bandArea  = diff(rectArea);
    % divide each density by the area of the band 
    density = (lengthIn./bandArea);