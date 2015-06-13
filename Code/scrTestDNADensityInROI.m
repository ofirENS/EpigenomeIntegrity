% scrTestDNADensityInROI
% this is a second version of the calculation of DNA denity in ROI. it is
% used for comparison and validation of the framework calculation, as an
% indipendent code
% parameter values are set proportional to (dimension*diffusionConst/b^2)
close all 
%% parameters
numParticles     = 800;
dimension        = 3;
% relxation time 
% (numParticles*b)^2 / (3*diffusionConst*pi^2)
numSteps         = 800;
beamTime         = 200; 
fixTime          = 600;
dt               = 0.001;
angle0           = pi;
b                = sqrt(3);
diffusionConst   = 1;
springConst      = 2*(dimension*diffusionConst./b^2)*ones(numParticles);
bendingConst     = 2*dimension*diffusionConst./b^2;
particlePosition = cumsum(sqrt(2*diffusionConst*dt)*randn(numParticles,dimension));

connectivityMap  = (diag(ones(1,numParticles-1),1)+diag(ones(1,numParticles-1),-1))~=0;
minParticleDist  = 0.5;
fixedParticleNum = [];

gyrationRadius   = sqrt(numParticles/6)*b; % for large numParticles
beamRad          = gyrationRadius/25;
roiRadius        = gyrationRadius/20;
roiRes           = 20; % number of rows and columns in roi
affectedBeads    = false(numParticles,1);

% flags
diffusionFlag    = true;
springsFlag      = true;
bendingFlag      = true;

%% Create Graphics

% Create graphics
mainFig       = figure('Units','norm');
densityFigure = figure('Units','norm');
projFig       = figure('Units','norm');  
mainAxes      = axes('Parent',mainFig,'Units','norm','Color','k','XLim',gyrationRadius.*[-1 1],...
                     'YLim',gyrationRadius.*[-1 1],'ZLim',gyrationRadius.*[-1 1],'FontSize',25);
densityAxes   = axes('Parent',densityFigure,'Units','norm','FontSize',25,'YLim',[0 1]);
projAxes      = axes('Parent',projFig,'Units','norm','FontSize',25,'XLim',get(mainAxes,'XLim'),'YLim',get(mainAxes,'YLim'));

xlabel(densityAxes,'Time'); ylabel(densityAxes,'Num. Beads In ROI')

cameratoolbar(mainFig);
daspect(mainAxes,[1 1 1]);
daspect(projAxes,[1 1 1]);

particleHandle = line('XData',particlePosition(:,1),...
                      'Ydata',particlePosition(:,2),...
                      'Zdata',particlePosition(:,3),...
                      'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
                      'LineStyle','-','Color','w','LineWidth',1,'Parent',mainAxes);                  

projBeadHandle = line('XData',particlePosition(:,1),...
                      'Ydata',particlePosition(:,2),...                      
                      'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
                      'LineStyle','-','Color','k','LineWidth',1,'Parent',projAxes);                  
                  
                  
affectedBeadsHandle = line('XData',NaN,...
                           'Ydata',NaN,...
                           'Zdata',NaN,...
                           'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g',...
                           'LineStyle','none','Color','w','LineWidth',1,'Parent',mainAxes);
affectedBeadsProjHandle   = line('XData',particlePosition(affectedBeads,1),...
                           'Ydata',particlePosition(affectedBeads,2),...                           
                           'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                           'LineStyle','none','Color','r','LineWidth',1,'Parent',projAxes,'Visible','off');
                       
numBeadsInHandle    = line('XData',NaN,'YData',NaN,'Color','b','Parent',densityAxes,'LineWidth',4); 

cm    = mean(particlePosition,1); % center of mass

% create the beam 
numCirclesInBeam = 40;
circlezPos       = linspace(-2*gyrationRadius,2*gyrationRadius,numCirclesInBeam);
theta            = linspace(0,2*pi,30);
cosTheta = cos(theta);
sinTheta = sin(theta);
for cIdx = 1:numCirclesInBeam;
    beamCircleHandle(cIdx)= line('XData',beamRad*cosTheta+cm(1),'YData',beamRad*sinTheta+cm(2),'ZData',ones(1,numel(theta))*circlezPos(cIdx),...
                                'Parent',mainAxes,'Color','y','Visible','off');
end    
        

 % projection plane in 3D
roiHandle           = patch([cm(1)-roiRadius,cm(1)+roiRadius,cm(1)+roiRadius,cm(1)-roiRadius],...
                            [cm(2)-roiRadius,cm(2)-roiRadius,cm(2)+roiRadius,cm(2)+roiRadius],'y',...
                                'Parent',mainAxes,'FaceAlpha',0.5);
projRoiHandle       = patch([cm(1)-roiRadius,cm(1)+roiRadius,cm(1)+roiRadius,cm(1)-roiRadius],...
                            [cm(2)-roiRadius,cm(2)-roiRadius,cm(2)+roiRadius,cm(2)+roiRadius],'y',...
                                'Parent',projAxes,'FaceAlpha',0.5);  
                            
baseLineDensity = ConcentricDensityInRoi(particlePosition,[cm(1)-roiRadius,cm(2)-roiRadius,2*roiRadius,2*roiRadius],roiRes);

%% Run Simulation                            
for sIdx = 1:numSteps
    
    particleDist = ForceManager.GetParticleDistance(particlePosition);
    
    % Calculate forces
    bendingForce   = ForceManager.GetBendingElasticityForce(bendingFlag,particlePosition,particleDist,connectivityMap,bendingConst,fixedParticleNum,angle0);
    bendingForce(~affectedBeads,:) = 0;% zero out bending forces for unaffected beads 
    
    diffusionForce = ForceManager.GetDiffusionForce(diffusionFlag,particlePosition,diffusionConst,dt,fixedParticleNum);
    springForce    = ForceManager.GetSpringForce(springsFlag,particlePosition,particleDist,springConst,connectivityMap,minParticleDist,fixedParticleNum);
    
    % Update position
    particlePosition = particlePosition +(bendingFlag*bendingForce +springsFlag*springForce)*dt+ diffusionFlag*diffusionForce;
    
    cm               = mean(particlePosition,1); % center of mass
    
    density = ConcentricDensityInRoi(particlePosition,[cm(1)-roiRadius,cm(2)-roiRadius,2*roiRadius,2*roiRadius],roiRes);
    
%     inROI = (particlePosition(:,1)<(cm(1)+roiRadius)) & (particlePosition(:,1)>(cm(1)-roiRadius)) & ...
%             (particlePosition(:,2)<(cm(2)+roiRadius)) & (particlePosition(:,2)>(cm(2)-roiRadius));
            
    % Update graphics
    set(particleHandle,'XData',particlePosition(:,1),...
                       'YData',particlePosition(:,2),...
                       'ZData',particlePosition(:,3));
    set(projBeadHandle,'XData',particlePosition(:,1),...
                       'YData',particlePosition(:,2));
    set(affectedBeadsProjHandle,'XData',particlePosition(affectedBeads,1),...
                                'YData',particlePosition(affectedBeads,2));
                            
    set(affectedBeadsHandle,'XData', particlePosition(affectedBeads,1),...
                            'YData',particlePosition(affectedBeads,2),...
                            'ZData',particlePosition(affectedBeads,3));
    set(numBeadsInHandle,'XData',(0:(roiRes-2)),...
                         'YData',(density./baseLineDensity(1)))
                     
    set([roiHandle,projRoiHandle],'XData',[cm(1)-roiRadius,cm(1)+roiRadius,cm(1)+roiRadius,cm(1)-roiRadius],...
                  'YData',[cm(2)-roiRadius,cm(2)-roiRadius,cm(2)+roiRadius,cm(2)+roiRadius]);
              
    for cIdx = 1:numCirclesInBeam
            set(beamCircleHandle(cIdx),'XData',beamRad*cosTheta+cm(1),'YData',beamRad*sinTheta+cm(2))
    end          
%     set(get(beamHandleGroup,'Children'),'XData',cellfun(get(get(beamHandleGroup,'Children'),'XData')-cm(1)),'YData',get(beamHandleGroup,'Children','YData')-cm(2));          
    drawnow
    
    if sIdx ==beamTime
        % Apply beam on particles
        % find all particles around center of mass in radius beamRad
        affectedBeads       = sum(bsxfun(@minus,particlePosition(:,1:2),cm(1:2)).^2,2)<=beamRad^2;
        affectedBeadsHandle = line('XData',particlePosition(affectedBeads ,1),...
                                   'Ydata',particlePosition(affectedBeads ,2),...
                                   'Zdata',particlePosition(affectedBeads ,3),...
                                   'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                                   'LineStyle','none','Color','w','LineWidth',2,'Parent',mainAxes); 
        for cIdx = 1:numCirclesInBeam
            set(beamCircleHandle(cIdx),'Visible','on')
        end
        
        set(affectedBeadsProjHandle,'Visible','on','XData',particlePosition(affectedBeads,1),...                         
                                    'YData',particlePosition(affectedBeads,2));
        
        % calculate baseLine density
        baseLineDensity = ConcentricDensityInRoi(particlePosition,[cm(1)-roiRadius,cm(2)-roiRadius,2*roiRadius,2*roiRadius],roiRes);
        
        % turn off diffusion 
        diffusionFlag   = false;
        minParticleDist = minParticleDist*2;
%         for aIdx = 1:numParticles            
%             if affectedBeads(aIdx)
%                springConst(aIdx,aIdx) = springConst(aIdx,aIdx)/2;
%             end
%         end
        sprintf('%s%f','Beam fired at time: ', sIdx*dt)
        
    end
    
    if sIdx ==fixTime
        % fixind time 
%       for aIdx = 1:numParticles            
%             if affectedBeads(aIdx)
%                springConst(aIdx,aIdx) = springConst(aIdx,aIdx)*4;
%             end
%       end
        minParticleDist = minParticleDist/2;
        affectedBeads   = false(numParticles,1);
        sprintf('%sf','Fixed at time: ', sIdx*dt); 
    end
end
%% save Results
save particlePosition particlePosition
