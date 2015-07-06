classdef BeamDamageSimulation<handle %[UNFINISHED]
    % This test function moves histones on a Rouse chain
    
    % % create chain and domain and register them in the ObjectManager
    
    % -Relaxation time for the Rouse chain defined by the longest relaxation time--
    % relaxation times - 300 beads ~= 1000 steps
    %                    400 beads ~= 2000 steps
    %                    500 beads ~= 3000 steps
    
    % (numBeads*b)^2 / (3*D*pi^2) %[Doi &  Edwards p.96 eq. 4.37)
    % using: b  = sqrt(3)~=1.7
    %        dt = 0.01
    %        D  = 1; diffusion const.
    % (500*sqrt(3))^2 /(3*pi^2 * 1)
    properties
        handles
        params
        results 
        roiPosition % dynamically updated [x,y,width,height]
        simulationRound = 0;
        simulation      = 0;
        step            = 0;
    end
    
    methods
                        
        function obj = BeamDamageSimulation(params)% define loadFullconfiguration
            obj.params = params;
            % set result structure             
            obj.results           = MoveHistonesOnChainResultStruct(obj.params.numRounds,obj.params.numSimulationsPerRound);
            
        end
                                
        function PrepareResultStruct(obj)
            % Preallocate the arrays in result struct for the current
            % simulation 
            cl = clock;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).date               = sprintf('%s',[num2str(cl(3)),'/',num2str(cl(2)),'/',num2str(cl(1))]);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).round              = obj.simulationRound;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).simulation         = obj.simulation;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).dimension          = obj.params.dimension;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numRelaxationSteps = obj.params.numRelaxationSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numRecordingSteps  = obj.params.numRecordingSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeamSteps       = obj.params.numBeamSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeads           = obj.params.numMonomers;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).bendingConst       = obj.params.bendingConst;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).springConst        = obj.params.springConst;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).openingAngle       = obj.params.bendingOpeningAngle;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).dt                 = obj.params.dt;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).ROI                = obj.roiPosition;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).params             = obj.params;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn         = zeros(1,obj.params.numRecordingSteps+obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex       = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity  = zeros(obj.params.numRecordingSteps+obj.params.numBeamSteps,obj.params.numConcentricBandsInROI-1);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentDNALoss     = zeros(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentHistoneLoss = zeros(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension    = zeros(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension = zeros(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPos = [];
                                    
        end
        
        function Run(obj)
            obj.simulationRound =0;            
            for rIdx = 1:obj.params.numRounds
                obj.simulation = 0;
                obj.simulationRound = obj.simulationRound+1;% increase counter
                for sIdx =1:obj.params.numSimulationsPerRound
                    obj.simulation        = obj.simulation+1;% increase counter
                    % initialize framework
                    obj.handles.framework = RouseSimulatorFramework(obj.params.simulatorParams);
                    obj.PrepareResultStruct
                    obj.RunRelaxationSteps;
                    obj.RunRecordingSteps;
                    obj.RunBeamSteps;                    
                end
            end            
        end
        
        function RunRelaxationSteps(obj)
            
            if obj.params.showSimulation
            % Turn the visibility of the beam off
             set(obj.handles.framework.simulationGraphics.handles.graphical.domain(2).mesh,'Visible','off');
            end
            
            obj.handles.framework.Run;% run the framework
            obj.step = obj.params.numRelaxationSteps;
        end
        
        function RunRecordingSteps(obj)
            % Start Recording 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps;
             stepIdx = 1;
            for sIdx =obj.handles.framework.simulationData.step:obj.handles.framework.params.simulator.numSteps
                obj.handles.framework.Step;
                [inROI,inRoiInds,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric;
                stepIdx = stepIdx+1;
                
            end
        end
                
        function RunBeamSteps(obj)
            
            % increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps;  
            if obj.params.showSimulation                                              
            % Turn on beam visibility 
              set(obj.handles.framework.simulationGraphics.handles.graphical.domain(2).mesh,'Visible','on');
            end
            % acivate the UVC beam
            obj.ApplyDamageEffect; 
            stepIdx = obj.params.numRecordingSteps+1;
            for sIdx =obj.handles.framework.simulationData.step:obj.handles.framework.params.simulator.numSteps
                obj.handles.framework.Step
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric;                
                [affectedMSDcm,nonAffectedMSDcm] = CalculateBeadsRadiusOfExpension(obj,obj.handles.framework.objectManager.curPos,...
                                                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                stepIdx = stepIdx+1;
            end
        end 
        
        function [inBeam,inBeamInds,inBeamNN] = FindMonomersInBeam(obj)
            % find the monomers in the beam 
            % Create DNA damages in the ray
            % monomers affected are those inside the 2D section of the beam
            % after projection. inBeam is a logical array indicating which monomers falls in the beam section.
            % inBEamInds are the indices of those monomers. 
            % Affected beads' nearest-neighbors are 
            % assigned bending force, therefore the collective index list
            % is given in inBeamNN
            
            chainPos    = obj.handles.framework.objectManager.curPos;
            obj.UpdateBeamPosition(chainPos) 
            inBeam      = obj.handles.framework.handles.classes.domain.InDomain(chainPos,2);
            inBeamInds  = find(inBeam);
%             numInBeam0  = sum(inBeam);% number of beads at time 0

            % Choose a fraction of inBeam to exclude
            frInBeam      = randperm(sum(inBeam));       
            inBeam(inBeamInds(frInBeam(1:round(numel(frInBeam)*(1-obj.params.fracOfMonomersAffected)))))= false;

            % Assign bending to neighbors of affected beads
            inBeamNN = inBeam;
            for bnIdx=1:size(inBeam,1)-1
                if inBeam(bnIdx)
                    inBeamNN(bnIdx+1) = true;
                end
            end

            inBeamNN = find(inBeamNN);
    
        end
        
        function [inROI,inRoiInds,numMonomersIn,monomersInConcentric] = GetMonomerDensityInROI(obj)
                % Find the number of monomers in the rectangular ROI
                chainPosition = obj.handles.framework.objectManager.curPos;
                centerOfMass  = obj.GetPolymerCenterOfMass(chainPosition);
                obj.UpdateROIPosition(centerOfMass);
                [inROI, inRoiInds,numMonomersIn] = obj.FindMonomersInROI(chainPosition);
                [monomersInConcentric] = obj.GetMonomersInRoiConcentric(chainPosition);
                              
        end
        
        function cm = GetPolymerCenterOfMass(obj,chainPos)
            % polymer center of mass in any dimension 
            cm = mean(chainPos);
        end
        
        function UpdateROIPosition(obj,polymerCenterOfMass)
            % Update the square roi position around the polymer center off mass
            obj.roiPosition = [polymerCenterOfMass(1)-obj.params.roiRadius,...
                               polymerCenterOfMass(2)-obj.params.roiRadius,...
                               2*obj.params.roiRadius,...
                               2*obj.params.roiRadius];
        end                        
        
        function [inROI, inRoiInds,numMonomersIn]= FindMonomersInROI(obj,chainPosition)
            % get the number and indices of the monomers inside the ROI
            inROI = chainPosition(:,1)<(obj.roiPosition(1)+obj.roiPosition(3)) &...
                    chainPosition(:,1)>(obj.roiPosition(1)) &...
                    chainPosition(:,2)<(obj.roiPosition(2)+obj.roiPosition(4)) &...
                    chainPosition(:,2)>(obj.roiPosition(2));
                
                inRoiInds     = find(inROI);
                numMonomersIn = sum(inROI);
        end
        
        function [monomersInConcentric] = GetMonomersInRoiConcentric(obj,chainPosition)
             % Calculate the density as a function of the distance from the roi center
              monomersInConcentric = ConcentricDensityInRoi(chainPosition,obj.roiPosition,obj.params.numConcentricBandsInROI);
%               densityInConcentric = densityInConcentric./baseLineDensity;
        end
        
        function ApplyDamageEffect(obj)
            % Apply damage effect to monomers and their nearest neighbors            
            [inBeam, inBeamInds,inBeamAffected] = obj.FindMonomersInBeam;% get indices and indicators for affected monomers and their nearest-neighbors
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex = inBeamInds;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam       = inBeam;
            obj.handles.framework.objectManager.handles.chain(1).params.forceParams.bendingElasticityForce   = true;
            obj.handles.framework.objectManager.handles.chain(1).params.forceParams.bendingAffectedParticles = inBeamAffected;
            obj.results
            
        end
        
        function UpdateBeamPosition(obj,chainPos)
            cm = mean(chainPos);
            obj.handles.framework.handles.classes.domain.params(2).domainCenter = cm;
        end                                
        
        function [affectedMSDcm,nonAffectedMSDcm]= CalculateBeadsRadiusOfExpension(obj,chainPos,inBeam)
            % calculate the mean radius of expension for affected beads and non
            % affected beads. chainPos is an numBedsXdimension vector of position
            % inBeam is a numBeadsX1 binary vector indicating the beads affected
            % (true) and non affected (false)
            cmInBeam         = mean(chainPos(inBeam,:));% center of mass for particles in beam
            cmOutBeam        = mean(chainPos(~inBeam,:)); % center of mass for paarticles out of the beam
            affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), cmInBeam).^2,2)));
            nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), cmOutBeam).^2,2)));
        end
        
        function RecordResults(obj,numMonomersIn,monomersInConcentric)
            obj.results(obj.simulationRound,obj.simulation).resultStruct
        end
                                        
        
        %================                        

                                
                
        function [histHandle,pAxes,pHistHandle,pPolyHandle,dAxes,projPlane2D,projPlane3D,dnaDensityHandle,numBeadsHandle,histoneDensityHandle]= ...
                InitializeGraphics(obj,rectX,rectY,rectWidth,rectHeight,mAxes,histonePosition,initialChainPosition)
            %     mAxes = r.simulationGraphics.handles.graphical.mainAxes;
            
            % initialize histone graphics %TODO: incorporate histone graphics in the simulationGraphics class
            
            obj.handles.histHandle = line('XData',histonePosition(:,1),...
                'YData',histonePosition(:,2),...
                'ZData',histonePosition(:,3),...
                'marker','o',...
                'MarkerFaceColor','y',...
                'MarkerSize',10,...
                'Parent',obj.handles.mAxes,...
                'LineStyle','none');
            
            % create figure for the projection in the x-y plane
            pFigure = figure;
            pAxes   = subplot(1,2,1);
            set(axes,'Parent',pFigure,...
                'XLim',get(mAxes,'XLim'), ...
                'YLim',get(mAxes,'YLim'),...
                'FontSize',30);
            daspect(pAxes,[1 1 1])
            
            % projected histone
            pHistHandle  = line('XData',histonePosition(:,1),...
                'YData',histonePosition(:,2),...
                'Marker','o',...
                'MarkerFaceColor','y',...
                'MarkerSize',10,...
                'Parent',pAxes,...
                'LineStyle','none');
            % projected Polymer
            pPolyHandle  = line('XData',initialChainPosition(:,1),...
                'YDAta',initialChainPosition(:,2),...
                'Marker','o',...
                'MarkerFaceColor','none',...
                'MarkerEdgeColor','b',...
                'markerSize',7,...
                'Parent',pAxes,...
                'LineStyle','-');
            
            % create density axes
            % dFig   = figure;
            dAxes = subplot(1,2,2);
            set(dAxes,'Parent',pFigure,'NextPlot','add','Color','none','FontSize',30,'YLim',[0 size(initialChainPosition,1)]);
            xlabel(dAxes,'Time [sec]','FontSize',40);
            ylabel(dAxes,'Density','FontSize',40)
            
            
            % insert the projection to the 3d axes
            projPlane3D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                [rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
                'r', 'Parent',mAxes, 'FaceAlpha',0.5);
            % insert the projection to the projection axes
            projPlane2D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
                'r', 'Parent',pAxes, 'FaceAlpha',0.5);
            dnaDensityHandle     = line('XData',0,'YData',NaN,'Parent',dAxes,'Color','b','LineWidth',4,'DisplayName','DNA length in ROI');
            histoneDensityHandle = line('XData',0,'YData',NaN,'Parent',dAxes,'Color','y','Linewidth',4,'DisplayName','HistoneDensity');
            numBeadsHandle       = line('XData',0,'YData',NaN,'Parent',dAxes,'Color','r','Linewidth',4,'DisplayName','num. Beads in ROI');
            legend(dAxes,get(dAxes,'Children'))
            
            
        end
        
        function [rectX,rectY]= UpdateProjectionPlanePositionByCM(chainPos,rectWidth, rectHeight)
            % set x and y to be at the cm
            
            rectX  = mean(chainPos(:,1));
            rectY  = mean(chainPos(:,2));
            rectX  = rectX- rectWidth/2;
            rectY  = rectY- rectHeight/2;
        end
        
        function UpdateGraphics(step,dt,histPosition,chainPos,histoneDensity,dnaDensity,numBeadsIn,rectX,rectY,rectWidth, rectHeight,...
                dnaDensityHandle,histoneDensityHandle,numBeadsHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
            
            UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,numBeadsIn,dnaDensityHandle,histoneDensityHandle,numBeadsHandle)
            UpdateHistoneGraphics(histHandle,histPosition)
            UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
            UpdateProjectedHistoneGraphics(pHistHandle,histPosition);
            UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
            drawnow
        end
        
        function UpdateHistoneGraphics(histHandle,curPos)
            set(histHandle,'XData',curPos(:,1),'YData',curPos(:,2),'ZData',curPos(:,3));
        end
        
        function UpdateProjectedHistoneGraphics(pHistHandle,curPos)
            set(pHistHandle,'XData',curPos(:,1),'YData',curPos(:,2));
        end
        
        function UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
            set(pPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2))
        end
        
        function UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,numBeadsIn,dnaDensityHandle,histoneDensityHandle,numBeadsHandle)
            dnaDensityDataX     = get(dnaDensityHandle,'XData');
            dnaDensityDataY     = get(dnaDensityHandle,'YData');
            histoneDensityDataX = get(histoneDensityHandle,'XData');
            histoneDensityDataY = get(histoneDensityHandle,'YData');
            numBeads            = get(numBeadsHandle,'YData');
            set(dnaDensityHandle,'XData',[dnaDensityDataX, step*dt], 'YData',[dnaDensityDataY,dnaDensity]);
            set(numBeadsHandle,'XData',[dnaDensityDataX, step*dt], 'YData',[numBeads,numBeadsIn]);
            % set(histoneDensityHandle,'XData',[histoneDensityDataX, step*dt], 'YData',[histoneDensityDataY,histoneDensity]);
            %        line('XData',step*dt,'YData',histoneDensity,'Marker','o','markerFaceColor','y','Parent',dAxes)
            %        line('XData',step*dt,'YData',dnaDensity,'Marker','o','Parent',dAxes,'markerFaceColor','b')
        end
        
        function UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
            % projPlane3d,projPlane2D are graphical handles
            set(projPlane3D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
            set(projPlane2D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
        end
        
        function UpdateGyrationSphere(simulationFrameworkHandle,domainNumber, gyrationSphereData,chainPos)
            % update the sphere graphics
            ccm  = mean(chainPos,1);% chain center of mass
            xPos = gyrationSphereData.points.x;
            yPos = gyrationSphereData.points.y;
            zPos = gyrationSphereData.points.z;
            
            xPos = xPos -xPos(1,1)+ccm(1);
            yPos = yPos -yPos(1,1)+ccm(2);
            zPos = zPos -zPos((size(zPos,1)-1)/2,1)+ccm(3);
            
            set(gyrationSphereData.mesh,'XData',xPos,...
                'YData',yPos,...
                'ZData',zPos);
            simulationFrameworkHandle.handles.classes.domain.params(3).domainCenter = ccm;
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.x = xPos;
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.y = yPos;
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.z = zPos;
        end
        
        function UpdateBeamGraphics(simulationFrameworkHandle,domainNumber)
            domainCenter = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainCenter;
            %     domainRad    = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainWidth;
            beamHandle   = simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).mesh;
            xPos         = get(beamHandle,'XData');
            yPos         = get(beamHandle,'YData');
            set(beamHandle,'XData',xPos+domainCenter(1),...
                'YData',yPos+domainCenter(2));
            
        end
        
        function [histoneDensity, dnaDensity,numBeads,densityInConcentric] = CalculateDensitiesInROI(chainPos,histonePos,rectX,rectY, rectWidth, rectHeight,roiRes,numHistones,initialDnaInLength,baseLineDensity)
            % calculate histone density
            histoneDensity =  sum((histonePos(:,1)<=(rectX+rectWidth) & histonePos(:,1)>=rectX &...
                histonePos(:,2)<=(rectY+rectHeight) & histonePos(:,2)>=rectY))/numHistones;
            
            numBeads =  sum((chainPos(:,1)<=(rectX+rectWidth) & chainPos(:,1)>=rectX &...
                chainPos(:,2)<=(rectY+rectHeight) & chainPos(:,2)>=rectY));
            
            % Calculate DNA density
            % [dnaLengthIn,totalDNALength] = PolygonLengthInRoi(chainPos(:,1:2),rectX,rectY,rectWidth, rectHeight);
            dnaLengthIn = 1;
            dnaDensity = dnaLengthIn;%./initialDnaInLength;%totalDNALength;
            
            
            % Calculate the density as a function of the distance from the roi center
            densityInConcentric = ConcentricDensityInRoi(chainPos,[rectX,rectY,rectWidth,rectHeight],roiRes);
            densityInConcentric = densityInConcentric./baseLineDensity;
            
        end
        
        function LoadFullConfiguration(obj)% finish loading 
        end
        
        function LoadRelaxationConfiguration(obj)% need to overwrite params
                       
            [fName,fPath] = uigetfile('*.mat');
            load(fullfile(fPath,fName)); % simulation framework, saved as r
            f = whos;
            cInd = strcmpi({f.class},'RouseSimulatorFramework');
            if ~any(cInd)% class indicator
                error('this is not a file of class rouseSimulatorFramework')
            end
            simFramework = eval(f(cInd).name);
           
        end
        
        function SaveRelaxationConfiguration(obj,r)
            % save RouseSimulatorFramework class            
              uisave('r')            
        end% remove dependancy on variable r
        
        
    end
end