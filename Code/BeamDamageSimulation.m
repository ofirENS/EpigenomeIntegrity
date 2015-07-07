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
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn         = nan(1,obj.params.numRecordingSteps+obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex       = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity  = nan(obj.params.numRecordingSteps+obj.params.numBeamSteps,obj.params.numConcentricBandsInROI-1);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentDNALoss     = nan(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentHistoneLoss = nan(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension    = nan(1,obj.params.numBeamSteps);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension = nan(1,obj.params.numBeamSteps);
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
                    obj.PreSimulationActions
                    obj.handles.framework = RouseSimulatorFramework(obj.params.simulatorParams);
                    
                    obj.PrepareResultStruct
                    obj.RunRelaxationSteps;
                    obj.InitializeGraphics;
                    obj.ShutDownDiffusion;
                    obj.RunRecordingSteps;
                    obj.RunBeamSteps;
                    obj.RunRepairSteps;
                end
            end            
        end
        
        function RunRelaxationSteps(obj)
            
            if obj.params.show3D
            % Turn the visibility of the beam off
             set(obj.handles.framework.simulationGraphics.handles.graphical.domain(obj.params.domainNumbers.beam).mesh,'Visible','off');
            end            
            obj.handles.framework.Run;% run the framework            
            obj.UpdateROIPosition;
            obj.UpdateBeamPosition;
            obj.step = obj.handles.framework.simulationData.step;
        end
        
        function RunRecordingSteps(obj)
            % Start Recording 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps;
                        
            for sIdx =1:obj.params.numRecordingSteps
                obj.handles.framework.Step;
                stepIdx = obj.handles.framework.simulationData.step;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric;
                
                obj.UpdateGraphics;                                
                obj.step = stepIdx;
            end
        end
         
        function ShutDownDiffusion(obj)
             % shut sown diffusion be
            obj.handles.framework.handles.classes.domain(obj.params.domainNumbers.sphere).params.forceParams.diffusionForce = false;
        end
        
        function PreSimulationActions(obj)
            % add connected beads 
            numBeads = obj.params.numMonomers;
            perc     = linspace(0,0.9,obj.params.numRounds);
            perc     = perc(obj.simulationRound);
            % connect perc of the monomers
            r        = randperm(numBeads);
            numPairs = round(0.5*numBeads*perc);
            if mod(numPairs,2)~=0
                numPairs = numPairs-1;
            end
            cBeads = [];
            for nIdx = 1:numPairs
               cBeads(nIdx,1) = r(2*nIdx-1);
               cBeads(nIdx,2) = r(2*nIdx);
            end
                
            obj.params.simulatorParams.chain.connectedBeads = cBeads;
        end
        
        function RunBeamSteps(obj)
            
            % increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps;  
            obj.UpdateBeamPosition;
            obj.UpdateBeamGraphics;

            
            % activate the UVC beam
            obj.ApplyDamageEffect; 
            obj.HeighlightAffectedBeads
            
            
            for sIdx =1:obj.params.numBeamSteps
                obj.handles.framework.Step
                stepIdx  = obj.handles.framework.simulationData.step;
                chainPos = obj.GetChainPosition;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(chainPos,...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                obj.UpdateGraphics
                obj.step = stepIdx;
                
            end
        end 
        
        function RunRepairSteps(obj)
                        % increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps+obj.params.numRepairSteps;  
%             obj.UpdateBeamPosition;
%             obj.UpdateBeamGraphics;

            
            % activate the UVC beam
            obj.RepairDamageEffect; 
            obj.TurnOffAffectedBeadsGraphics
            
            
            for sIdx =1:obj.params.numRepairSteps
                obj.handles.framework.Step
                stepIdx  = obj.handles.framework.simulationData.step;
                chainPos = obj.GetChainPosition;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(chainPos,...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                obj.UpdateGraphics
                obj.step = stepIdx;
                
            end
        end
        
        function [inBeam,inBeamInds] = FindMonomersInBeam(obj)
            % find the monomers in the beam 
            % Create DNA damages in the ray
            % monomers affected are those inside the 2D section of the beam
            % after projection. inBeam is a logical array indicating which monomers falls in the beam section.
            % inBEamInds are the indices of those monomers. 
            % Affected beads' nearest-neighbors are 
            % assigned bending force, therefore the collective index list
            % is given in inBeamNN
            chainPos    = obj.GetChainPosition; 
            inBeam      = obj.handles.framework.handles.classes.domain.InDomain(chainPos,obj.params.domainNumbers.beam);
            inBeamInds  = find(inBeam);            
            
    
        end
        
        function [inROI,inRoiInds,numMonomersIn,monomersInConcentric] = GetMonomerDensityInROI(obj)
                % Find the number of monomers in the rectangular ROI                
                obj.UpdateROIPosition;
                [inROI, inRoiInds,numMonomersIn] = obj.FindMonomersInROI;
                [monomersInConcentric]           = obj.GetMonomersInRoiConcentric;                              
        end
        
        function cm = GetPolymerCenterOfMass(obj)
            % polymer center of mass in any dimension 
            chainPos = obj.handles.framework.objectManager.curPos;
            cm       = mean(chainPos);
        end
        
        function UpdateROIPosition(obj)
            % Update the square roi position such that it geometrical center is at the polymer center of mass
            polymerCenterOfMass = obj.GetPolymerCenterOfMass;
            obj.roiPosition = [polymerCenterOfMass(1)-obj.params.roiWidth/2,...
                               polymerCenterOfMass(2)-obj.params.roiHeight/2,...
                               obj.params.roiWidth,...
                               obj.params.roiHeight];
        end                        
        
        function [inROI, inRoiInds,numMonomersIn]= FindMonomersInROI(obj)
            % get the number and indices of the monomers inside the ROI
            chainPosition = obj.GetChainPosition;
            inROI = chainPosition(:,1)<(obj.roiPosition(1)+obj.roiPosition(3)) &...
                    chainPosition(:,1)>(obj.roiPosition(1)) &...
                    chainPosition(:,2)<(obj.roiPosition(2)+obj.roiPosition(4)) &...
                    chainPosition(:,2)>(obj.roiPosition(2));
                
                inRoiInds     = find(inROI);
                numMonomersIn = sum(inROI);
        end
        
        function chainPosition = GetChainPosition(obj)
            % Get chain position in the last step 
            chainPosition= obj.handles.framework.objectManager.curPos;
        end
        
        function [monomersInConcentric] = GetMonomersInRoiConcentric(obj)
             % Calculate the density as a function of the distance from the roi center
             chainPosition = obj.GetChainPosition;
              monomersInConcentric = ConcentricDensityInRoi(chainPosition,obj.roiPosition,obj.params.numConcentricBandsInROI);
%               densityInConcentric = densityInConcentric./baseLineDensity;
        end
        
        function ApplyDamageEffect(obj)
            % Apply damage effect to monomers and their nearest neighbors    
            chainPos    = obj.GetChainPosition;
            obj.UpdateBeamPosition 
            [inBeam, ~] = obj.FindMonomersInBeam;% get indices and indicators for affected monomers and their nearest-neighbors
            
            % calculate the distance of each monomer in the beam to the
            % beam's center 
            r          = pdist2(chainPos,obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter);
            % calculate the probability of each monomer to be affected
            damageProb = exp(-obj.params.beamDamageSlope.*(r-obj.params.beamDamagePeak).^2);
            % use theresholding to define the actual monomers affected 
            inBeam     = (damageProb>rand(size(r,1),1) & inBeam);%obj.params.beamDamageProbThresh);
            inBeamInds = find(inBeam);
            
            % Assign bending to neighbors of affected beads
            inBeamNN = inBeam;
            for bnIdx=1:size(inBeam,1)-1
                if inBeam(bnIdx)
                    inBeamNN(bnIdx+1) = true;
                end
            end

            inBeamAffected = find(inBeamNN);
            % save 
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex = inBeamInds;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam       = inBeam;
            
            % Activate bending for affected monomers
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = true;
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = inBeamAffected;            
            
        end
        
        function RepairDamageEffect(obj)
            % Turn-off bending for affected monomers         
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = false;
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = [];
        end
        
        function UpdateBeamPosition(obj)
            cm = obj.GetPolymerCenterOfMass;
            obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter = cm;
        end                                                                 
        
        function InitializeGraphics(obj)
            %     mAxes = r.simulationGraphics.handles.graphical.mainAxes;
            
            
            % create figure for the projection in the x-y plane
            
            [rectX,rectY,rectWidth, rectHeight] = obj.GetROI;
            
            if obj.params.show3D
            % insert the projection to the 3d axes
            mainAxes = obj.handles.framework.simulationGraphics.handles.graphical.mainAxes;
            obj.handles.projPlane3D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                [rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
                'r', 'Parent',mainAxes, 'FaceAlpha',0.5);
            obj.handles.affectedBeads3D = line('XDAta',NaN,'YData',NaN,'ZData',NaN,...
                    'Marker','o',...
                    'MarkerFaceColor','r',...
                    'MarkerEdgeColor','r',...
                    'MarkerSize',7,...
                    'Parent',mainAxes,...
                    'LineStyle','none');
            end
                        
            if obj.params.show2D
                obj.handles.projectionFigure = figure;
                obj.handles.projectionAxes   = axes('Parent',obj.handles.projectionFigure,...
                    'FontSize',30);
                daspect(obj.handles.projectionAxes ,[1 1 1])
                                
                % projected Polymer
                initialChainPosition = obj.GetChainPosition;
                obj.handles.projectedPolyHandle  = line('XData',initialChainPosition(:,1),...
                    'YDAta',initialChainPosition(:,2),...
                    'Marker','o',...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor','b',...
                    'markerSize',7,...
                    'Parent',obj.handles.projectionAxes,...
                    'LineStyle','-');
                
                obj.handles.affectedBeads2D = line('XDAta',NaN,'YData',NaN,'Marker','o',...
                    'MarkerFaceColor','r',...
                    'MarkerEdgeColor','r',...
                    'MarkerSize',7,...
                    'Parent',obj.handles.projectionAxes,...
                    'LineStyle','none');
                
                % insert the projection to the projection axes
                obj.handles.projPlane2D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
                    'r', 'Parent',obj.handles.projectionAxes, 'FaceAlpha',0.5);
            end
                        
            % create density in ROI axes
            if obj.params.showDensity
                obj.handles.densityFigure = figure;
                
                obj.handles.densityAxes = axes('Parent',obj.handles.densityFigure,'NextPlot','add',...
                    'Color','none',...
                    'FontSize',30,...
                    'YLim',[0 size(initialChainPosition,1)]);
                title(obj.handles.densityAxes ,'Number of monomers in the ROI');
                xlabel(obj.handles.densityAxes,'Time [sec]','FontSize',40);
                ylabel(obj.handles.densityAxes,'Density','FontSize',40)
                obj.handles.numMonomersIn = line('XData',0,'YData',NaN,'Parent',obj.handles.densityAxes ,'Color','r',...
                    'Linewidth',4,'DisplayName','num. Beads in ROI');
                legend(obj.handles.densityAxes,get(obj.handles.densityAxes ,'Children'))
            end
            
            % show concentric density in ROI
            if obj.params.showConcentricDensity
                obj.handles.concentricDensityFigure = figure;
                obj.handles.concentricDensityAxes = axes('Parent',obj.handles.concentricDensityFigure,'NextPlot','add',...
                    'Color','none',...
                    'FontSize',30);
                title(obj.handles.concentricDensityAxes,' Concentric density around UVC center')    
                xlabel(obj.handles.concentricDensityAxes,'Dist. from CM','FontSize',40);
                ylabel(obj.handles.concentricDensityAxes,'radius','FontSize',40);
                obj.handles.concentricDensity  = line('XData',0,'YData',NaN,'Parent',obj.handles.concentricDensityAxes ,'Color','b',...
                    'LineWidth',4,'DisplayName','Density ');
                legend(obj.handles.concentricDensityAxes,get(obj.handles.concentricDensityAxes ,'Children'))
            end
            
            % show MSD of expansion after UVC
            if obj.params.showExpensionMSD
                obj.handles.expensionFigure = figure;
                obj.handles.expensionAxes = axes('Parent',obj.handles.expensionFigure,'NextPlot','add',...
                    'Color','none',...
                    'FontSize',30);
                 title(obj.handles.expensionAxes,'affected monomers expension MSD')    
                xlabel(obj.handles.expensionAxes,'Time [sec]','FontSize',40);
                ylabel(obj.handles.expensionAxes,'MSD from CM','FontSize',40)
                obj.handles.msdAffected  = line('XData',0,'YData',NaN,'Parent',obj.handles.expensionAxes,'Color','b',...
                    'LineWidth',4,'DisplayName','MSD Affected');
                obj.handles.msdNonAffected = line('XData',0,'YData',NaN,'Parent',obj.handles.expensionAxes,'Color','r',...
                    'Linewidth',4,'DisplayName','MSD Non Affected');
                legend(obj.handles.expensionAxes,get(obj.handles.expensionAxes,'Children'))
            end
                        
        end

        function [rectX,rectY,rectWidth, rectHeight] = GetROI(obj)
            % get roi x, y, width and height
            rectX      = obj.roiPosition(1);
            rectY      = obj.roiPosition(2);
            rectWidth  = obj.roiPosition(3);
            rectHeight = obj.roiPosition(4);
        end
        
        function UpdateGraphics(obj)
            % update graphics
            obj.UpdateROIGraphics            
            obj.UpdateDensityGraphics;            
            obj.UpdateProjectedPolymerGraphics
            obj.UpdateConcentricDensityGraphics
            obj.UpdateExpensionGraphics
            obj.HeighlightAffectedBeads
            drawnow
        end
        
        function UpdateProjectedPolymerGraphics(obj)
            % Update the position of the projected polymer in 2D 
            if obj.params.show2D
                chainPos = obj.GetChainPosition;
                set(obj.handles.projectedPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2))
                bIn = obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex;
                set(obj.handles.affectedBeads2D,'XData',chainPos(bIn,1), 'YDAta',chainPos(bIn,2));
            end
        end
        
        function UpdateDensityGraphics(obj)
            % Show monomer density in the ROI
            if obj.params.showDensity                
                dnaDensityDataY     = obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn;%get(obj.handles.dnaDensity,'YData');            
                dnaDensityDataX     = 1:numel(dnaDensityDataY);
                set(obj.handles.numMonomersIn,'XData',dnaDensityDataX, 'YData',dnaDensityDataY);
            end
        end
                                   
        function UpdateROIGraphics(obj)
            % projPlane3d,projPlane2D are graphical handles
            if obj.params.show2D
               [rectX,rectY,rectWidth, rectHeight] = obj.GetROI;

                if obj.params.show3D
                set(obj.handles.projPlane3D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                    'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
                end
                if obj.params.show2D
                set(obj.handles.projPlane2D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                    'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
                end
            end
        end
        
        function HeighlightAffectedBeads(obj)
            % get affected beads numbers 
            bIn      = obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex;
            chainPos = obj.GetChainPosition;
            if obj.params.show2D                
                set(obj.handles.affectedBeads2D,'XData',chainPos(bIn,1), 'YDAta',chainPos(bIn,2));
            
            end
            
            if obj.params.show3D
              set(obj.handles.affectedBeads3D,'XData',chainPos(bIn,1),'YData',chainPos(bIn,2),...
                  'ZData',chainPos(bIn,3))
            end
            
        end
        
        function TurnOffAffectedBeadsGraphics(obj)
             % Turn-off affected beasd highlighting

            if obj.params.show2D                
                set(obj.handles.affectedBeads2D,'XData',NaN, 'YDAta',NaN);
            
            end
            
            if obj.params.show3D
              set(obj.handles.affectedBeads3D,'XData',NaN,'YData',NaN,...
                  'ZData',NaN)
            end
        end
        
        function UpdateBeamGraphics(obj)%TODO: show beam in 2D
            % update the position of the 3D beam 
            if obj.params.show3D
                set(obj.handles.framework.simulationGraphics.handles.graphical.domain(obj.params.domainNumbers.beam).mesh,'Visible','on');
                domainCenter = obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter;            
                beamHandle   = obj.handles.framework.simulationGraphics.handles.graphical.domain(obj.params.domainNumbers.beam).mesh;
                xPos         = get(beamHandle,'XData');
                yPos         = get(beamHandle,'YData');
                set(beamHandle,'XData',xPos+domainCenter(1),...
                    'YData',yPos+domainCenter(2));
            end
            
            if obj.params.show2D
                
            end
            
        end        
        
        function UpdateConcentricDensityGraphics(obj)
            % show concentric monomer density in the ROI
            if obj.params.showConcentricDensity
                concentric = obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(obj.step-obj.params.numRelaxationSteps+1,:);
                set(obj.handles.concentricDensity,'XData',1:obj.params.numConcentricBandsInROI-1,'YData',concentric);
            end
        end
        
        function UpdateExpensionGraphics(obj)
            % Show expension of affected and non-affected monomers
            if obj.params.showExpensionMSD
                msdAffected    = obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension;
                msdNonAffected = obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension;
                set(obj.handles.msdAffected,'XData',1:numel(msdAffected) ,'YDAta',msdAffected);
                set(obj.handles.msdNonAffected,'XData',1:numel(msdNonAffected),'YData',msdNonAffected);
            end
            
        end
        %====================
        
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
    
    methods (Static)
        function [affectedMSDcm,nonAffectedMSDcm]= CalculateBeadsRadiusOfExpension(chainPos,inBeam)
            % calculate the mean radius of expension for affected beads and non
            % affected beads. chainPos is an numBedsXdimension vector of position
            % inBeam is a numBeadsX1 binary vector indicating the beads affected
            % (true) and non affected (false)
            cmInBeam         = mean(chainPos(inBeam,:));% center of mass for particles in beam
            cmOutBeam        = mean(chainPos(~inBeam,:)); % center of mass for paarticles out of the beam
            affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), cmInBeam).^2,2)));
            nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), cmOutBeam).^2,2)));
        end
        
    end
end