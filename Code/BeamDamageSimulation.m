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
        date
        handles
        params
        results 
        description
        roiPosition % dynamically updated [x,y,width,height]
        beamCenterPosition
        simulationRound = 0;
        simulation      = 0;
        step            = 0;
    end
    
    methods
                        
        function obj = BeamDamageSimulation(params)% define loadFullconfiguration
            obj.params      = params;            
            cl              = clock;
            obj.date        = sprintf('%s',[num2str(cl(3)),'/',num2str(cl(2)),'/',num2str(cl(1))]);            
            obj.results     = obj.NewResultStruct;
            obj.description = params.description;
            
            for rIdx = 1:obj.params.numRounds
                for sIdx = 1:obj.params.numSimulationsPerRound
                 obj.results.resultStruct(rIdx,sIdx) = obj.NewResultStruct;                 
                end
            end
        end
                                
        function PrepareResultStruct(obj)
            % Preallocate the arrays in result struct for the current
            % simulation 
            cl = clock;
            numStepsToRecord = obj.params.numRecordingSteps+obj.params.numBeamSteps+obj.params.numRepairSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).date               = sprintf('%s',[num2str(cl(3)),'/',num2str(cl(2)),'/',num2str(cl(1))]);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).round              = obj.simulationRound;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).simulation         = obj.simulation;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).dimension          = obj.params.dimension;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numRelaxationSteps = obj.params.numRelaxationSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numRecordingSteps  = obj.params.numRecordingSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeamSteps       = obj.params.numBeamSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numRepairSteps     = obj.params.numRepairSteps;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeads           = obj.params.numMonomers;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).bendingConst       = obj.params.bendingConst;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).springConst        = obj.params.springConst;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).openingAngle       = obj.params.bendingOpeningAngle;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).dt                 = obj.params.dt;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).ROI                = obj.roiPosition;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).params             = obj.params;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn         = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex       = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity  = nan(numStepsToRecord,obj.params.numConcentricBandsInROI-1);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentDNALoss     = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentHistoneLoss = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension    = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPos = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads          = [];  
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentOfConnectedBeads = obj.params.percentOfConnectedMonomers;
            
            obj.results.resultStruct(obj.simulationRound,obj.simulation).description      = obj.params.description;
            
        end
        
        function Run(obj)
            obj.simulationRound =0;            
            for rIdx = 1:obj.params.numRounds
                obj.simulation = 0;
                obj.simulationRound = obj.simulationRound+1;% increase counter
                obj.PreRoundActions
                for sIdx =1:obj.params.numSimulationsPerRound
                    obj.simulation        = obj.simulation+1;% increase counter
                    % initialize framework
%                     obj.params           = BeamDamageParams;
                    obj.PreSimulationActions      
                    if isfield(obj.handles,'framework')
                        delete(obj.handles.framework)
                    end
                    obj.handles.framework = RouseSimulatorFramework(obj.params.simulatorParams);
                    
                    obj.PrepareResultStruct
                    obj.RunRelaxationSteps;
                    obj.InitializeGraphics;
                    obj.ShutDownDiffusion;
                    obj.RunRecordingSteps;
                    obj.RunBeamSteps;
                    obj.RunRepairSteps;
                    obj.PostSimulationActions
                end
                obj.PostRoundActions
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
            obj.step = obj.params.numRelaxationSteps+1;% handles.framework.simulationData.step;
        end
        
        function RunRecordingSteps(obj)
            % Start Recording 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps;
             
            for sIdx =1:obj.params.numRecordingSteps
                obj.handles.framework.Step;
                stepIdx  = obj.step -obj.params.numRelaxationSteps;
                chainPos = obj.GetChainPosition;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                obj.UpdateGraphics;                                
                obj.step = obj.step+1;%obj.handles.framework.simulationData.step;
            end
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
                stepIdx  = obj.step-obj.params.numRelaxationSteps;
                chainPos = obj.GetChainPosition;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(chainPos,...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                obj.UpdateGraphics
                obj.step = obj.step+1;
                
            end
        end
        
        function RunRepairSteps(obj)
                        % increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps+obj.params.numRepairSteps;  
                        
            obj.RepairDamageEffect; 
            obj.TurnOffAffectedBeadsGraphics
            
            
            for sIdx =1:obj.params.numRepairSteps
                obj.handles.framework.Step
                stepIdx  = obj.step-obj.params.numRelaxationSteps;
                chainPos = obj.GetChainPosition;
                [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(chainPos,...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                obj.UpdateGraphics
                obj.step = obj.step+1;                
            end
        end
                
        function ShutDownDiffusion(obj)
             % shut down diffusion be
             if obj.params.shutDownDiffusionAfterRelaxationSteps
               obj.handles.framework.handles.classes.domain(obj.params.domainNumbers.sphere).params.forceParams.diffusionForce = false;
%                obj.handles.framework.handles.classes.domain(obj.params.domainNumbers.sphere).params.forceParams.diffusionConst = 0.001;
             end
        end
        
        function PreRoundActions(obj)
            % actions performed before the begining of each round           
        end
        
        function PostRoundActions(obj)
            % Actions performed at the end of each round   
            
            if obj.params.saveAfterEachRound
                % save results to result folder     
                obj.SaveResults
            end
        end
        
        function PreSimulationActions(obj)
            % Add connected beads 
            if ~isempty(obj.params.tryConnectivity)
              obj.params.percentOfConnectedMonomers = obj.params.tryConnectivity(obj.simulationRound);
              obj.params.InitializeParamClasses;              
            end
            
            cl = clock;
            cl = [num2str(cl(4)),':' num2str(cl(5)), ':', num2str(cl(6))];
            sprintf('%s%s%s%s%s%s','Started Round ', num2str(obj.simulationRound), ' Simulation ', num2str(obj.simulation),' at time: ', cl)
        end
        
        function PostSimulationActions(obj)
     
            if obj.params.saveAfterEachSimulation
                % save results to result folder                 
                obj.SaveResults                
            end
            cl = clock;
            cl = [num2str(cl(4)),':' num2str(cl(5)), ':', num2str(cl(6))];
            sprintf('%s%s%s%s%s%s','Ended Round ', num2str(obj.simulationRound), ' Simulation ', num2str(obj.simulation),' at time: ', cl)  
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
        
        function [affectedMSDcm,nonAffectedMSDcm] = CalculateBeadsRadiusOfExpension(obj,chainPos,inBeam)
            % calculate the mean radius of expension for affected beads and non
            % affected beads. chainPos is an numBedsXdimension vector of position
            % inBeam is a numBeadsX1 binary vector indicating the beads affected
            % (true) and non affected (false)
         
            if obj.params.calculateMSDFromCenterOfMass
               cmInBeam         = mean(chainPos(inBeam,:));% center of mass for particles in beam
               cmOutBeam        = mean(chainPos(~inBeam,:)); % center of mass for paarticles out of the beam
               affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), cmInBeam).^2,2)));
               nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), cmOutBeam).^2,2)));
            elseif obj.params.calculateMSDFromBeamCenter
                affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), obj.beamCenterPosition).^2,2)));
                nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), obj.beamCenterPosition).^2,2)));
            else 
                error('unsupported option')
            end
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
            obj.beamCenterPosition = obj.GetPolymerCenterOfMass;
            chainPos               = obj.GetChainPosition;
            
            obj.UpdateBeamPosition 
            [inBeam, ~] = obj.FindMonomersInBeam;% get indices and indicators for affected monomers and their nearest-neighbors
            
            % Calculate the distance of each monomer in the beam to the
            % beam's center 
            r = pdist2mex(chainPos',obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter','euc',[],[],[]);
            
            % Calculate the probability of each monomer to be affected
            damageProb = exp(-obj.params.beamDamageSlope.*(r-obj.params.beamDamagePeak).^2);
            
            % use theresholding to define the actual monomers affected 
            affectedInBeam     = (damageProb>rand(size(r,1),1) & inBeam);
            inBeamInds = find(affectedInBeam);
                        
            inBeamNN = affectedInBeam;
            
            if obj.params.assignBendingToAffectedMonomers
             % Assign bending also to neighbors of damaged monomers, to
             % create continuity
            for bnIdx=1:size(affectedInBeam,1)-1
                if affectedInBeam(bnIdx)
                    inBeamNN(bnIdx+1) = true;
                end
            end
            
            elseif obj.params.assignBendingToNonAffectedMonomers
                % assign bending to non damaged monomers 
                inBeamNN = ~affectedInBeam;
            elseif obj.params.assignBendingToNonAffectedMonomersInBeam
                % take only monomer in the beam which are not damaged                
                inBeamNN = inBeam&(~affectedInBeam);                
            end
            
            inBeamAffected = find(inBeamNN); % monomers for which we apply bending elasticity
               
            % Save 
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex = inBeamInds;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam       = affectedInBeam;
            
            
            % Activate bending for affected monomers
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = true;
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = inBeamAffected;            
            
            obj.BreakConnections
            obj.FixDamageBeadsToPlace; 
        end
        
        function RepairDamageEffect(obj)
            % Turn-off bending for affected monomers         
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = false;
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = [];
%             obj.ReFormConnections;
        end
        
        function ReFormConnections(obj)
            % At repair time, re-form the broken connections due to UVC damage
            for cIdx = 1:size(obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads,1)
            obj.handles.framework.objectManager.ConnectParticles(...
                obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads(cIdx,1),...
                obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads(cIdx,2));
            end
        end
        
        function BreakConnections(obj)
            % Break all non-nearest neighbor connections
            cm  = obj.handles.framework.objectManager.GetMembersConnectedParticles(1,'offDiagonals');
            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads = cm;
            if obj.params.breakAllConnectors
            if ~isempty(cm)
                for bIdx = 1:size(cm,1)
                 obj.handles.framework.objectManager.DisconnectParticles(...
                cm(bIdx,1),cm(bIdx,2));
                end
            end
            
            end
            % break all connectors between beads in the beam 
            if obj.params.breakAllDamagedConnectorsInBeam
                  if ~isempty(cm)
                      inBeam = obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam;
                      inBeam = find(inBeam);
                      for bIdx = 1:size(cm,1)
                         if ismember(cm(bIdx,1),inBeam) || ismember(cm(bIdx,2),inBeam)
                         obj.handles.framework.objectManager.DisconnectParticles(...
                        cm(bIdx,1),cm(bIdx,2));
                         end
                      end
                  end
            end
            cm  = obj.handles.framework.objectManager.GetMembersConnectedParticles(1,'offDiagonals');

            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterBeam = cm;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numConnectionsLost =...
                size( obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads,1)-size(cm,1);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentOfConnectedBeadsAfterBeam = size(cm,1)/(obj.params.numMonomers/2);
        end
        
        function FixDamageBeadsToPlace(obj)
            % keep the damaged beads in their place after UVC
            if obj.params.fixDamageMonomersToPlaceAfterBeam
                inBeamInds=obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex;
                obj.handles.framework.objectManager.handles.chain.params.fixedBeadNum = inBeamInds;
            end
            
        end
        
        function UpdateBeamPosition(obj)
            % Update the center of the UVC beam to the polymer's center of mass
            cm = obj.GetPolymerCenterOfMass;
            obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter = cm;
        end                                                                 

        function [rectX,rectY,rectWidth, rectHeight] = GetROI(obj)
            % get roi x, y, width and height
            rectX      = obj.roiPosition(1);
            rectY      = obj.roiPosition(2);
            rectWidth  = obj.roiPosition(3);
            rectHeight = obj.roiPosition(4);
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
               if obj.params.showAdditionalPolymerConnectors
                    cb = obj.params.simulatorParams.chain.connectedBeads;
                    for cIdx =1:size(cb,1)
                    obj.handles.projectedPolymerAdditionalConnectors(cIdx) = line('XData',[initialChainPosition(cb(cIdx,1),1),initialChainPosition(cb(cIdx,2),1)],...
                        'YData',[initialChainPosition(cb(cIdx,1),2),initialChainPosition(cb(cIdx,2),2)],'Color','g','Parent',obj.handles.projectionAxes);
                    end
               end
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
                                initialChainPosition = obj.GetChainPosition;

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
                obj.handles.msdAffected  = line('XData',0,'YData',NaN,'Parent',obj.handles.expensionAxes,'Color','r',...
                    'LineWidth',4,'DisplayName','MSD Affected');
                obj.handles.msdNonAffected = line('XData',0,'YData',NaN,'Parent',obj.handles.expensionAxes,'Color','b',...
                    'Linewidth',4,'DisplayName','MSD Non Affected');
                legend(obj.handles.expensionAxes,get(obj.handles.expensionAxes,'Children'))
            end
                        
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
                % update additional connectors
                if obj.params.showAdditionalPolymerConnectors
                    cb = obj.params.simulatorParams.chain.connectedBeads;
                    for cIdx =1:size(cb,1)
                    obj.handles.projectedPolymerAdditionalConnectors(cIdx) = line('XData',[chainPos(cb(cIdx,1),1),chainPos(cb(cIdx,2),1)],...
                        'YData',[chainPos(cb(cIdx,1),2),chainPos(cb(cIdx,2),2)],'Color','g');
                    end
                end
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
        
        function SaveResults(obj)
            res = obj.results;
            save(fullfile(obj.params.resultsPath,obj.params.resultsFolder,obj.params.resultFileName),'res');
            % update the readme files             
            obj.CreateReadMeFile
        end
        
        function CreateReadMeFile(obj)
            fid = fopen(fullfile(obj.params.resultsPath,obj.params.resultsFolder,'ReadMe.txt'),'w');
            fprintf(fid,'%s\n',num2str(obj.date));
            fprintf(fid,'%s\n%','The list of parameters used in this simulation:');
            fprintf(fid,'%s\n\n','_________________________________________________');
            fnames = fieldnames(obj.params);
            for fIdx = 1:numel(fnames)
                c= class(obj.params.(fnames{fIdx}));
                if strcmpi(c,'double')                    
                        fprintf(fid,'%s%s%s\n',fnames{fIdx} ,': ',num2str(obj.params.(fnames{fIdx})));
                elseif strcmpi(c,'char')
                        fprintf(fid,'%s%s%s\n',fnames{fIdx}, ': ',(obj.params.(fnames{fIdx})));
                end
            end
                                                
            fclose(fid);
        end
    end
    
    methods (Static)
       
        
        function ns = NewResultStruct
            % create a result structure 
            ns = struct('date',[],...
                        'description',[],...
                        'round',[],...
                        'simulation',[],...
                        'dimension',[],...
                        'numRelaxationSteps',[],...
                        'numRecordingSteps',[],...
                        'numBeamSteps',[],...
                        'numRepairSteps',[],...
                        'numBeads',[],...
                        'bendingConst',[],...
                        'springConst',[],...
                        'openingAngle',[],...
                        'connectedBeads',[],...
                        'percentOfConnectedBeads',[],...                        
                        'connectedBeadsAfterBeam',[],...
                        'percentOfConnectedBeadsAfterBeam',[],...
                        'numConnectionsLost',[],...
                        'dt',[],...
                        'ROI',struct('rectX',[],'rectY',[],'rectWidth',[],'rectHeight',[]),...
                        'params',[],...
                        'numBeadsIn',[],...
                        'inBeam',[],...
                        'beadsInIndex',[],...
                        'concentricDensity',[],...                                                        
                        'percentDNALoss',[],...
                        'percentHistoneLoss',[],...
                        'affectedBeadsRadOfExpension',[],...
                        'nonAffectedBeadsRadOfExpension',[],...
                        'chainPos',[]);
        end
    end
end