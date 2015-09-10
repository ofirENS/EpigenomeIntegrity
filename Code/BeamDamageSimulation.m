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
        state     % simulation state coresponds to the snapshot gallery folder names
        chainPosition % record position of the chain at each step after Recording starts
    end
    
    properties (Access =private)
        snapshotMagazine % holds the indices of the frames in which snapshots are acquired
        roi     % dynamically updated [center, radius]
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
            
            % Prepare save folder with snapshot gallery directory tree
            obj.PrepareSaveDirectoryTree            
            
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
            obj.results.resultStruct(obj.simulationRound,obj.simulation).ROI                = obj.roi;
%             obj.results.resultStruct(obj.simulationRound,obj.simulation).params             = obj.params;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn         = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).lengthIn           = nan(1,numStepsToRecord);% polyymer length in ROI
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex       = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity  = nan(numStepsToRecord,obj.params.numConcentricBandsInROI);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentDNALoss     = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentHistoneLoss = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension    = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension = nan(1,numStepsToRecord);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition           = [];
            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads          = [];  
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentOfConnectedBeads = obj.params.percentOfConnectedMonomers;            
            obj.results.resultStruct(obj.simulationRound,obj.simulation).description             = obj.params.description;
            
        end
        
        function PrepareSaveDirectoryTree(obj)
            % create a directory for the simulation result, the readme file
            % and snapshot gallery
            if obj.params.saveAfterEachRound || obj.params.saveAfterEachSimulation || ...
                    any([obj.params.numSnapshotsDuringBeam,obj.params.numSnapshotsDuringRecording, obj.params.numSnapshotsDuringRepair]>0)
                [~]  = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder));% create the result diretory
                % create simulation Round directories inside, create simulation
                % folders with snapshot folder of recording beam and
                % repair
                
                % Add the snapshot folder tree
                for rIdx = 1:obj.params.numRounds
                    roundFolderName = sprintf('%s%s','Round',num2str(rIdx));
                    [~] = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName));% create recording snapshots
                                        
                    for sIdx = 1:obj.params.numSimulationsPerRound
                        simulationFolderName = sprintf('%s%s','Simulation',num2str(sIdx));  
                        % create the snapshot folder with Recording, Beam,
                        % and Repair 
                        [~] = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName));% create recording snapshots                                                                        
                        [~] = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName,'Snapshots','Recording'));% create recording snapshots
                        [~] = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName,'Snapshots','Beam'));% create recording snapshots 
                        [~] = mkdir(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName,'Snapshots','Repair'));% create recording snapshots

                    end
                end                
            end
        end
        
        function PrepareSnapshotMagazine(obj)
            % prpare the step indices for the snapshots
            st = obj.params.numRelaxationSteps+1:(obj.params.numRelaxationSteps+obj.params.numRecordingSteps);
            if ~isempty(st)
              obj.snapshotMagazine  = unique(round(linspace(st(1),st(end),obj.params.numSnapshotsDuringRecording)));
            end
            
            st = (obj.params.numRelaxationSteps+obj.params.numRecordingSteps+1):(obj.params.numRelaxationSteps+obj.params.numRecordingSteps+obj.params.numBeamSteps);
            if ~isempty(st)
              obj.snapshotMagazine  = [obj.snapshotMagazine unique(round(linspace(st(1),st(end),obj.params.numSnapshotsDuringBeam)))];
            end
            
            st = (obj.params.numRelaxationSteps+obj.params.numRecordingSteps+obj.params.numBeamSteps+1):(obj.params.numRelaxationSteps+obj.params.numRecordingSteps+obj.params.numBeamSteps+obj.params.numRepairSteps);
            if ~isempty(st)
             obj.snapshotMagazine = [obj.snapshotMagazine unique(round(linspace(st(1),st(end),obj.params.numSnapshotsDuringRepair)))];
            end
            obj.snapshotMagazine = unique(obj.snapshotMagazine);
        end
        
        function Run(obj)
            obj.simulationRound =0;            
            for rIdx = 1:obj.params.numRounds            
                obj.PreRoundActions
                for sIdx =1:obj.params.numSimulationsPerRound                    
                    obj.PreSimulationActions;                    
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
            obj.state = 'Relaxation';
            obj.handles.framework.Run;% run the framework            
%             obj.UpdateROIPosition;
%             obj.UpdateBeamPosition;
            obj.step = obj.params.numRelaxationSteps+1;% handles.framework.simulationData.step;
        end
        
        function RunRecordingSteps(obj)
            % Start Recording 
            obj.state = 'Recording';
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps;
            
            for sIdx =1:obj.params.numRecordingSteps
                obj.handles.framework.Step;
                stepIdx  = obj.step -obj.params.numRelaxationSteps;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx) = obj.GetChainPosition;% record position of the chain           
                obj.UpdateGraphics;                                
                obj.Snapshot
                obj.step = obj.step+1;%obj.handles.framework.simulationData.step;
            end
        end
         
        function RunBeamSteps(obj)
            obj.state = 'Beam';
            % Increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps;  
            
            % Activate the UVC beam
            obj.ApplyDamageEffect; 
            % record neighboring position for affected monomers 
            
            for sIdx =1:obj.params.numBeamSteps
                obj.ApplyExclusionByVolume
                obj.handles.framework.Step
%                 chainPos = obj.ApplyExclusionByVolume;
                
                stepIdx  = obj.step-obj.params.numRelaxationSteps;              
                % record position of the chain
                obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx) = obj.GetChainPosition;
                obj.UpdateGraphics
                obj.Snapshot               
                % calculate the radius of expansion for damaged and
                % non-damaged monomers
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx),...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);  
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                obj.step = obj.step+1;                                
            end
            obj.SetROI
        end
              
        function RunRepairSteps(obj)
            obj.state = 'Repair';
            % increase the step count 
            obj.handles.framework.params.simulator.numSteps = obj.params.numRelaxationSteps+obj.params.numRecordingSteps+...
                                                              obj.params.numBeamSteps+obj.params.numRepairSteps;  
                        
        
            obj.TurnOffAffectedBeadsGraphics
                        
            for sIdx =1:obj.params.numRepairSteps
                
                obj.handles.framework.Step
                obj.RepairDamageEffect; % perform repair 
                
                stepIdx  = obj.step-obj.params.numRelaxationSteps;                
                obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx) = obj.GetChainPosition;% record position of the chain   
                obj.Snapshot
                obj.UpdateGraphics
                  % calculate the radius of expansion for damaged and
                % non-damaged monomers
                [affectedMSDcm,nonAffectedMSDcm] = obj.CalculateBeadsRadiusOfExpension(obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx),...
                                                   obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);  
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
                
                obj.step = obj.step+1;                     
%                 [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI;
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                
%                 obj.CalculateBeadsRadiusOfExpension(chainPos,...
%                                                    obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam);
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;
           
            end
        end
        
        function SetROI(obj)
           % Set the ROI according to the expansion of the monomers after beam phase      
           % take the last 10% of the beam phase and calculate the
                % mean                 
            sIdx = round(0.9*(obj.params.numRecordingSteps+obj.params.numBeamSteps)):obj.params.numRecordingSteps+obj.params.numBeamSteps;
            if strcmpi(obj.params.calculateExpansionAccordingTo,'damaged')                
                obj.roi.radius = mean(obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(sIdx));
                obj.results.resultStruct(obj.simulationRound,obj.simulation).ROI = obj.roi.radius;
            elseif strcmpi(obj.params.calculateExpansionAccordingTo,'nondamaged')
                obj.roi.radius = mean(obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(sIdx));
                obj.results.resultStruct(obj.simulationRound,obj.simulation).ROI = obj.roi.radius;
            else
                error('Unsupported option. The options for calculateExpansionAccordingTo are "damaged" or "nondamaged" only')
            end
            % set the concentric ROi, the circle in which we calculate
            % concentric densities.
            % the concentric ROi is always determined by the non affected
            % monomers' radius of expansion 
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricROI = ...
                mean(obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(sIdx));
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
                obj.simulation = 0;
                obj.simulationRound = obj.simulationRound+1;% increase counter
        end
        
        function PostRoundActions(obj)
            % Actions performed at the end of each round   
            % calculate the number and density of monomers in the ROI                        
            if obj.params.saveAfterEachRound
                % Save results to result folder                
                obj.SaveResults
            end
        end
        
        function PreSimulationActions(obj)
            
            obj.simulation = obj.simulation+1;% increase counter                
%             obj.params     = BeamDamageParams;% initialize a new param class
            
            % Add connected beads 
            if ~isempty(obj.params.tryConnectivity)
              obj.params.percentOfConnectedMonomers = obj.params.tryConnectivity(obj.simulationRound);
              obj.params.InitializeParamClasses;            
            end
            
            if ~isempty(obj.params.tryMechanicalForceMagnitude)
              obj.params.mechanicalForceMagnitude = obj.params.tryMechanicalForceMagnitude(obj.simulationRound);
              obj.params.InitializeParamClasses;            
            end
            
            if ~isempty(obj.params.tryMechanicalForceCutoff)
              obj.params.mechanicalForceCutoff = obj.params.tryMechanicalForceCutoff(obj.simulationRound);
              obj.params.InitializeParamClasses;                
            end
            
            % crete/remove previous framework
            if isfield(obj.handles,'framework')
                delete(obj.handles.framework)
            end
             obj.results.resultStruct(obj.simulationRound,obj.simulation).params        = obj.params;
            obj.handles.framework = RouseSimulatorFramework(obj.params.simulatorParams);
            
            % prepare the result structure
            obj.PrepareResultStruct
            
            cl = clock;
            cl = [num2str(cl(4)),':' num2str(cl(5)), ':', num2str(cl(6))];
            sprintf('%s%s%s%s%s%s','Started Round ', num2str(obj.simulationRound), ' Simulation ', num2str(obj.simulation),' at time: ', cl)
            obj.PrepareSnapshotMagazine;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition = zeros(obj.params.numMonomers,3,obj.params.numRecordingSteps+obj.params.numBeamSteps+obj.params.numRepairSteps);
             
        end
        
        function PostSimulationActions(obj)
            % Because the ROi is determined only after beam phase, we calculate the following properties at the end of each simulation 
             
            % Calculate the number of monomers in the ROI before UVC
            for stepIdx = 1:obj.params.numRecordingSteps
                 [inROI,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI(stepIdx);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric;
                [dnaLengthIn] = obj.GetDnaLengthInROI(inROI,stepIdx);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).lengthIn(stepIdx) = dnaLengthIn;
            end
            
            % Calculate the densities in the ROI after UVC beam stated
            for stepIdx = obj.params.numRecordingSteps:obj.params.numBeamSteps+obj.params.numRecordingSteps+obj.params.numRepairSteps
                % get the polymer's center of mass at that point                 
               [~,~,numMonomersIn,monomersInConcentric] = obj.GetMonomerDensityInROI(stepIdx);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
                 [dnaLengthIn] = obj.GetDnaLengthInROI(inROI,stepIdx);
                obj.results.resultStruct(obj.simulationRound,obj.simulation).lengthIn(stepIdx) = dnaLengthIn;
            end
            
            % calculate the radius of expansion for the damaged monomers
            % before UVC 
            for stepIdx = 1:obj.params.numRecordingSteps
                [affectedMSDcm,nonAffectedMSDcm]=obj.CalculateBeadsRadiusOfExpension( obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx),...
                    obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam );
                obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx)    = affectedMSDcm;
                obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx) = nonAffectedMSDcm;                
            end
                        
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
            chainPos   = obj.GetChainPosition; 
            inBeam     = obj.handles.framework.handles.classes.domain.InDomain(chainPos,obj.params.domainNumbers.beam);
            inBeamInds = find(inBeam);            
                
        end
        
        function [affectedMSDcm,nonAffectedMSDcm] = CalculateBeadsRadiusOfExpension(obj,chainPos,inBeam)
            % Calculate the mean radius of expension for affected beads and non
            % affected beads. chainPos is an numBedsXdimension vector of position
            % inBeam is a numBeadsX1 binary vector indicating the beads affected
            % (true) and non affected (false)

            if obj.params.calculateExpansionFromCenterOfMass
                cmInBeam         = mean(chainPos(inBeam,:));  % center of mass for particles in beam
                cmOutBeam        = mean(chainPos(~inBeam,:)); % center of mass for paarticles out of the beam
                % calculate the radius of the expansion circle according to
                % the percentage of monomers inside
                r   = sort(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), cmInBeam).^2,2)));
                if ~isempty(r)
                    affectedMSDcm    = r(round(numel(r)*obj.params.percentOfMonomersIncludedInROI/100));
                    
                else
                    affectedMSDcm = nan;
                end
                
                r                    = sort(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), cmOutBeam).^2,2)));
                if ~isempty(r)
                    nonAffectedMSDcm = r(round(numel(r)*obj.params.percentOfMonomersIncludedInROI/100));
                else
                    nonAffectedMSDcm =nan;
                end
               
%                affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), cmInBeam).^2,2)));
%                nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), cmOutBeam).^2,2)));
            elseif obj.params.calculateExpansionFromBeamCenter
                affectedMSDcm    = mean(sqrt(sum(bsxfun(@minus, chainPos(inBeam,:), obj.beamCenterPosition).^2,2)));
                nonAffectedMSDcm = mean(sqrt(sum(bsxfun(@minus, chainPos(~inBeam,:), obj.beamCenterPosition).^2,2)));
            else 
                error('unsupported option')
            end
             stepIdx  = obj.step-obj.params.numRelaxationSteps;

        end
        
        function [inROI,inRoiInds,numMonomersIn,monomersInConcentric] = GetMonomerDensityInROI(obj,stepIdx)% finish concentric
                % Find the number of monomers in the rectangular ROI                
%                 obj.UpdateROIPosition;
                [inROI, inRoiInds,numMonomersIn] = obj.FindMonomersInROI(stepIdx);
                [monomersInConcentric]           = obj.GetMonomersInRoiConcentric(stepIdx);  
                % calculate the radius from the center of mass                                
                %----
                obj.results.resultStruct(obj.simulationRound,obj.simulation).numBeadsIn(stepIdx)          = numMonomersIn;            
        end
        
        function dnaLength = GetDnaLengthInROI(obj,inROI, stepIdx)
            % find the length of the bonds in the ROI
            chainPos = obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx);
            cm       = mean(chainPos,1);
            chainPos(:,1) = chainPos(:,1)-cm(1);
            chainPos(:,2) = chainPos(:,2)-cm(2);
            chainPos(:,3) = chainPos(:,3)-cm(3);
            
            % get all monomers in ROI 
            dnaLength = 0;
            for rIdx =2:numel(inROI)
                if inROI(rIdx-1)&& inROI(rIdx)
                    dnaLength = dnaLength+norm(chainPos(rIdx,:)-chainPos(rIdx-1,:));
                elseif ~inROI(rIdx-1)&& inROI(rIdx)
                    % find intersection with circle
                   a    = sum((chainPos(rIdx,:)-chainPos(rIdx-1,:)).^2);
                   b    = 2*dot(chainPos(rIdx-1,:),chainPos(rIdx,:)-chainPos(rIdx-1,:));
                   c    = sum(chainPos(rIdx-1,:).^2 -obj.roi.radius^2);
                   t(1) = (-b +sqrt(b^2 -4*a*c))/(2*a);
                   t(2) = (-b -sqrt(b^2 -4*a*c))/(2*a);
                   % take the positive t
                   t= min(t(t>0));
                   % calculate the distance from the intersection point to
                   % the monomer
                   dnaLength = dnaLength + norm((chainPos(rIdx,:)-chainPos(rIdx-1,:))*(1-t));
                elseif inROI(rIdx-1,:) && ~inROI(rIdx,:)
                    a = sum((chainPos(rIdx,:)-chainPos(rIdx-1,:)).^2);
                   b = 2*dot(chainPos(rIdx-1,:),chainPos(rIdx,:)-chainPos(rIdx-1,:));
                   c = sum(chainPos(rIdx-1,:).^2 -obj.roi.radius^2);
                   t(1) = (-b +sqrt(b^2 -4*a*c))/(2*a);
                   t(2) = (-b -sqrt(b^2 -4*a*c))/(2*a);
                   % take the positive t
                   t= min(t(t>0));
                   % calculate the distance from the intersection point to
                   % the monomer
                   if ~isempty(t)% needs checkup 
                    dnaLength = dnaLength + norm((chainPos(rIdx,:)-chainPos(rIdx-1,:))*t);
                   end
                end
            end
        end
        
        function cm = GetPolymerCenterOfMass(obj)
            % polymer center of mass in any dimension 
            chainPos = obj.handles.framework.objectManager.curPos;
            cm       = mean(chainPos);
        end
        
        function UpdateROIPosition(obj)
            % Update the square roi position such that its geometrical center is at the polymer center of mass
%             polymerCenterOfMass = obj.GetPolymerCenterOfMass;
%             obj.roiPosition = [polymerCenterOfMass(1)-obj.params.roiWidth/2,...
%                                polymerCenterOfMass(2)-obj.params.roiHeight/2,...
%                                obj.params.roiWidth,...
%                                obj.params.roiHeight];           
        end                        
        
        function [inROI, inRoiInds,numMonomersIn]= FindMonomersInROI(obj,stepIdx)
            % Get the number and indices of the monomers inside the ROI at
            % step sIdx from the end of relaxation step
            
            chainPos = obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx);

            cm = mean(chainPos,1);
            r  = pdist2mex(chainPos',cm','euc',[],[],[]);

            inROI         = r<obj.roi.radius;   
            inRoiInds     = find(inROI);
            numMonomersIn = sum(inROI);
        end
        
        function chainPosition = GetChainPosition(obj)
            % Get chain position in the last step 
            chainPosition= obj.handles.framework.objectManager.curPos;
        end
        
        function [monomersInConcentric] = GetMonomersInRoiConcentric(obj,stepIdx)
             % Calculate the density as a function of the distance from the roi center
                chainPos = obj.results.resultStruct(obj.simulationRound,obj.simulation).chainPosition(:,:,stepIdx);
                cm = mean(chainPos,1);
                r  = pdist2mex(chainPos',cm','euc',[],[],[]);
                % divide the ragion into concentric circles 
                monomersInConcentric = zeros(1,obj.params.numConcentricBandsInROI);
                R = linspace(0,obj.roi.radius,obj.params.numConcentricBandsInROI+1);
                
                for rIdx = 1:numel(R)-1
                    bandArea = pi*(R(rIdx+1)^2 - R(rIdx)^2);
                    monomersInConcentric(rIdx) = sum(r>=R(rIdx) & r<R(rIdx+1))./bandArea;
                end
            
            obj.results.resultStruct(obj.simulationRound,obj.simulation).concentricDensity(stepIdx,:) = monomersInConcentric; 
            
        end
        
        function ApplyDamageEffect(obj)
            % Apply damage effect to monomers and their nearest neighbors
            obj.beamCenterPosition = obj.GetPolymerCenterOfMass;
            chainPos               = obj.GetChainPosition;
            
            obj.UpdateBeamPosition 
            obj.UpdateBeamGraphics;
            obj.HeighlightDamagedMonomers;
            
            [inBeam, ~] = obj.FindMonomersInBeam;% get indices and indicators for affected monomers and their nearest-neighbors
            
            % Calculate the distance of each monomer in the beam to the
            % beam's center 
            r = pdist2mex(chainPos',obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter','euc',[],[],[]);
            
            % Calculate the probability of each monomer to be affected
            damageProb = exp(-obj.params.beamDamageSlope.*(r-obj.params.beamDamagePeak).^2);
            
            % use theresholding to define the actual monomers affected 
            affectedInBeam = (damageProb>rand(size(r,1),1) & inBeam);
            inBeamInds     = find(affectedInBeam); 
            
            affectedMonomers = false(size(affectedInBeam));
            
            if obj.params.assignBendingToAffectedMonomers
             % Assign bending also to neighbors of damaged monomers, to
             % create continuity
             affectedMonomers= affectedInBeam;
            for bnIdx=1:size(affectedInBeam,1)-1
                if affectedInBeam(bnIdx)
                    affectedMonomers(bnIdx+1) = true;
                end
            end
            
            elseif obj.params.assignBendingToNonAffectedMonomers
                % assign bending to non damaged monomers 
                affectedMonomers = ~inBeam;
                
            elseif obj.params.assignBendingToNonAffectedMonomersInBeam
                % take only monomer in the beam which are not damaged   
                affectedMonomers = inBeam&(~affectedInBeam);
%                 inBeamNN = inBeam&(~affectedInBeam);       
                for bnIdx=1:size(affectedInBeam,1)-1
                  if affectedInBeam(bnIdx)&& ~inBeam(bnIdx+1)
                    affectedMonomers (bnIdx+1) = true;
                  end
                end
            end
            
            affectedMonomersInds = find(affectedMonomers); % monomers for which we apply bending elasticity
               
            % Save 
            obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex = inBeamInds;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam       = affectedInBeam;
            
            obj.ApplyBending(affectedMonomersInds)
            obj.ApplyExclusionByVolume;
            obj.BreakCrosslinks
            obj.FixDamagedMonomersToPlace; 
        end
        
        function ApplyBending(obj,affectedMonomersInds)
            % Activate bending for affected monomers
            if obj.params.assignBendingToNonAffectedMonomers || obj.params.assignBendingToNonAffectedMonomersInBeam ||...
                    obj.params.assignBendingToAffectedMonomers
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = true;
            obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = affectedMonomersInds; 
            end
        end
        
        function ApplyExclusionByVolume(obj)
            % create an exclusion region around affected monomer with
            % harmonic pushing force (half spring)
            if obj.params.excludeMonomersAroundAffected
             obj.handles.framework.objectManager.handles.chain.params.forceParams.mechanicalForce = true;
             obj.handles.framework.objectManager.handles.chain.params.forceParams.mechanicalForceCenter = ...
                obj.handles.framework.objectManager.handles.chain.position.cur(...
                 obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam,:);
            end
        end
        
        function RepairDamageEffect(obj)
            % the effect after repair stage is over (15 min post UVC)
            if obj.params.turnOffBendingAfterRepair
            % Turn-off bending for affected monomers         
             obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingElasticityForce   = false;
             obj.handles.framework.objectManager.handles.chain.params.forceParams.bendingAffectedParticles = [];
            end
            
            if obj.params.removeExclusionVolumeAfterRepair
                % remove the exclusion volume around affected monomers
                obj.handles.framework.objectManager.handles.chain.params.forceParams.mechanicalForce = false;
            end
            
            if obj.params.repairBrokenCrosslinks
                obj.ReformConnections
            end      
        end
        
        function ReformConnections(obj)
            % At repair time, re-create links 
            if obj.params.addCrosslinksByDistance
                % add cross links up to the initial level by adding them to
                % close monomers in the ROI
                initialConnectivity = obj.params.percentOfConnectedMonomers;
                currentConnectivity = ceil(2*size(obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair,1)./obj.params.numMonomers *100); 

                % reconnect beads initially in beam 
             if initialConnectivity>ceil(currentConnectivity)
                inBeam  = obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex;
                % find the distances to those monomers  
                partDist = obj.handles.framework.objectManager.particleDist;               
                % check connectivity for particles in beam
               for ibIdx = 1:numel(inBeam)
                % get particle distance
               neighbPartDist = (partDist(inBeam(ibIdx),:));
                % assign Inf to nearest neighbors
               neighbPartDist([max([1,inBeam(ibIdx)-1]),inBeam(ibIdx), min([obj.params.numMonomers,inBeam(ibIdx)+1])])= Inf;                
                % find those particles who are
                % at a distance lower than the threshold for connecting
               f = find(neighbPartDist<obj.params.distanceTheresholdToCrosslink);
                if ~isempty(f)
                    % connect to the smaller one 
                    [~,m] = min(neighbPartDist(f));
                    % connect f(mf) and inBeam(ibIdx) if not already
                    % connected
                    if ~obj.handles.framework.objectManager.connectivity(f(m),inBeam(ibIdx))
                    obj.handles.framework.objectManager.ConnectParticles(f(m),inBeam(ibIdx));% connnect
                    obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair = ...
                        [obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair; [f(m),inBeam(ibIdx)]];% update connection list
                    currentConnectivity = ceil(2*size(obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair,1)./obj.params.numMonomers *100); 
                    end
                end
                if currentConnectivity>=initialConnectivity 
                    break
                end
               end
            end
            
%             for cIdx = 1:size(obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads,1)
%             obj.handles.framework.objectManager.ConnectParticles(...
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads(cIdx,1),...
%                 obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads(cIdx,2));            
%             end
            
            cm  = obj.handles.framework.objectManager.GetMembersConnectedParticles(1,'offDiagonals');
            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair = cm;
            end
        end
        
        function BreakCrosslinks(obj)
            % Break all non-nearest neighbor connections after UVC beam-shot
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
            obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeadsAfterRepair = cm;
            obj.results.resultStruct(obj.simulationRound,obj.simulation).numConnectionsLost =...
                size( obj.results.resultStruct(obj.simulationRound,obj.simulation).connectedBeads,1)-size(cm,1);
            obj.results.resultStruct(obj.simulationRound,obj.simulation).percentOfConnectedBeadsAfterBeam = size(cm,1)/(obj.params.numMonomers/2);            
        end
        
        function FixDamagedMonomersToPlace(obj)
            % keep the damaged beads in their place after UVC
            if obj.params.fixDamagedMonomersToPlaceAfterBeam
                inBeamInds=obj.results.resultStruct(obj.simulationRound,obj.simulation).beadsInIndex;
                obj.handles.framework.objectManager.handles.chain.params.fixedBeadNum = inBeamInds;
            end
            
        end
        
        function UpdateBeamPosition(obj)
            % Update the center of the UVC beam to the polymer's center of mass
            cm = obj.GetPolymerCenterOfMass;
            obj.handles.framework.handles.classes.domain.params(obj.params.domainNumbers.beam).domainCenter = cm;
        end                                                                 

        function [rectX,rectY,rectWidth, rectHeight] = GetROI(obj)% replace
            % get roi x, y, width and height
            rectX      = obj.roiPosition(1);
            rectY      = obj.roiPosition(2);
            rectWidth  = obj.roiPosition(3);
            rectHeight = obj.roiPosition(4);
        end
        
        function InitializeGraphics(obj)
            
            % create figure for the projection in the x-y plane            
           
            
            if obj.params.show3D
            % insert the projection to the 3d axes
%             [rectX,rectY,rectWidth, rectHeight] = obj.GetROI;
            mainAxes = obj.handles.framework.simulationGraphics.handles.graphical.mainAxes;
            t        = linspace(0,2*pi,20);            
%             obj.handles.projPlane3D = patch([rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
%                 [rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)],...
%                 'r', 'Parent',mainAxes, 'FaceAlpha',0.5);
            obj.handles.projPlane3D = patch(cos(t),sin(t),...
                'r', 'Parent',mainAxes, 'FaceAlpha',0.5,'Visible','off');
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
%             obj.UpdateROIGraphics
%             obj.UpdateDensityGraphics;
            obj.UpdateProjectedPolymerGraphics
%             obj.UpdateConcentricDensityGraphics
            obj.UpdateExpensionGraphics
            obj.HeighlightDamagedMonomers
            obj.ShowRadiusOfExpansion;
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
                                   
        function UpdateROIGraphics(obj)%obsolete
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
        
        function Snapshot(obj)
            % in 2D only for now
            if ~isempty(obj.snapshotMagazine)
            if obj.step == obj.snapshotMagazine(1)
                obj.snapshotMagazine = obj.snapshotMagazine(2:end);% remove the first index
                
                f        = figure('Visible','off','Units','norm');
                a        = axes('Parent',f);
                chainPos = obj.GetChainPosition;
                % plot the chain 
                line('XData',chainPos(:,1),...
                    'YDAta',chainPos(:,2),...
                    'Marker','o',...
                    'MarkerFaceColor','none',...
                    'MarkerEdgeColor','b',...
                    'markerSize',7,...
                    'Parent',a,...
                    'LineStyle','-');
              
                    cb = obj.params.simulatorParams.chain.connectedBeads;
                    for cIdx =1:size(cb,1)
                    line('XData',[chainPos(cb(cIdx,1),1),chainPos(cb(cIdx,2),1)],...
                         'YData',[chainPos(cb(cIdx,1),2),chainPos(cb(cIdx,2),2)],...
                         'Color','g','Parent',a);
                    end
              % plot additional connectors
               damagedMonomers = chainPos(obj.results.resultStruct(obj.simulationRound,obj.simulation).inBeam ,:);
               line('XDAta',damagedMonomers(:,1),'YData',damagedMonomers(:,2),'Marker','o',...
                    'MarkerFaceColor','r',...
                    'MarkerEdgeColor','r',...
                    'MarkerSize',7,...
                    'Parent',a,...
                    'LineStyle','none');
                fileName = sprintf('%s',['snapshotStep',num2str(obj.step),'.jpg']);
                
                % Save the snapshot image 
                roundFolderName      = sprintf('%s%s','Round',num2str(obj.simulationRound));
                simulationFolderName = sprintf('%s%s','Simulation',num2str(obj.simulation));
                saveas(f,fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName,'Snapshots',obj.state,fileName));
                   
                close(f);% close the hidden figure to save memory
            end
            end
        end
        
        function HeighlightDamagedMonomers(obj)
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
        
        function ShowRadiusOfExpansion(obj)
            if obj.params.showExpansionCircle && (obj.params.dimension==2) && obj.params.show3D
                    stepIdx          = obj.step-obj.params.numRelaxationSteps;
                    t                = 0:(2*pi/30):2*pi;                    
                    if obj.params.calculateExpansionFromCenterOfMass                        
                        circleCenter = obj.GetPolymerCenterOfMass;
                    elseif obj.params.calculateExpansionFromBeamCenter
                        circleCenter = obj.beamCenterPosition;
                    end
                        
                    damagedExpRad    = obj.results.resultStruct(obj.simulationRound,obj.simulation).affectedBeadsRadOfExpension(stepIdx);
                    xDamaged         = damagedExpRad*cos(t)+circleCenter(1);
                    yDamaged         = damagedExpRad*sin(t)+circleCenter(2); 
                    nonDamagedExpRad = obj.results.resultStruct(obj.simulationRound,obj.simulation).nonAffectedBeadsRadOfExpension(stepIdx);
                    xNonDamaged      = nonDamagedExpRad*cos(t)+circleCenter(1);
                    yNonDamaged      = nonDamagedExpRad*sin(t)+circleCenter(2); 
                    
                if ~isfield(obj.handles,'expansionCircleDamaged')% if circles do not exist yet
                
                    % create expansion circles
                    
                    % damaged monomers
                    obj.handles.expansionCircleDamaged = line('XData',xDamaged,'YData',yDamaged,'Color','m',...
                        'Parent',obj.handles.framework.simulationGraphics.handles.graphical.mainAxes,'LineWidth',2);
                    % non-damaged monomers           
                    obj.handles.expansionCircleNonDamaged = line('XData',xNonDamaged,'YData',yNonDamaged,'Color','b',...
                        'Parent',obj.handles.framework.simulationGraphics.handles.graphical.mainAxes,'LineWidth',2);
                else % if they exists 
                    % update their position 
                    set(obj.handles.expansionCircleDamaged,'XData',xDamaged,'YData',yDamaged,...
                        'Parent',obj.handles.framework.simulationGraphics.handles.graphical.mainAxes);
                    set(obj.handles.expansionCircleNonDamaged,'XData',xNonDamaged,'YData',yNonDamaged,...
                        'Parent',obj.handles.framework.simulationGraphics.handles.graphical.mainAxes);                    
                end                                                
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
            save(fullfile(obj.params.resultsPath,obj.params.resultsFolder,obj.params.resultFileName),'res','-v7.3');
            % update the readme files             
            obj.CreateReadMeFile
        end
        
        function CreateReadMeFile(obj)
            roundFolderName      = sprintf('%s%s','Round',num2str(obj.simulationRound));
            simulationFolderName = sprintf('%s%s','Simulation',num2str(obj.simulation)); 
            fid                  = fopen(fullfile(obj.params.resultsPath,obj.params.resultsFolder,roundFolderName,simulationFolderName,'ReadMe.txt'),'w');
            fprintf(fid,'%s\n',num2str(obj.date));
            fprintf(fid,'%s\n%','The list of parameters used in this simulation:');
            fprintf(fid,'%s\n\n','_________________________________________________');
            fnames = fieldnames(obj.params);
            for fIdx = 1:numel(fnames)
                c= class(obj.params.(fnames{fIdx}));
                if strcmpi(c,'double')       
                        if size(obj.params.(fnames{fIdx}),1)>1
                            fprintf(fid,'%s%s\n',fnames{fIdx} ,': [');
                            for sIdx = 1:size(obj.params.(fnames{fIdx}),1)
                                fprintf(fid,'%s\n',num2str(obj.params.(fnames{fIdx})(sIdx,:)));
                            end
                            fprintf(fid,'%s\n',']');
                        else
                            
                        fprintf(fid,'%s%s%s\n',fnames{fIdx} ,': ',num2str(obj.params.(fnames{fIdx})));
                        end
                elseif strcmpi(c,'char')
                        fprintf(fid,'%s%s%s\n',fnames{fIdx}, ': ',(obj.params.(fnames{fIdx})));
                elseif strcmpi(c,'logical')        
                       fprintf(fid,'%s%s%s\n',fnames{fIdx}, ': ',num2str(double(obj.params.(fnames{fIdx}))));
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
                        'dt',[],...
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
                        'connectedBeadsAfterRepair',[],...
                        'percentOfConnectedBeadsAfterBeam',[],...
                        'numConnectionsLost',[],...                        
                        'ROI',struct('center',[],'radius',[]),...
                        'concentricROI',[],...
                        'params',[],...
                        'numBeadsIn',[],...
                        'lengthIn',[],...
                        'inBeam',[],...
                        'beadsInIndex',[],...
                        'concentricDensity',[],...
                        'percentDNALoss',[],...
                        'percentHistoneLoss',[],...
                        'affectedBeadsRadOfExpension',[],...
                        'nonAffectedBeadsRadOfExpension',[],...
                        'chainPosition',[]);
        end
    end
end