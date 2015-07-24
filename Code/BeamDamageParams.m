classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties% Do not edit this section manually
        % These properties are place-holders for parameter values. Do not edit this section
        % manually. To set parameter values see BeamDamageParams function
        % in the method section of this class
        
        % ___ General ___
        simulatorParams = SimulationFrameworkParams
        description@char              % description of the simulation performed
        numRelaxationSteps@double     % number of relaxation steps before recording starts
        numRecordingSteps@double      % number of recording steps before UVC beam shot
        numBeamSteps@double           % number of steps after UVC beam shot
        numRepairSteps@double         % number of steps after repair mechanism finished
        numRounds@double              % number of simulation rounds
        numSimulationsPerRound@double % number of simulations in each round
        dimension@double              % dimension 
        dt@double                     % time step
                        
        %___Simulation trials ___
        tryOpeningAngles@double  % opening angle values to simulate
        tryConnectivity@double   % percent of connected monomers to simulate 
        tryNumMonomers@double    % number of monomers to simulate
        tryBendingConst@double   % bending constant to simulate
        trySpringConst@double    % spring constant to simulate 
        
        %__Chain and chain forces parameters ___
        numMonomers@double % num monomers in chain
        b@double           % bead distance STD
        connectedMonomers@double
        percentOfConnectedMonomers@double % the percentage of connected monomers (pairs)
        diffusionForce@logical % on/off [true/false]
        diffusionConst@double
        springForce@logical    % on/off [true/false]
        springConst@double
        bendingForce@logical   % on/off [true/false]
        bendingConst@double
        bendingOpeningAngle@double     
        minParticleEqDistance@double % particle eq distance (springs)
        gyrationRadius@double
        lennardJonesForce@logical
        LJPotentialDepth@double
        LJPotentialWidth@double
        LJPotentialType@char
        morseForce@logical
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char
        mechanicalForce@logical
        mechanicalForceCenter@double
        mechanicalForceDirection@char
        mechanicalForceMagnitude@double
        
        %__Domain params ___
        domainRadius@double
        domainCenter@double
        shutDownDiffusionAfterRelaxationSteps@logical
        
        %__ ROI params ___
        roiWidth@double
        roiHeight@double
        numConcentricBandsInROI@double% the number of concentric bands in which densities are calculated
        
        %___ Beam  and beam damage params____
        beamRadius@double                     % radius of the UVC beam 
        beamHeight@double                     % for display purposes 
        beamDamagePeak@double                 % the distance from center the highest pro. for damage
        beamDamageSlope@double                % the alpha term in exp(-alpha(r-beamDamagePeak)^2)
        beamDamageProbThresh@double           % threshold below the monomer is damaged [unused]
        assignBendingToAffectedMonomers@logical          % assign bending to the affectd monomers after beam
        assignBendingToNonAffectedMonomers@logical       % assign bending to the non-affected monomers after beam
        assignBendingToNonAffectedMonomersInBeam@logical % assign bending to the non-affected monomers after beam
        breakAllDamagedConnectorsInBeam@logical      % break all connectors in Beam after UVC
        breakAllConnectors@logical            % break all connectors after UVC
        fixDamageMonomersToPlaceAfterBeam@logical % keep damaged monomers in place after UVC
        calculateMSDFromCenterOfMass@logical  % calculate expansion relative to the affected/nonaffected monomers c.m
        calculateMSDFromBeamCenter@logical    % calculate expansion relative to the beam's center
        
        %__Save and load___
        loadRelaxationConfiguration@logical
        loadFullConfiguration@logical
        saveRelaxationConfiguration@logical
        saveEndConfiguration@logical
        resultsFolder@char
        resultsPath@char
        resultFileName@char
        saveAfterEachSimulation@logical
        saveAfterEachRound@logical
        
        %__ Snapshots__
        numSnapshotsDuringRelaxation@double
        numSnapshotsDuringRecording@double
        numSnapshotsDuringBeam@double
        numSnapshotsDuringRepair@double
        
        %__Display params ___
        show3D@logical
        show2D@logical
        showDensity@logical
        showConcentricDensity@logical
        showExpensionMSD@logical   % show a graph of the expansion mean square radius
        showExpansionCircle@logical % show the circle representing the mean expansion of the damaged and non-damaged monomers (plot on the simulation)
        showAdditionalPolymerConnectors@logical % show non nearest neighbors connections (if exist) slows down display
        
        %__Misc. ___
        domainNumbers % indices for the sphere and beam 
        
        % relxation time
        % (numParticles*b)^2 / (3*diffusionConst*pi^2)
    end
    
    methods
        
        function obj = BeamDamageParams()% Edit parametes here
            
            % Simulation trials 
            % variables to simulate 
            obj.description      = 'Test the expansion of the damaged and non-damaged monomers, with crosslinking. crosslinks not removed. damaged monomers are assigned bending. Expansion is relative to beam center. No repair steps. Simulation in 2D';
            obj.tryOpeningAngles = [];
            obj.tryConnectivity  = linspace(0,1,11);
            obj.tryNumMonomers   = [];
            obj.tryBendingConst  = [];
            obj.trySpringConst   = [];
                        
            
            %___Simulation parameters___
            obj.numRounds              = 1;%numel(obj.tryConnectivity); 
            obj.numSimulationsPerRound = 1;
            obj.numRelaxationSteps     = 10;   % initialization step (burn-in time)
            obj.numRecordingSteps      = 100;  % start recording before UVC beam
            obj.numBeamSteps           = 10; % the steps until repair
            obj.numRepairSteps         = 10;    % repair and relaxation of the fiber
            obj.dt                     = 0.1;
            obj.dimension              = 2;
                                    
            %__Polymer parameters and forces___
            obj.numMonomers           = 100;
            obj.b                     = sqrt(obj.dimension);                            
            obj.diffusionForce        = true;
            obj.diffusionConst        = 1;
            obj.shutDownDiffusionAfterRelaxationSteps = true;
            obj.springForce           = true;
            obj.springConst           = 1*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.connectedMonomers     = [];
            obj.percentOfConnectedMonomers = 0.7;
            obj.minParticleEqDistance = 1;     % sqrt(obj.dimension); % for spring force
            obj.bendingForce          = false; % (only at initialization)
            obj.bendingConst          = 0.1*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.bendingOpeningAngle   = pi;
            obj.gyrationRadius        = sqrt(obj.numMonomers/6)*obj.b;
            obj.morseForce            = false;
            obj.morsePotentialDepth   = 0.01;
            obj.morsePotentialWidth   = 0.01;
            obj.morseForceType        = 'repulsive';
            obj.lennardJonesForce     = false;
            obj.LJPotentialWidth      = 0.1;
            obj.LJPotentialDepth      = 1;
            obj.LJPotentialType       = 'repulsive';
            obj.mechanicalForce       = false;
            obj.mechanicalForceCenter = [0 0 0];
            obj.mechanicalForceDirection = 'out';
            obj.mechanicalForceMagnitude = 0.3;
                        
            %___Domain parameters____
            obj.domainRadius          = obj.gyrationRadius/2;
            obj.domainCenter          = [0 0 0];
            
            %__Beam parameters/damage effect___
            obj.beamRadius                         = obj.gyrationRadius/10;
            obj.beamDamagePeak                     = 0;     % [mu/m/]  zero value coresponds to the focus of the beam 
            obj.beamDamageSlope                    = 1.5;   % slope of the Gaussian shape beam [unitless]
            obj.beamDamageProbThresh               = 1/100; % threshold to determine affected monomers in the UVC beam (obsolete)
            obj.beamHeight                         = 70;    % for 3d graphics purposes
            obj.breakAllDamagedConnectorsInBeam    = false;  % break all connections between damaged monomers in UVC beam  
            obj.breakAllConnectors                 = false; % break all connections in the polymer after UVC
            obj.assignBendingToAffectedMonomers    = true; % assign bending elasticity for affected monomers after UVC
            obj.assignBendingToNonAffectedMonomers = false;  % assign bending elasticity for non affected monomers after UVC
            obj.assignBendingToNonAffectedMonomersInBeam = false; % assign bending elasticity for non affected monomers located in the beam after UVC
            obj.fixDamageMonomersToPlaceAfterBeam  = false; % keep the damaged beads in place after UVC
            obj.calculateMSDFromCenterOfMass       = true; % calculate expansion dynamically from the center of mass
            obj.calculateMSDFromBeamCenter         = false;  % calculate expansion from a fixed point (beam center)
            
            
            % ___ROI parameters___
            obj.roiWidth                = obj.gyrationRadius/6;
            obj.roiHeight               = obj.roiWidth;
            obj.numConcentricBandsInROI = 15;
            
            %___Save and load___
            obj.saveRelaxationConfiguration  = false;
            obj.saveEndConfiguration         = false;            
            obj.loadRelaxationConfiguration  = false;
            obj.loadFullConfiguration        = false;
            obj.resultsPath                  = fullfile('/home/ofir/Copy/ENS/EpigenomicIntegrity/SimulationResults/');
            obj.resultsFolder                = 'Test';%'ExpansionTest/CrossLinked/NoLennardJones/BendingDamaged/00';           
            cl                               = clock;            
            obj.resultFileName               = sprintf('%s',[num2str(cl(3)),'_',num2str(cl(2)),'_',num2str(cl(1))]); 
            obj.saveAfterEachSimulation      = false;
            obj.saveAfterEachRound           = false;
            
            %___Snapshots [unfinished]____
            obj.numSnapshotsDuringRelaxation = 0;
            obj.numSnapshotsDuringRecording  = 0;
            obj.numSnapshotsDuringBeam       = 0;
            obj.numSnapshotsDuringRepair     = 0;
            
            %__Display real-time parameters___
            obj.show3D                = false;
            obj.show2D                = false;
            obj.showDensity           = false;
            obj.showConcentricDensity = false;
            obj.showExpensionMSD      = false;
            obj.showExpansionCircle   = false; % show expansion circles (works for 2D only) 
            obj.showAdditionalPolymerConnectors = false; % false speeds up display (for 2D projection image only)
            
            %___Initializ the classes___
            obj.InitializeParamClasses
        end
        
        function InitializeParamClasses(obj)%Do not edit this section manually
            
            % Initialize parameter classes 
            
            % Initialize simulator framework parameters
            obj.simulatorParams = SimulationFrameworkParams('showSimulation',obj.show3D,...
                                                            'numSteps',obj.numRelaxationSteps,...
                                                            'dimension',obj.dimension,...
                                                            'dt',obj.dt,...
                                                            'diffusionConst',obj.diffusionConst,...
                                                            'objectInteraction',false);

            % create an open domain
            domainForces = ForceManagerParams('lennardJonesForce',obj.lennardJonesForce,...
                                             'LJPotentialWidth',obj.LJPotentialWidth,...
                                             'LJPotentialDepth',obj.LJPotentialDepth,...
                                             'LJPotentialType',obj.LJPotentialType,...
                                             'morseForce',obj.morseForce,...
                                             'morsePotentialDepth',obj.morsePotentialDepth,...
                                             'morsePotentialWidth',obj.morsePotentialWidth,...
                                             'morseForceType',obj.morseForceType,...
                                             'diffusionForce',obj.diffusionForce,...
                                             'diffusionConst',obj.diffusionConst,...
                                             'mechanicalForce',obj.mechanicalForce,...
                                             'mechanicalForceCenter',[0 0 0],...
                                             'mechanicalForceDirection','out',...
                                             'mechanicalForceMagnitude',0.1,...
                                             'dt',obj.dt,...
                                             'minParticleEqDistance',obj.minParticleEqDistance);
            
            dp(1)         = DomainHandlerParams('domainShape','sphere',...
                                                'reflectionType','preserveEnergy',...
                                                'dimension',obj.dimension,...
                                                'domainCenter',obj.domainCenter,...
                                                'forceParams',domainForces,...
                                                'domainWidth',obj.domainRadius,...
                                                'dimension',obj.dimension);
            obj.domainNumbers.sphere = 1; 
            % % create the polymer
            polymerForces  = ForceManagerParams('dt',obj.dt,...
                                                'springForce',obj.springForce,...
                                                'bendingElasticityForce',obj.bendingForce,...
                                                'bendingConst',obj.bendingConst,...
                                                'springConst', obj.springConst,...
                                                'bendingOpeningAngle',obj.bendingOpeningAngle,...
                                                'minParticleEqDistance',obj.minParticleEqDistance);
                                            
            if ~isempty(obj.percentOfConnectedMonomers)
                n    = obj.numMonomers;
                rp   = randperm(n);
                perc = floor(n*obj.percentOfConnectedMonomers);
            % make sure it is divisible by 2
            if mod(perc,2)~=0
                perc = max([perc-1,0]);
            end
            for pIdx = 1:(perc)/2
                nnFlag = false; % exit flag
                while ~nnFlag 
                     rp1   = randperm(n);   
                     if abs(rp(pIdx)-rp1(1))~=1
                         nnFlag = true;
                     end
                end
                obj.connectedMonomers(pIdx,:) = [rp(pIdx), rp1(1)];
            end
            
            end
            cp          = ChainParams('numBeads',obj.numMonomers,...
                                      'dimension',obj.dimension,...
                                      'connectedBeads',obj.connectedMonomers,...
                                      'initializeInDomain',1,...
                                      'minBeadDistance',obj.minParticleEqDistance,...
                                      'forceParams',polymerForces,...
                                      'b',obj.b);
            
            % create a cylindrical Beam as a domain
            beamForces = ForceManagerParams('diffusionForce',false,...
                                            'lennardJonesForce',false,...
                                            'morseForce',false);
            
            dp(2)     = DomainHandlerParams('domainShape','cylinder',...
                                            'reflectionType','off',...
                                            'domainCenter',[0 0 0],...
                                            'dimension',3,...
                                            'domainWidth',obj.beamRadius,...
                                            'domainHeight',obj.beamHeight,...
                                            'forceParams',beamForces);
           
            obj.domainNumbers.beam   = 2;
            % register the object parameters in the simulator framework
            obj.simulatorParams.SetDomainParams(dp);
            obj.simulatorParams.SetChainParams(cp);
            
        end
        
    end
end
