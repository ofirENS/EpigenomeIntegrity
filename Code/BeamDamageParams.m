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
        numMonomers@double % num monomers in polymer
        b@double           % neighboring monomers' distance STD
        connectedMonomers@double  % indicate which monomers are connected (aside from nearest neighbors)
        percentOfConnectedMonomers@double % the percentage of connected monomers (pairs)
        diffusionForce@logical            % on/off [true/false]
        diffusionConst@double
        springForce@logical               % on/off [true/false]
        springConst@double
        bendingForce@logical              % on/off [true/false]
        bendingConst@double
        bendingOpeningAngle@double        % opening angle for bending elasticity potential
        minParticleEqDistance@double      % particle equilibrium distance (for springs)
        gyrationRadius@double
        lennardJonesForce@logical         % on/off [true/false]
        LJPotentialDepth@double
        LJPotentialWidth@double
        LJPotentialType@char
        morseForce@logical                % on/off [true/false]
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char
        mechanicalForce@logical           % on/off [true/false]
        mechanicalForceCenter@double
        mechanicalForceDirection@char
        mechanicalForceMagnitude@double
        mechanicalForceCutoff@double
        
        %___Nucleosomes and forces___
        
        %__Domain params ___
        domainRadius@double
        domainCenter@double
        shutDownDiffusionAfterRelaxationSteps@logical
        
        %__ ROI params ___
        roiWidth@double
        roiHeight@double
        numConcentricBandsInROI@double              % the number of concentric bands in which densities are calculated
        calculateExpansionFromCenterOfMass@logical  % calculate expansion relative to the affected/nonaffected monomers c.m
        calculateExpansionFromBeamCenter@logical    % calculate expansion relative to the beam's center
        calculateExpansionAccordingTo@char          % either 'damaged' or 'nondamaged'
        percentOfMonomersIncludedInROI@double
        
        %___ Beam and beam damage params____
        beamRadius@double                                % radius of the UVC beam 
        beamHeight@double                                % for display purposes 
        beamDamagePeak@double                            % the distance from center the highest pro. for damage
        beamDamageSlope@double                           % the alpha term in exp(-alpha(r-beamDamagePeak)^2)
        beamDamageProbThresh@double                      % threshold below the monomer is damaged [unused]
        assignBendingToAffectedMonomers@logical          % assign bending to the affectd monomers after beam
        assignBendingToNonAffectedMonomers@logical       % assign bending to the non-affected monomers after beam
        assignBendingToNonAffectedMonomersInBeam@logical % assign bending to the non-affected monomers after beam
        breakAllDamagedConnectorsInBeam@logical          % break all connectors in Beam after UVC
        breakAllConnectors@logical                       % break all connectors after UVC
        fixDamagedMonomersToPlaceAfterBeam@logical       % keep damaged monomers in place after UVC
        excludeMonomersAroundAffected@logical             % create a ball around damaged monomers and exclude other monomers from it after UVC
        
        
                
        %__Save and load___
        loadRelaxationConfiguration@logical % unused in this version
        loadFullConfiguration@logical       % unused in this version
        saveRelaxationConfiguration@logical % unused in this version
        saveEndConfiguration@logical        % unused in this version
        resultsFolder@char
        resultsPath@char
        resultFileName@char
        saveAfterEachSimulation@logical
        saveAfterEachRound@logical
        saveClassInstance@logical % save the class at the end of simulations
        
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
        showExpensionMSD@logical    % show a graph of the expansion mean square radius
        showExpansionCircle@logical % show the circle representing the mean expansion of the damaged and non-damaged monomers (plot on the simulation)
        showAdditionalPolymerConnectors@logical % show non nearest neighbors connections (if exist) slows down display
        
        %__Misc. ___
        domainNumbers % indices for the sphere and beam (used for initialization of the polymers)
        
        % relxation time
        % (numParticles*b)^2 / (3*diffusionConst*pi^2)
    end
    
    methods
        
        function obj = BeamDamageParams()% Edit parametes here
            
            % Simulation trials 
            % variables to simulate 
            obj.description      = 'Test the expansion of the damaged and non-damaged monomers with crosslinking. Damaged crosslinks are removed after UVC. Expansion is assigned to both damaged and non-damaged monomers. No Lennard-Jones. Expansion is relative center of mass. No repair steps. Simulation in 2D';
            obj.tryOpeningAngles = [];
            obj.tryConnectivity  = [];%linspace(0,100,6);
            obj.tryNumMonomers   = [];
            obj.tryBendingConst  = [];
            obj.trySpringConst   = [];
                        
            
            %___Simulation parameters___
            obj.numRounds              = 1;% numel(obj.tryConnectivity); 
            obj.numSimulationsPerRound = 1;
            obj.numRelaxationSteps     = 20;  % initialization step (burn-in time)
            obj.numRecordingSteps      = 20;  % start recording before UVC beam
            obj.numBeamSteps           = 1000; % the steps until repair
            obj.numRepairSteps         = 0;  % repair and relaxation of the fiber
            obj.dt                     = 0.1;
            obj.dimension              = 2;
                                    
            %__Polymer parameters and forces___
            obj.numMonomers           = 200;
            obj.b                     = sqrt(obj.dimension);                            
            obj.diffusionForce        = false;
            obj.diffusionConst        = 1;
            obj.shutDownDiffusionAfterRelaxationSteps = true;
            obj.springForce           = true;
            obj.springConst           = 1*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.connectedMonomers     = [];
            obj.percentOfConnectedMonomers = 80; % range: 0 to 100
            obj.minParticleEqDistance = 1;       % sqrt(obj.dimension); % for spring force
            obj.bendingForce          = false;   % (only at initialization)
            obj.bendingConst          = 1*obj.dimension*obj.diffusionConst/obj.b^2;
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
            obj.mechanicalForce       = true;
            obj.mechanicalForceCenter = [0 0 0];
            obj.mechanicalForceDirection = 'out';
            obj.mechanicalForceMagnitude = 10;
            obj.mechanicalForceCutoff    = 5;
                        
            %___Domain parameters____
            obj.domainRadius          = obj.gyrationRadius;
            obj.domainCenter          = [0 0 0];
            
            %__Beam parameters/damage effect___
            obj.beamRadius                         = obj.gyrationRadius/10;
            obj.beamDamagePeak                     = 0;     % [mu/m/]  zero value coresponds to the focus of the beam 
            obj.beamDamageSlope                    = 1.5;   % slope of the Gaussian shape beam [unitless]
            obj.beamDamageProbThresh               = 1/100; % threshold to determine affected monomers in the UVC beam (obsolete)
            obj.beamHeight                         = 70;    % for 3d graphics purposes
            obj.breakAllDamagedConnectorsInBeam    = false;  % break all connections between damaged monomers in UVC beam  
            obj.breakAllConnectors                 = false; % break all connections in the polymer after UVC
            obj.assignBendingToAffectedMonomers    = false; % assign bending elasticity for affected monomers after UVC
            obj.assignBendingToNonAffectedMonomers = false;  % assign bending elasticity for non affected monomers after UVC
            obj.assignBendingToNonAffectedMonomersInBeam = false; % assign bending elasticity for non affected monomers located in the beam after UVC
            obj.fixDamagedMonomersToPlaceAfterBeam       = false; % keep the damaged beads in place after UVC            
            obj.excludeMonomersAroundAffected             = true;
            
            % ___ROI parameters___
            obj.roiWidth                           = obj.gyrationRadius/6; % obsolete used for graphics
            obj.roiHeight                          = obj.roiWidth;         % obsolete used for graphics
            obj.numConcentricBandsInROI            = 15;
            obj.calculateExpansionAccordingTo      = 'damaged'; % either 'damaged' or 'nondamaged'
            obj.calculateExpansionFromCenterOfMass = true;      % calculate expansion dynamically from the center of mass
            obj.calculateExpansionFromBeamCenter   = false;     % calculate expansion from a fixed point (beam center)
            obj.percentOfMonomersIncludedInROI     = 90;        % calculate the ROI such that x percent are in the ROI (circular ROI), this is calculated relatively to the damaged or nondamaged, depending on the value of calculateExpansionAccordingTo
            
            
            %___Save and load___
            obj.saveRelaxationConfiguration  = false; % unused
            obj.saveEndConfiguration         = false; % unused           
            obj.loadRelaxationConfiguration  = false; % unused
            obj.loadFullConfiguration        = false; % unused
            obj.resultsPath                  = fullfile('/home/ofir/Work/ENS/OwnCloud/EpigenomicIntegrity/SimulationResults/'); % top level folder name
            obj.resultsFolder                = 'ROIPostExpansion/ROIByDamaged/BendingDamagedAndNonDamaged/NoLennardJones/BreakDamagedCrosslinks/00';% result sub-folder name
            cl                               = clock;            
            obj.resultFileName               = sprintf('%s',[num2str(cl(3)),'_',num2str(cl(2)),'_',num2str(cl(1))]); % result file name is given the current time 
            obj.saveAfterEachSimulation      = false;  % save results and create a Readme file after each simulation
            obj.saveAfterEachRound           = false;  % save results and create a Readme file after each simulation round
            obj.saveClassInstance            = false;  % save an instance of the class with the results at the end of operation (usually big files ~50-100Mb).
            
            %___Snapshots____ (for 2D display only)
            obj.numSnapshotsDuringRelaxation = 0;  % unused 
            obj.numSnapshotsDuringRecording  = 0; % how many snapshots during recording phase 
            obj.numSnapshotsDuringBeam       = 0; % how many snapshots during beam phase 
            obj.numSnapshotsDuringRepair     = 0;  % how many snapshots during repair phase
            
            %__Display real-time parameters___
            obj.show3D                          = true;
            obj.show2D                          = false;
            obj.showDensity                     = false;
            obj.showConcentricDensity           = false;
            obj.showExpensionMSD                = false;
            obj.showExpansionCircle             = false; % show expansion circles (works for 2D only) 
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
                                             'mechanicalForce',false,...
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
                                                'minParticleEqDistance',obj.minParticleEqDistance,...
                                                'mechanicalForce',obj.mechanicalForce,...
                                                'mechanicalForceCenter',[0 0 0],...
                                                'mechanicalForceDirection','out',...
                                                'mechanicalForceMagnitude',0.1,...
                                                'mechanicalForceCutoff',obj.mechanicalForceCutoff);
                                            
            if ~isempty(obj.percentOfConnectedMonomers)
                if (obj.percentOfConnectedMonomers>100 ||obj.percentOfConnectedMonomers<0)
                    error('percentage of connected monomer must be positive and <100')
                end
                obj.connectedMonomers = []; % remove the previous list of connected monomers, if exist
                n    = obj.numMonomers;
                rp   = randperm(n);
                perc = floor(n*obj.percentOfConnectedMonomers/100);
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
                                            'dimension',obj.dimension,...
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
