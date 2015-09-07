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
        tryMechanicalForceMagnitude@double % mechanical force magnitude to simulate
        tryMechanicalForceCutoff@double   % mechanical spring point force curoff to simulate 
    
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
        incorporateNucleusomes@logical   % simulate nucleusome dynamics  [false true]
        maxNumNucleusomesPerBond@ double 
        
        
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
        
        %__ Repair 
        repairBrokenCrosslinks           % reatach damaged monomers crosslinks after repair stage
        addCrosslinksByDistance          % reattach monomers only if they are located at a distance less than distanceTheresholdToCrosslink
        distanceTheresholdToCrosslink    % at what distance should the monomers re-connect
        turnOffBendingAfterRepair        % turn bending elastivity force off for affected monomers
        removeExclusionVolumeAfterRepair % revmove the volume of exclusion around damaged monomers
                
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
            simulationState        = 'debug'; % options: [debug | simulation]
            
            % Simulation trials 
            % variables to simulate 
            obj.description      = 'Test the expansion and repair of the damaged monomers with 90% crosslinking. Damaged crosslinks are broken after UVC. A volume of exclusion is placed around each damaged monomer with radius 0.68. We test values of exclusion rangin 0.1 to 1. with repair and crosslinks repaires. No Lennard-Jones. Simulation in 2D';
            obj.tryOpeningAngles = []; % obsolete
            obj.tryConnectivity  = [];
            obj.tryNumMonomers   = [];
            obj.tryBendingConst  = [];
            obj.trySpringConst   = [];
            obj.tryMechanicalForceMagnitude = [];%linspace(0.1, 0.5,3);
            obj.tryMechanicalForceCutoff    = [];%linspace(0.5, sqrt(2),3);
            
            %___Simulation parameters___
            obj.numRounds              = 1;%numel(obj.tryMechanicalForceCutoff);
            obj.numSimulationsPerRound = 1;
            obj.numRelaxationSteps     = 10;  % initialization step (burn-in time)
            obj.numRecordingSteps      = 20;  % start recording before UVC beam
            obj.numBeamSteps           = 30; % the steps until repair
            obj.numRepairSteps         = 30;    % repair and relaxation of the fiber
            obj.dt                     = 0.1;
            obj.dimension              = 3;
                                    
            %__Polymer parameters and forces___
            obj.numMonomers           = 1000;
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
            obj.morsePotentialDepth   = 0;
            obj.morsePotentialWidth   = 0;
            obj.morseForceType        = 'repulsive';
            obj.lennardJonesForce     = false;
            obj.LJPotentialWidth      = 0;
            obj.LJPotentialDepth      = 0;
            obj.LJPotentialType       = 'repulsive';
            obj.mechanicalForce       = false;
            obj.mechanicalForceCenter = [];
            obj.mechanicalForceDirection = 'out';
            obj.mechanicalForceMagnitude = 1*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.mechanicalForceCutoff    = 2;
                        
            %___Domain parameters____
            obj.domainRadius          = obj.gyrationRadius;
            obj.domainCenter          = [0 0 0];
            
            %__Beam parameters/damage effect___
            obj.beamRadius                         = obj.gyrationRadius/10;
            obj.beamDamagePeak                     = 0;     % [mu/m/]  zero value coresponds to the focus of the beam 
            obj.beamDamageSlope                    = 1.5;   % slope of the Gaussian shape beam [unitless]
            obj.beamDamageProbThresh               = 1/100; % threshold to determine affected monomers in the UVC beam (obsolete)
            obj.beamHeight                         = 70;    % for 3d graphics purposes
            obj.breakAllDamagedConnectorsInBeam    = true;  % break all connections between damaged monomers in UVC beam  
            obj.breakAllConnectors                 = false; % break all connections in the polymer after UVC
            obj.assignBendingToAffectedMonomers    = false; % assign bending elasticity for affected monomers after UVC
            obj.assignBendingToNonAffectedMonomers = false; % assign bending elasticity for non affected monomers after UVC
            obj.assignBendingToNonAffectedMonomersInBeam = false; % assign bending elasticity for non affected monomers located in the beam after UVC
            obj.fixDamagedMonomersToPlaceAfterBeam       = false; % keep the damaged beads in place after UVC            
            obj.excludeMonomersAroundAffected            = true;
            
            %_ Repair__
            obj.repairBrokenCrosslinks            = true;
            obj.addCrosslinksByDistance           = true;
            obj.distanceTheresholdToCrosslink     = 1;    
            obj.turnOffBendingAfterRepair         = false;
            obj.removeExclusionVolumeAfterRepair  = true;
            
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
            obj.resultsFolder                = 'ROIPostExpansion/ROIByDamaged/ExcludeAroundDamagedMonomers/NoLennardJones/BreakDamagedCrosslinks/10';% result sub-folder name
            cl                               = clock;            
            obj.resultFileName               = sprintf('%s',[num2str(cl(3)),'_',num2str(cl(2)),'_',num2str(cl(1))]); % result file name is given the current time 
            obj.saveAfterEachSimulation      = false;  % save results and create a Readme file after each simulation
            obj.saveAfterEachRound           = true;  % save results and create a Readme file after each simulation round
            obj.saveClassInstance            = false;  % save an instance of the class with the results at the end of operation (usually big files ~50-100Mb).
            
            %___Snapshots____ (for 2D display only)
            obj.numSnapshotsDuringRelaxation = 0;  % unused 
            obj.numSnapshotsDuringRecording  = 5;  % how many snapshots during recording phase 
            obj.numSnapshotsDuringBeam       = 15; % how many snapshots during beam phase 
            obj.numSnapshotsDuringRepair     = 0;  % how many snapshots during repair phase
            
            %__Display real-time parameters___
            obj.show3D                          = false;
            obj.show2D                          = false;
            obj.showDensity                     = false;
            obj.showConcentricDensity           = false;
            obj.showExpensionMSD                = false;
            obj.showExpansionCircle             = false; % show expansion circles (works for 2D only) 
            obj.showAdditionalPolymerConnectors = false; % false speeds up display (for 2D projection image only)
            
            obj.SetSimulationState(simulationState)
            
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
                                                'mechanicalForceCenter',obj.mechanicalForceCenter,...
                                                'mechanicalForceDirection','out',...
                                                'mechanicalForceMagnitude',0.68,...obj.mechanicalForceMagnitude,...
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
                                            'morseForce',false,...
                                            'mechanicalForce',false,...
                                            'springForce',false);
            
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
        
        function SetSimulationState(obj,simulationState)
            % preset fot the debug state
            % set values for debug or simulation 
            if strcmpi(simulationState,'debug')                           
                obj.saveAfterEachSimulation      = false;  % save results and create a Readme file after each simulation
                obj.saveAfterEachRound           = false;  % save results and create a Readme file after each simulation round
                obj.saveClassInstance            = false;
                obj.numSnapshotsDuringRelaxation = 0;  % unused 
                obj.numSnapshotsDuringRecording  = 0; % how many snapshots during recording phase 
                obj.numSnapshotsDuringBeam       = 0; % how many snapshots during beam phase 
                obj.numSnapshotsDuringRepair     = 0;  % how many snapshots during repair phase
                obj.saveRelaxationConfiguration  = false; % unused
                obj.saveEndConfiguration         = false; % unused           
                obj.loadRelaxationConfiguration  = false; % unused
                obj.loadFullConfiguration        = false; % unused
            elseif strcmpi(simulationState,'simulation')
                % keep input parameters
            else
                error('unsupported simulationState option')                
            end
        end
    end
end
