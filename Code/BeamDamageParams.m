classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties% Do not edit this section
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
        minCrossLinkedParticlesEqDistance@double % cross-linked particle equilibrium distance
        crossLinkedParticlesSpringConst@double  % spring constant for cross-links
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
        beamDamageSTD@double                             % the std term in exp(-(1/std)^2(r-beamDamagePeak)^2)      
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
        maxCrossLinksPerMonomer           % maximum number of cross links in addition to the linear nearest neighbors connection
        turnOffBendingAfterRepair        % turn bending elastivity force off for affected monomers
        removeExclusionVolumeAfterRepair % revmove the volume of exclusion around damaged monomers
        repairProbRate                   % rate of repair (Poissonian)
        encounterDistance                % the distance at which two monomers are considerd to have encounterd. used for comparision of structure before and after UV
        
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
        
        function obj = BeamDamageParams()% Edit parameters here
            simulationState        = 'debug'; % options: [debug | simulation]
            
            % Simulation trials            
            obj.description      = 'Test the expansion and repair of the damaged monomers with 90% crosslinking. Damaged crosslinks are broken after UVC. A volume of exclusion is placed around each damaged monomer with radius 0.55. We test values of exclusion rangin 0.1 to 1. with repair and crosslinks repaires. No Lennard-Jones. Simulation in 2D';
            obj.tryOpeningAngles            = []; % obsolete
            obj.tryConnectivity             = [];
            obj.tryNumMonomers              = [];
            obj.tryBendingConst             = [];
            obj.trySpringConst              = [];
            obj.tryMechanicalForceMagnitude = [];
            obj.tryMechanicalForceCutoff    = [];
            
            %___Simulation parameters___
            obj.numRounds              = 1;   % numel(obj.tryConnectivity);
            obj.numSimulationsPerRound = 1;
            obj.numRelaxationSteps     = 50; % initialization step (burn-in time)
            obj.numRecordingSteps      = 150; % start recording before UVC beam
            obj.numBeamSteps           = 400; % the steps until repair
            obj.numRepairSteps         = 300; % repair and relaxation of the fiber
            obj.dt                     = 0.1;
            obj.dimension              = 2;
                                    
            %__Polymer parameters and forces___
            obj.numMonomers                       = 500;
            obj.percentOfConnectedMonomers        = 90; % range: 0 to 100            
            obj.b                                 = sqrt(obj.dimension);                            
            obj.diffusionForce                    = true;
            obj.diffusionConst                    = 1;
            obj.shutDownDiffusionAfterRelaxationSteps = true;
            obj.springForce                       = true;
            obj.springConst                       = obj.dimension*obj.diffusionConst/obj.b^2;
            obj.connectedMonomers                 = [];            
            obj.minParticleEqDistance             = 0.5*obj.b.*ones(obj.numMonomers); % sqrt(obj.dimension); % for spring force
            obj.minCrossLinkedParticlesEqDistance = 0;%0.5*obj.minParticleEqDistance(1);       % minimal distance for cross-linked particles   
            obj.crossLinkedParticlesSpringConst   = obj.springConst; % spring constant for crosslinks
            % Assign minimal contraction distance for connected monomers'
            % springs
            for cIdx = 1:size(obj.connectedMonomers,1)
                obj.minParticleEqDistance(obj.connectedMonomers(cIdx,1),obj.connectedMonomers(cIdx,2)) = obj.minCrossLinkedParticlesEqDistance;
                obj.minParticleEqDistance(obj.connectedMonomers(cIdx,2),obj.connectedMonomers(cIdx,1)) = obj.minCrossLinkedParticlesEqDistance;                
            end
            
            obj.bendingForce             = false;   % (only at initialization)
            obj.bendingConst             = obj.dimension*obj.diffusionConst/obj.b^2;
            obj.bendingOpeningAngle      = pi;
            obj.gyrationRadius           = sqrt(obj.numMonomers/6)*obj.b;
            obj.morseForce               = false;
            obj.morsePotentialDepth      = 0;
            obj.morsePotentialWidth      = 0;
            obj.morseForceType           = 'repulsive';
            obj.lennardJonesForce        = false;
            obj.LJPotentialWidth         = 0;
            obj.LJPotentialDepth         = 0;
            obj.LJPotentialType          = 'repulsive';
            obj.mechanicalForce          = false;
            obj.mechanicalForceCenter    = [];
            obj.mechanicalForceDirection = 'out';
            obj.mechanicalForceMagnitude = 0.5*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.mechanicalForceCutoff    = 1*obj.minParticleEqDistance(1);
                        
            %___Domain parameters____
            obj.domainRadius          = obj.gyrationRadius;
            obj.domainCenter          = [0 0 0];
            
            %__Beam parameters___
            obj.beamRadius           = obj.gyrationRadius/20;% defines a unit of 1 mu/m in true space
            obj.beamDamagePeak       = 0;     % [mu/m/]  zero value coresponds to the focus of the beam 
            obj.beamDamageSTD        = obj.beamRadius/5;   % slope of the Gaussian shape beam [unitless]            
            obj.beamHeight           = 70;    % for 3d graphics purposes
            
            %__Damage effect ___
            obj.breakAllDamagedConnectorsInBeam          = true;  % break all connections between damaged monomers in UVC beam  
            obj.breakAllConnectors                       = false; % break all connections in the polymer after UVC
            obj.assignBendingToAffectedMonomers          = false; % assign bending elasticity for affected monomers after UVC
            obj.assignBendingToNonAffectedMonomers       = false; % assign bending elasticity for non affected monomers after UVC
            obj.assignBendingToNonAffectedMonomersInBeam = false; % assign bending elasticity for non affected monomers located in the beam after UVC
            obj.fixDamagedMonomersToPlaceAfterBeam       = false; % keep the damaged beads in place after UVC            
            obj.excludeMonomersAroundAffected            = true;  % create an exclusion sphere around damaged monomers with spring pushing force
            
            %__ Repair__
            obj.repairBrokenCrosslinks            = true; % at repair stage, do we reintroduce crosslinks
            obj.addCrosslinksByDistance           = true; % if yes, do we do it by distance
            obj.distanceTheresholdToCrosslink     = obj.minParticleEqDistance(1)/5;% 0.1*obj.b; % what is the distance between monomers for which we cross-link them after repair  
            obj.maxCrossLinksPerMonomer           = 2;  % in addition to the linear NN connection
            obj.turnOffBendingAfterRepair         = false;
            obj.removeExclusionVolumeAfterRepair  = true; % the pushing force around damaged monomers
            obj.repairProbRate                    = 5; % (Poissonian) the number of event is determined with rate repairProbRate*dt
            obj.encounterDistance                 = obj.minParticleEqDistance(1)/5;% 0.1*obj.b; % the distance at which monomers are considered to have encounter
            
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
            obj.resultsFolder                = 'TestRepair/ROIByDamaged/ExcludeAroundDamagedMonomers/BreakDamagedCrosslinks/02';% result sub-folder name
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
            obj.showExpansionCircle             = true; % show expansion circles (works for 2D only) 
            obj.showAdditionalPolymerConnectors = true; % false speeds up display (for 2D projection image only)
            
            % set state to debug or simulation
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
            % the spherical domain
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
                                                'mechanicalForceMagnitude',obj.mechanicalForceMagnitude,...
                                                'mechanicalForceCutoff',obj.mechanicalForceCutoff);
                                                       
            % randomize connected monomers other than linear connectivity           
            if obj.percentOfConnectedMonomers>0
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
                     if abs(rp(pIdx)-rp1(1))>1
                         nnFlag = true;
                     end
                end
                obj.connectedMonomers(pIdx,:) = [rp(pIdx), rp1(1)];
            end            
            end
            
            % assign min distance for connected monomers
            for cIdx = 1:size(obj.connectedMonomers,1)
                obj.minParticleEqDistance(obj.connectedMonomers(cIdx,1),obj.connectedMonomers(cIdx,2)) = obj.minCrossLinkedParticlesEqDistance;
                obj.minParticleEqDistance(obj.connectedMonomers(cIdx,2),obj.connectedMonomers(cIdx,1)) = obj.minCrossLinkedParticlesEqDistance;
            end
            
            
            cp          = ChainParams('numBeads',obj.numMonomers,...
                                      'dimension',obj.dimension,...
                                      'connectedBeads',obj.connectedMonomers,...
                                      'initializeInDomain',obj.domainNumbers.sphere,...
                                      'minBeadDistance',obj.minParticleEqDistance,...
                                      'forceParams',polymerForces,...
                                      'b',obj.b);
            
            % create a cylindrical Beam as a domain
            beamForces = ForceManagerParams('diffusionForce',false,...
                                            'lennardJonesForce',false,...
                                            'morseForce',false,...
                                            'mechanicalForce',false,...
                                            'springForce',false);
            % Representing the laser beam
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
