classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties
        simulatorParams = SimulationFrameworkParams
        numRelaxationSteps@double
        numRecordingSteps@double
        numBeamSteps@double
        numRepairSteps@double
        numRounds@double
        numSimulationsPerRound@double
        dimension@double
        dt@double
                        
        % simulation trials 
        tryOpeningAngles@double 
        tryConnectivity@double  
        tryNumMonomers@double   
        tryBendingConst@double  
        trySpringConst@double   
        
        % Chain parameters and forces
        numMonomers@double % num monomers in chain
        b@double % bead distance std
        connectedMonomers@double
        diffusionForce@logical
        diffusionConst@double
        springForce@logical
        springConst@double
        bendingForce@logical
        bendingConst@double
        bendingOpeningAngle@double     
        minParticleEqDistance@double % particle eq distance (springs)
        gyrationRadius@double
        lennardJonesForce@logical
        LJPotentialDepth@double
        LJPotentialWidth@double
        morseForce@logical
        morsePotentialDepth@double
        morsePotentialWidth@double
        morseForceType@char
        
        % Domain 
        domainRadius@double
        domainCenter@double
        
        % ROI 
        roiWidth@double
        roiHeight@double
        numConcentricBandsInROI@double% the number of concentric bands in which densities are calculated
        
        % Beam  
        beamRadius@double           % radius of the UVC beam 
        beamHeight@double           % for display purposes 
        beamDamagePeak@double       % the distance from center the highest pro. for damage
        beamDamageSlope@double      % the alpha term in exp(-alpha(r-beamDamagePeak)^2)
        beamDamageProbThresh@double % threshold below the monomer is damaged
        
        % Save and load        
        loadRelaxationConfiguration@logical
        loadFullConfiguration@logical
        saveRelaxationConfiguration@logical
        saveEndConfiguration@logical
        
        % Display
        show3D@logical
        show2D@logical
        showDensity@logical
        showConcentricDensity@logical
        showExpensionMSD@logical
        
        % Misc. 
        domainNumbers % indices for the sphere and beam 
        
        % relxation time
        % (numParticles*b)^2 / (3*diffusionConst*pi^2)
    end
    
    methods
        
        function obj = BeamDamageParams()
            
            % Simulation trials 
            % variables to simulate 
            obj.tryOpeningAngles = 1;
            obj.tryConnectivity  = 1;
            obj.tryNumMonomers   = [];
            obj.tryBendingConst  = [];
            obj.trySpringConst   = [];
            
            % Simulation parameters
            obj.numRounds              = max([numel(obj.tryOpeningAngles),1]); 
            obj.numSimulationsPerRound = max([numel(obj.tryConnectivity),1]);
            obj.numRelaxationSteps     = 100; % initialization step (burn-in time)
            obj.numRecordingSteps      = 500; % start recording before UVC beam
            obj.numBeamSteps           = 300;% the steps until repair
            obj.numRepairSteps         = 200;% repair and relaxation of the fiber
            obj.dt                     = 1/10;
            obj.dimension              = 3;
                                    
            % Polymer parameters and forces
            obj.numMonomers           = 400;
            obj.b                     = sqrt(obj.dimension);                            
            obj.diffusionForce        = true;
            obj.diffusionConst        = 1;
            obj.springForce           = true;
            obj.springConst           = obj.dimension*obj.diffusionConst/obj.b^2;
            obj.connectedMonomers     = [];
            obj.minParticleEqDistance = obj.b; % for spring force
            obj.bendingForce          = false; % (only at initialization)
            obj.bendingConst          = obj.dimension*obj.diffusionConst/obj.b^2;
            obj.bendingOpeningAngle   = pi;            
            obj.gyrationRadius        = sqrt(obj.numMonomers/6)*obj.b;
            obj.morseForce            = true;
            obj.morsePotentialDepth   = 0.01;
            obj.morsePotentialWidth   = 0.01;
            obj.morseForceType        = 'repulsive';
            obj.lennardJonesForce     = false;
            obj.LJPotentialWidth      = 0.05;
            obj.LJPotentialDepth      = 0.05; 
            
            % Domain parameters
            obj.domainRadius          = obj.gyrationRadius/2;
            obj.domainCenter          = [0 0 0];
            
            % Beam parameters
            obj.beamRadius           = obj.gyrationRadius/6;
            obj.beamDamagePeak       = 0;%obj.beamRadius/5 ; % in mu/m
            obj.beamDamageSlope      = 0.01; % unitless
            obj.beamDamageProbThresh = 1/100;% threshold to determine affected monomers in the UVC beam (obsolete)
            obj.beamHeight           = 70; % for 3d graphics purposes
            
            % ROI parameters                        
            obj.roiWidth                = obj.gyrationRadius/4;
            obj.roiHeight               = obj.roiWidth;
            obj.numConcentricBandsInROI = 15;
            
            % Save and load
            obj.saveRelaxationConfiguration  = false;
            obj.saveEndConfiguration         = false;            
            obj.loadRelaxationConfiguration  = false;
            obj.loadFullConfiguration        = false;
            
            % Display on-line parameters
            obj.show3D                = true;
            obj.show2D                = true;
            obj.showDensity           = false;
            obj.showConcentricDensity = false;
            obj.showExpensionMSD      = true;
            
            % initializ the classes
            obj.InitializeParamClasses
        end
        
        function InitializeParamClasses(obj)
            
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
                                             'morseForce',obj.morseForce,...
                                             'morsePotentialDepth',obj.morsePotentialDepth,...
                                             'morsePotentialWidth',obj.morsePotentialWidth,...
                                             'morseForceType',obj.morseForceType,...
                                             'diffusionForce',obj.diffusionForce,...
                                             'diffusionConst',obj.diffusionConst,...
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
            
            cp          = ChainParams('numBeads',obj.numMonomers,...
                                      'dimension',obj.dimension,...
                                      'initializeInDomain',1,...
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
