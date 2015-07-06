classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties
        simulatorParams = SimulationFrameworkParams
        numRelaxationSteps@double
        numRecordingSteps@double
        numBeamSteps@double
        numRounds@double
        numSimulationsPerRound@double
        dimension@double
        dt@double
                        
        
        
        % Chain parameters and forces
        numMonomers@double % num monomers in chain
        b@double % bead distance std
        diffusionForce@logical
        diffusionConst@double
        springForce@logical
        springConst@double
        bendingForce@logical
        bendingConst@double
        bendingOpeningAngle@double
        connectedMonomers@logical        
        minParticleEqDistance@double % particle eq distance (springs)
        gyrationRadius@double
        lennardJonesForce@logical
        LJPotentialDepth@double
        LJPotentialWidth@double
        
        % Domain 
        domainRadius@double
        domainCenter@double
        
        % ROI 
        roiRadius@double
        numConcentricBandsInROI@double% the number of concentric bands in which densities are calculated
        
        % Beam  
        beamRadius@double
        beamHeight@double 
        fracOfMonomersAffected@double 
        
        % Save and load        
        loadRelaxationConfiguration@logical
        loadFullConfiguration@logical
        saveRelaxationConfiguration@logical
        saveEndConfiguration@logical
        
        % Display
        showSimulation@logical

        
        % relxation time
        % (numParticles*b)^2 / (3*diffusionConst*pi^2)
    end
    
    methods
        
        function obj = BeamDamageParams()
            
            % Simulation parameters
            obj.numRounds              = 1; 
            obj.numSimulationsPerRound = 1;
            
            obj.numRelaxationSteps = 100;
            obj.numRecordingSteps  = 100;
            obj.numBeamSteps       = 500;
            obj.dt                 = 0.01;
            obj.dimension          = 2;
            
            % Polymer parameters and forces
            obj.numMonomers           = 500;
            obj.b                     = sqrt(obj.dimension);                            
            obj.diffusionForce        = false;
            obj.diffusionConst        = 1;
            obj.springForce           = true;
            obj.springConst           = obj.dimension*obj.diffusionConst/obj.b^2;
            obj.minParticleEqDistance = obj.b; % for spring force
            obj.bendingForce          = false;% (only at initialization)
            obj.bendingConst          = 2*obj.dimension*obj.diffusionConst/obj.b^2;
            obj.bendingOpeningAngle   = pi;            
            obj.gyrationRadius        = sqrt(obj.numMonomers/6)*obj.b;
            obj.lennardJonesForce     = false;
            obj.LJPotentialWidth      = 0.1;
            obj.LJPotentialDepth      = 0.1;            
            % Domain parameters

            obj.domainRadius          = obj.gyrationRadius/2;
            obj.domainCenter          = [0 0 0];
            
            % Beam parameters
            obj.beamRadius = obj.gyrationRadius/4;
            obj.beamHeight = 70;
            obj.fracOfMonomersAffected = 1; % [0-1] fraction of the monomers to affected by the beam 
            
            % ROI parameters                        
            obj.roiRadius               = obj.gyrationRadius/4;
            obj.numConcentricBandsInROI = 10;
            
            % Save and load
            obj.saveRelaxationConfiguration  = false;
            obj.saveEndConfiguration         = false;
            
            obj.loadRelaxationConfiguration  = false;
            obj.loadFullConfiguration        = false;
            
            % Display parameters
            obj.showSimulation  = false;
            
            % initializ the classes
            obj.InitializeParamClasses
        end
        
        function InitializeParamClasses(obj)
            
            % Initialize parameter classes 
            
            % Initialize simulator framework parameters
            obj.simulatorParams = SimulationFrameworkParams('showSimulation',obj.showSimulation,...
                                                            'numSteps',obj.numRelaxationSteps,...
                                                            'dimension',obj.dimension,...
                                                            'dt',obj.dt,...
                                                            'diffusionConst',obj.diffusionConst,...
                                                            'objectInteraction',false);

            % create an open domain
            domainForces = ForceManagerParams('lennardJonesForce',obj.lennardJonesForce,...
                'LJPotentialWidth',obj.LJPotentialWidth,...
                'LJPotentialDepth',obj.LJPotentialDepth,...
                'diffusionForce',obj.diffusionForce,...
                'diffusionConst',obj.diffusionConst,...
                'dt',obj.dt);
            
            dp(1)         = DomainHandlerParams('domainShape','sphere',...
                                                'reflectionType','preserveEnergy',...
                                                'dimension',obj.dimension,...
                                                'domainCenter',obj.domainCenter,...
                                                'forceParams',domainForces,...
                                                'domainWidth',obj.domainRadius,...
                                                'dimension',obj.dimension);
            
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
            
            
%             % Create a sphere for visualization
%             gSphereForce = ForceManagerParams('lennardJonesForce',false,...
%                 'diffusionForce',false,...
%                 'morseForce',false,...
%                 'mechanicalForce',false);
%             
%             dp(3)        = DomainHandlerParams('domainWidth',obj.gyrationRadius,...
%                 'dimension',obj.dimension,...
%                 'domainCenter',[0 0 0],...
%                 'reflectionType','preserveEnergy',...
%                 'forceParams',gSphereForce);
            
            % register the object parameters in the simulator framework
            obj.simulatorParams.SetDomainParams(dp);
            obj.simulatorParams.SetChainParams(cp);
            
        end
        
    end
end
