% scrRunBeamDamagesimulation
% script to run the beam damage simulation
numRelaxationSteps = 500;
numRecordingSteps  = 200;
numBeamSteps       = 300;

saveConfiguration  = false;
loadConfiguration  = false;

% Figures
show3D                = true;
show2D                = true;
showConcentricDensity = true;


% Initialize simulator framework parameters
simulatorParams = SimulationFrameworkParams('showSimulation',show3D,...
                                            'numSteps',numRelaxationSteps,...
                                            'dimension',3,...
                                            'dt',0.001,...
                                            'objectInteraction',false);

% create an open domain
openSpaceForces = ForceManagerParams('lennardJonesForce',false,...
                                     'LJPotentialWidth',0.1,...
                                     'LJPotentialDepth',0.1,...
                                     'diffusionForce',true,...
                                     'diffusionConst',1,...
                                     'mechanicalForce',false,...
                                     'mechanicalForceDirection','out',...
                                     'mechanicalForceCenter',[0 0 0],...
                                     'mechanicalForceMagnitude',0,...
                                     'dt',simulatorParams.simulator.dt);

dp(1)         = DomainHandlerParams('domainShape','open',...
                                    'reflectionType','off',...
                                    'dimension',simulatorParams.simulator.dimension,...
                                    'domainCenter',[0 0 0],...
                                    'forceParams',openSpaceForces,...
                                    'domainWidth',100,...
                                    'dimension',simulatorParams.simulator.dimension);

% % create a chain
chainForces = ForceManagerParams('dt',simulatorParams.simulator.dt,...
                                 'springForce',true,...
                                 'bendingElasticityForce',false,...
                                 'bendingConst',5*simulatorParams.simulator.dimension*openSpaceForces.diffusionConst/(sqrt(3))^2,...
                                 'springConst', 1*simulatorParams.simulator.dimension*openSpaceForces.diffusionConst/(sqrt(3))^2,...
                                 'bendingOpeningAngle',pi,...
                                 'minParticleEqDistance',1);

cp          = ChainParams('numBeads',600,...
    'dimension',simulatorParams.simulator.dimension,...
    'initializeInDomain',3,...
    'forceParams',chainForces,...
    'b',sqrt(3));

% create a cylindrical Beam as a domain
cylinderForces = ForceManagerParams('diffusionForce',false,...
    'lennardJonesForce',false,...
    'morseForce',false);

dp(2)          = DomainHandlerParams('domainShape','cylinder',...
    'reflectionType','off',...
    'domainCenter',[0 0 0],...
    'dimension',simulatorParams.simulator.dimension,...
    'domainWidth',sqrt(cp.numBeads/6)*cp.b/8,...
    'domainHeight',70,...
    'forceParams',cylinderForces);


% Create a sphere for visualization
gSphereForce = ForceManagerParams('lennardJonesForce',false,...
    'diffusionForce',false,...
    'morseForce',false,...
    'mechanicalForce',false);

dp(3)        = DomainHandlerParams('domainWidth',sqrt(cp.numBeads/6)*cp.b,...% radius of Gyration
    'dimension',simulatorParams.simulator.dimension,...
    'domainCenter',[0 0 0],...
    'reflectionType','preserveEnergy',...
    'forceParams',gSphereForce);

% register the object parameters in the simulator framework
simulatorParams.SetDomainParams(dp);
simulatorParams.SetChainParams(cp);

beamDamageParams = BeamDamageParams();

% Initialize simulator framework
r = RouseSimulatorFramework(simulatorParams);
% Run until relaxation time
r.Run
% save configuration by name as chainPos

