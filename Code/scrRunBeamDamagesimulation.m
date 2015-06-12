% scrRunBeamDamagesimulation
% script to run the beam damage simulation 
            numParticles     = 800;
            dimension        = 3;
            % relxation time 
            % (numParticles*b)^2 / (3*diffusionConst*pi^2)
            numSteps         = 400;
            beamTime         = 200; 
            dt               = 0.01;
            angle0           = pi;
            b                = sqrt(3);
            diffusionConst   = 1;
            springConst      = 1*dimension*diffusionConst./b^2;
            bendingConst     = 2*dimension*diffusionConst./b^2;
            particlePosition = cumsum(sqrt(2*diffusionConst*dt)*randn(numParticles,dimension));

            connectivityMap  = (diag(ones(1,numParticles-1),1)+diag(ones(1,numParticles-1),-1))~=0;
            minParticleDist  = 0;
            fixedParticleNum = [];

            gyrationRadius      = sqrt(numParticles/6)*b; % for large numParticles
            beamRad             = gyrationRadius/8;
            roiRadius           = gyrationRadius/5;
            affectedBeads       = false(numParticles,1);

            % flags
            diffusionFlag      = true;
            springsFlag        = true;
            bendingFlag        = true;

            
params = BeamDamageParams('dimension',3,...
                          'numSteps',400,...
                          'beamTime',200,...
                          'dt',0.01,...
                          'angle0',pi);
                                      
chainParams = ChainParams('numBeads',800,...
                          'dimension',3,...
                          'b',sqrt(3),...
                          'fixedBeamNum',[]);
                      
chainForces = ForceManagerParams('bendingForce',false,...
                                 'springForce',true,...
                                 'bendingConst',2*params.dimension*diffusionConst./chainParams.b^2,...
                                 'springConst',1*params.dimension*diffusionConst./chainParams.b^2,...
                                 'minParticleEqDistance',0);
                             
chainParams.forceParams = chainForces;

domainForces = ForceManagerParams('diffusionForce',true,...
                                  'diffusionConst',1);
                                  
domainParams = DomainHandlerParams('shape','open',...
                                   'dimension',3,...
                                   'forceParams',domainForces);
                          
                                                   