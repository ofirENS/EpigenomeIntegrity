classdef BeamDamageSimulation<handle %[UNFINISHED]
    properties
        params
        handles
        particlePosition
    end
    
    methods
        
        function obj = BeamDamageSimulation(params)
            % parameters
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
            
                                    
        end
        
    end
end