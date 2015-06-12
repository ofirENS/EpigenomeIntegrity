classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties
        numParticles@double        
        numRelaxationSteps@double
        numRecordingSteps@double
        numBeamSteps@double
        dimension@double
        dt@double
        connectivityMat@logical
        forceParams = ForceManagerParams;
        roiRadius@double
        beamRad@double
        gyrationRadius@double
        particlePosition@double 
        loadConfiguration@logical
        saveAfterRelaxationTime@logical
        saveAfterBeamTime@logical
        showSimulation@logical
% relxation time 
% (numParticles*b)^2 / (3*diffusionConst*pi^2)
    end
    
    methods         
        function obj = BeamDamageParams(varargin)
             obj.SetInputParams(vrargin)
        end
        
        function SetInputParams(obj,varargin)
            if mod(numel(varargin),2)~=0
                error(' input argument must come in name value pairs')
            else
                for vIdx = 1:numel(varargin)
                    obj.(varargin{2*vIdx-1}) = varargin{2*vIdx};
                end
            end
        end
    end
    
end
