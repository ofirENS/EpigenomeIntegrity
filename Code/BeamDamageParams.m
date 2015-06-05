classdef BeamDamageParams<handle %{UNFINISHED}
    
    properties
        numRelaxationSteps@double
        numRecordingSteps@double
        numBeamSteps@double
        roiParams =struct('rectX',[],'rectY',[],'rectWidth',[],'rectHeight');
    end
    
    methods         
        function obj = BeamDamageParams(varargin)

        end
    end
    
end
