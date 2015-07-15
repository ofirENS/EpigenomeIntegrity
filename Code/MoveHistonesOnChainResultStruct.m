classdef MoveHistonesOnChainResultStruct<handle
    % collect results from MoveHistonesOnChain
    properties 
        date
        numRounds
        numSimulationsPerRound
        resultStruct=struct;
    end
    
    methods
        function obj = MoveHistonesOnChainResultStruct(numRounds,numSimulationsPerRound)
            
            cl       = clock;
            obj.date = sprintf('%s',[num2str(cl(3)),'/',num2str(cl(2)),'/',num2str(cl(1))]);
            obj.numRounds              = numRounds;
            obj.numSimulationsPerRound = numSimulationsPerRound;
            obj.resultStruct = obj.NewStruct;
            for rIdx = 1:numRounds
                for sIdx = 1:numSimulationsPerRound
            obj.resultStruct(rIdx,sIdx) = obj.NewStruct;
                 
                end
            end

        end
        
        function mbIn= CalculateMeanNumBeadIn(obj)
            
            for rIdx = 1:obj.numRounds
                for sIdx = 1:obj.numSimulationsPerRound
                    bIn(sIdx,:) = obj.resultStruct(rIdx,sIdx).numBeadsIn;
                end
                mbIn(rIdx,:) = mean(bIn);
            end
            
        end
        
        function expRad= CalcualteMeanAffectedRadOfExpansion(obj)
           for rIdx = 1:obj.numRounds
                for sIdx = 1:obj.numSimulationsPerRound
                    expIn(sIdx,:) = obj.resultStruct(rIdx,sIdx).affectedBeadsRadOfExpension;
                end
                expRad(rIdx,:) = mean(expIn);
            end
        end
        
        function expRad= CalcualteMeanNonAffectedRadOfExpansion(obj)
                  for rIdx = 1:obj.numRounds
                for sIdx = 1:obj.numSimulationsPerRound
                    expIn(sIdx,:) = obj.resultStruct(rIdx,sIdx).nonAffectedBeadsRadOfExpension;
                end
                expRad(rIdx,:) = mean(bIn);
            end
        end
        
    end
    
    methods (Static)
        function ns = NewStruct()
            ns = struct('date',[],...
                        'round',[],...
                        'simulation',[],...
                        'dimension',[],...
                        'numRelaxationSteps',[],...
                        'numRecordingSteps',[],...
                        'numBeamSteps',[],...
                        'numREpairSteps',[],...
                        'numBeads',[],...
                        'bendingConst',[],...
                        'springConst',[],...
                        'openingAngle',[],...
                        'connectedBeads',[],...
                        'connectedBeadsAfterBeam',[],...
                        'numConnectionsLost',[],...
                        'dt',[],...
                        'ROI',struct('rectX',[],'rectY',[],'rectWidth',[],'rectHeight',[]),...
                        'params',[],...
                        'numBeadsIn',[],...
                        'inBeam',[],...
                        'beadsInIndex',[],...
                        'concentricDensity',[],...                                                        
                        'percentDNALoss',[],...
                        'percentHistoneLoss',[],...
                        'affectedBeadsRadOfExpension',[],...
                        'nonAffectedBeadsRadOfExpension',[],...
                        'chainPos',[]);
        end
    end
end
    