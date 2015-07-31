classdef BeamDamageResultsAnalysis<handle
    % This class analyzes the results from BeamDamageSimulation
    % it take as an input the class BeamDamageSimulation and allows
    % presentation of the resultStruct in the results property of
    % BeamDamageSimulation class
    properties
        results
        numRounds       % the number of simulationi rounds 
        numSimulations  % the number of simulations per round 
        numSteps
        numRelaxationSteps
        numRecordingSteps
        numBeamSteps 
        numRepairSteps
        recordingTime
        beamTime
        repairTime
        dt
        meanNumBeadsIn
        meanNumBeadsBeforeBeam
        meanNumBeadsAfterBeam
        meanNumBeadsAfterRepair
        
        meanPercentBeadsBeforeBeam
        meanPercentBeadsAfterBeam
        meanPercentBeadsAfterRepair
        
        meanRadiusOfExpansionAffected
        meanRadiusOfExpansionNonAffected
    end
    
    methods
        function obj = BeamDamageResultsAnalysis(beamDamageResults)
            obj.results            = beamDamageResults.resultStruct;                        
            obj.numRounds          = size(obj.results,1);
            obj.numSimulations     = size(obj.results,2);
            obj.dt                 = obj.results(1,1).dt;
            obj.numRelaxationSteps = obj.results(1,1).numRelaxationSteps;
            obj.numRecordingSteps  = obj.results(1,1).numRecordingSteps;
            obj.numBeamSteps       = obj.results(1,1).numBeamSteps;
            obj.numRepairSteps     = obj.results(1,1).numRepairSteps;
            obj.recordingTime      = (obj.numRelaxationSteps-obj.numRelaxationSteps)*obj.dt;
            obj.beamTime           = (obj.numRecordingSteps)*obj.dt;
            obj.repairTime         = (obj.numRecordingSteps+obj.numBeamSteps)*obj.dt;
        end
        
        function MeanNumBeadsIn(obj,plotFig)
            % Display the mean number of beads in the ROI
            % plotFig is a flag for plotting, leave empty or set to true
            % for plotting 
            obj.meanNumBeadsIn = zeros(obj.numRounds,numel(obj.results(1,1).numBeadsIn));
            for rIdx = 1:obj.numRounds
                nbIn = zeros(obj.numSimulations,numel(obj.results(1,1).numBeadsIn));
                for sIdx = 1:obj.numSimulations
                    nbIn(sIdx,:) = obj.results(rIdx,sIdx).numBeadsIn;
                end
                 obj.meanNumBeadsIn(rIdx,:) = mean(nbIn,1);                
            end
            
            if ~exist('plotFig','var')
                plotFig = true;
            end
                        
            if plotFig
            stepsToPlot = 1:(obj.numRecordingSteps+obj.numBeamSteps+obj.numRepairSteps);
            
            f= figure;             
            a = axes('Parent',f);
            title(a,'Number of Monomers in the ROI','FontSize',20);
            ylabel(a,'Mean Num. monomers in ROI')
            lineWidth =  4;                                    
           
            [cr,cg,cb] = obj.GetColors(obj.numRounds);
                        
            for rIdx = 1:obj.numRounds
              line('XData',stepsToPlot*obj.dt,'YData',  obj.meanNumBeadsIn(rIdx,:),'displayName',sprintf('%s',['experiment ',num2str(rIdx)]),...
                  'Color',[cr(rIdx),cg(rIdx),cb(rIdx)],'LineWidth',lineWidth,'Parent',a) 
            end
            
            %  lines indicating recording time, beam time  and recoverty
            % time 
            obj.PlotTimeline(a)
            
            end
        end
        
        function MeanLostAndRecovery(obj,plotFig)
           % Calculate the mean bead in if not already calculated
            if  isempty(obj.meanNumBeadsIn)
               obj.MeanNumBeadsIn(false)
            end
            if ~exist('plotFig','var')
                plotFig = true;
            end
            % for each experiment, calculate the mean number of monomers
           % before beam, mean number after beam, and mean number at repair
           obj.meanNumBeadsBeforeBeam  = zeros(1,obj.numRounds);
           obj.meanNumBeadsAfterBeam   = zeros(1,obj.numRounds);
           obj.meanNumBeadsAfterRepair = zeros(1,obj.numRounds);           
            for rIdx = 1:obj.numRounds
                obj.meanNumBeadsBeforeBeam (rIdx)  = mean(obj.meanNumBeadsIn(rIdx,obj.numRecordingSteps));
                obj.meanNumBeadsAfterBeam   (rIdx) = mean(obj.meanNumBeadsIn(rIdx,obj.numRecordingSteps+1:(obj.numRecordingSteps+obj.numBeamSteps)));
                obj.meanNumBeadsAfterRepair (rIdx) = mean(obj.meanNumBeadsIn(rIdx,(obj.numRecordingSteps+obj.numBeamSteps)+1:(obj.numRecordingSteps+obj.numBeamSteps+obj.numRepairSteps)));                
            end
           
            if plotFig
                figure, 
                if obj.numRounds==1
                   barVals =  [obj.meanNumBeadsBeforeBeam nan;obj.meanNumBeadsAfterBeam nan;obj.meanNumBeadsAfterRepair nan]';
                else 
                   barVals = [obj.meanNumBeadsBeforeBeam;obj.meanNumBeadsAfterBeam;obj.meanNumBeadsAfterRepair]';
                end
                b=bar(barVals);
                
                set(b(1),'DisplayName','num. monomers in ROI before Beam')
                set(b(2),'DisplayName','num. monomers in ROI during repair')
                set(b(3),'DisplayName','num. monomers in ROI after repair')
            end
            legend('show')
           title('Mean number of monomers in ROI at 3 stages')
           xlabel('Experiment number')
           ylabel('Mean number of Monomers in ROI')
           set(gca,'FontSize',20);
        end
        
        function MeanLostAndRecoveryPrecentage(obj,plotFig)
           % Calculate the mean bead in if not already calculated
            if  isempty(obj.meanNumBeadsBeforeBeam)
               obj.MeanLostAndRecovery(false)
            end
            
            if ~exist('plotFig','var')
                plotFig = true;
            end
            
           % For each experiment, calculate the mean number of monomers
           % before beam, mean number after beam, and mean number at repair                  
            for rIdx = 1:obj.numRounds
                obj.meanPercentBeadsBeforeBeam(rIdx)  = 100;
                obj.meanPercentBeadsAfterBeam(rIdx)   = 100*obj.meanNumBeadsAfterBeam(rIdx)/obj.meanNumBeadsBeforeBeam(rIdx);
                obj.meanPercentBeadsAfterRepair(rIdx) = 100*obj.meanNumBeadsAfterRepair(rIdx)/obj.meanNumBeadsBeforeBeam(rIdx);                
            end
           
            if plotFig
                figure, 
                if obj.numRounds==1
                   barVals =  [obj.meanPercentBeadsBeforeBeam nan;obj.meanPercentBeadsAfterBeam nan;obj.meanPercentBeadsAfterRepair nan]';
                else 
                   barVals = [obj.meanPercentBeadsBeforeBeam;obj.meanPercentBeadsAfterBeam;obj.meanPercentBeadsAfterRepair]';
                end
                b=bar(barVals);
                
                set(b(1),'DisplayName','percent monomers in ROI before Beam')
                set(b(2),'DisplayName','percent monomers in ROI during repair')
                set(b(3),'DisplayName','percent monomers in ROI after repair')
            end
            legend('show')
           title('Percent of monomers in ROI at 3 stages')
           xlabel('Experiment number')
           ylabel('Percent of Monomers in ROI')
           set(gca,'FontSize',20);
        end
        
        function ConcentricDensity(obj,plotFlag)% unfinished
        end
        
        function RadiusOfExpansion(obj,plotFig)
            % Calculate the mean radius of expansion for affected and non
            % affected monomers 
           if ~exist('plotFig','var')
               plotFig = true;
           end
               obj.meanRadiusOfExpansionAffected = zeros(obj.numRounds,numel(obj.results(1,1).affectedBeadsRadOfExpension));
               obj.meanRadiusOfExpansionNonAffected = zeros(obj.numRounds,numel(obj.results(1,1).nonAffectedBeadsRadOfExpension));
            for rIdx = 1:obj.numRounds
                ar  = zeros(obj.numSimulations,numel(obj.results(1,1).affectedBeadsRadOfExpension));
                nar = zeros(obj.numSimulations,numel(obj.results(1,1).nonAffectedBeadsRadOfExpension));
                for sIdx = 1:obj.numSimulations
                    ar(sIdx,:)  = obj.results(rIdx,sIdx).affectedBeadsRadOfExpension;
                    nar(sIdx,:) = obj.results(rIdx,sIdx).nonAffectedBeadsRadOfExpension;
                end
                 obj.meanRadiusOfExpansionAffected(rIdx,:) = mean(ar,1);     
                 obj.meanRadiusOfExpansionNonAffected(rIdx,:) = mean(nar,1);
            end
            
            if plotFig
                
            stepsToPlot = 1:(obj.numRecordingSteps+obj.numBeamSteps+obj.numRepairSteps);
            
            % plot affected beads radius of expansion
            f  = figure;            
            ax = axes('Parent',f);
            title(ax,'Mean Radius Of Expansion- affected monomers','FontSize',20)
            lineWidth =  4;                                               
            [cr,cg,cb] = obj.GetColors(obj.numRounds);
                        
            for rIdx = 1:obj.numRounds
              line('XData',stepsToPlot*obj.dt,'YData',  obj.meanRadiusOfExpansionAffected(rIdx,:),'displayName',sprintf('%s',['experiment ',num2str(rIdx)]),...
                  'Color',[cr(rIdx),cg(rIdx),cb(rIdx)],'LineWidth',lineWidth,'Parent',ax) 
            end
            %  lines indicating recording time, beam time  and recoverty
            % time 
            obj.PlotTimeline(ax);
            ylabel(ax,'Radius')
            set(ax,'FontSize',30,'LineWidth',4);
            legend(ax,'show')
            
            % Plot non-affected beads radius of expansion
            f1 = figure;                       
            ax1 = axes('Parent',f1);
            title(ax1,'Mean Radius Of Expansion-non affected monomers','FontSize',20);
            lineWidth =  4;                                               
            [cr,cg,cb] = obj.GetColors(obj.numRounds);
                        
            for rIdx = 1:obj.numRounds
              line('XData',stepsToPlot*obj.dt,'YData',  obj.meanRadiusOfExpansionNonAffected(rIdx,:),'displayName',sprintf('%s',['experiment ',num2str(rIdx)]),...
                  'Color',[cr(rIdx),cg(rIdx),cb(rIdx)],'LineWidth',lineWidth,'Parent',ax1) 
            end
            %  lines indicating recording time, beam time  and recoverty
            % time 
            obj.PlotTimeline(ax1);            
            ylabel(ax1,'Radius')
            set(ax1,'FontSize',30,'LineWidth',4);
            legend(ax1,'show')
            end            
        end
        
        function PlotTimeline(obj,axisHandle)
            % Show the line for recording time, beam time, and repair time 
            hold(axisHandle,'on');
            ylim      = get(axisHandle,'YLim');
            textYpos  = ylim(1)+(ylim(2)-ylim(1))*0.05;
            fontSize  = 20;           
            line('XData',[obj.recordingTime, obj.recordingTime], 'YData',ylim,'Color','b','HandleVisibility','off','Parent',axisHandle,'DisplayName','Recording')
            % add text 
            text(obj.recordingTime+obj.dt,textYpos,'Recording','Color','b','LineStyle','-.','FontSize',fontSize,'HandleVisibility','off','Parent',axisHandle)
            % beam time 
            line('XData',[obj.beamTime, obj.beamTime], 'YData',ylim,'Color','r','HandleVisibility','off','Parent',axisHandle,'DisplayName','Beam Shot')
            text(obj.beamTime+obj.dt,textYpos,'Beam Shot','Color','r','LineStyle','--','FontSize',fontSize,'HandleVisibility','off','Parent',axisHandle)

            line('XData',[obj.repairTime, obj.repairTime], 'YData',ylim,'Color','g','HandleVisibility','off','Parent',axisHandle,'DisplayName','Repair')
            text(obj.repairTime+obj.dt,textYpos,'Repair','Color','g','LineStyle','--','FontSize',fontSize,'HandleVisibility','off','Parent',axisHandle)            
            xlabel(axisHandle,'Time [sec]');            
            set(axisHandle,'FontSize',20)
        end
    end
    
    methods (Static)
        
        function [cr,cg,cb] = GetColors(numColors)
             % sample the color channel
            c  = colormap; % the current colormap 
            cr = interp1(1:size(c,1),c(:,1)',linspace(1,size(c(:,1),1),numColors));
            cg = interp1(1:size(c,1),c(:,2)',linspace(1,size(c(:,2),1), numColors));
            cb = interp1(1:size(c,1),c(:,3)',linspace(1,size(c(:,3),1), numColors));
        end        

    end
end