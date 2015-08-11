classdef BeamDamageSimulationViewer<handle
    properties
        simulator
        chainPosition 
        affectedMonomers
        radiusOfExpansion
        numMonomersInROI
        connectedMonomers
        connectedMonomersAfterBeam
        connectedMonomersAfterRepair
        concentricDensity
        roi
        params
        handles
        dimension
    end
        
    methods
        function obj = BeamDamageSimulationViewer(resultStruct)
%             obj.simulator        = bSimulationClass;
            obj.params           = resultStruct.params;
            obj.dimension        = obj.params.dimension;
            obj.chainPosition    = resultStruct.chainPosition;
            obj.affectedMonomers = resultStruct.inBeam;
            obj.radiusOfExpansion.nonAffected = resultStruct.nonAffectedBeadsRadOfExpension;
            obj.radiusOfExpansion.affected    = resultStruct.affectedBeadsRadOfExpension;
            obj.roi                           = resultStruct.ROI;
            obj.numMonomersInROI              = resultStruct.numBeadsIn;
            obj.connectedMonomers             = resultStruct.connectedBeads;
            obj.connectedMonomersAfterBeam    = resultStruct.connectedBeadsAfterBeam;%true(size(resultStruct.connectedBeads,1),1);
            if isfield('connectedBeadsAfterRepair',resultStruct) % backward compatibility 
            obj.connectedMonomersAfterRepair  = resultStruct.connectedBeadsAfterRepair;
            end
            % make a list of indices for the difference between the
            % connected before and after beam
            cmInd = false(size(obj.connectedMonomers,1),1);
            for sIdx = 1:size(obj.connectedMonomers,1)
                if any([resultStruct.connectedBeadsAfterBeam(:,1)==obj.connectedMonomers(sIdx,1) & ...
                       resultStruct.connectedBeadsAfterBeam(:,2)==obj.connectedMonomers(sIdx,2)])
                    cmInd(sIdx) = true; % keep the connection after beam 
                end
            end         
            obj.connectedMonomersAfterBeam = cmInd;
            obj.concentricDensity          = resultStruct.concentricDensity;
            if ~isempty(obj.chainPosition)
                obj.InitializeGraphics;
            end
        end
        
        function InitializeGraphics(obj)
            % main figure
            obj.handles.mainFigure    = figure('Units','norm','ToolBar','none');
            cameratoolbar(obj.handles.mainFigure);
            % panels
            obj.handles.mainPanel     = uipanel('Parent',obj.handles.mainFigure,'Units','norm','Position', [0 0.1 0.5 0.9]);
            obj.handles.controlPanel  = uipanel('Parent',obj.handles.mainFigure,'Units','norm','Position', [0 0 1 0.1]);
            obj.handles.mainAxes      = axes('Parent',obj.handles.mainPanel,'Units','norm','Position',[0 0 1 1],'XTick',[],'YTick',[],'ZTick',[],'Color',[0.5  0.72 0.53]);
            box(obj.handles.mainAxes,'on');
            obj.handles.infoPanel     = uipanel('Parent',obj.handles.mainFigure,'Units','norm','Position',[0.5,0.1,0.5,0.9]);
            
            % expansion axes
            obj.handles.expansionAxes = axes('Parent',obj.handles.infoPanel,'Units','norm','Position',[0.1 0.72 0.85 0.25],'FontSize',8); 
            ylabel(obj.handles.expansionAxes,'Radius'), xlabel(obj.handles.expansionAxes,'Step');  
            title(obj.handles.expansionAxes,'Radius of expansion','FontSize',13)
            set(obj.handles.expansionAxes,'FontSize',13);            
            
            % monomer density (num monomers in ROI)
            obj.handles.densityAxes   = axes('Parent',obj.handles.infoPanel,'Units','norm','Position',[0.1 0.38 0.85 0.25],'FontSize',8);
            title(obj.handles.densityAxes,'Number of monomers in ROI','FontSize',13)
            ylabel(obj.handles.densityAxes,'Num. monomers in ROI'); xlabel(obj.handles.densityAxes,'Step')
            set(obj.handles.densityAxes,'FontSize',13);
            
            % concentric density 
            obj.handles.concentricDensityAxes   = axes('Parent',obj.handles.infoPanel,'Units','norm','Position',[0.1 0.05 0.85 0.25],'FontSize',8);
            title(obj.handles.concentricDensityAxes,'Number of monomers in band','FontSize',13)
            ylabel(obj.handles.concentricDensityAxes,'concentric num. monomers in ROI','FontSize',11); xlabel(obj.handles.concentricDensityAxes,'dist')
            set(obj.handles.concentricDensityAxes,'FontSize',13);
            
            daspect(obj.handles.mainAxes,[1 1 1]);
            mx = (obj.chainPosition(:,1,:));
            mx = mx(:);
            my = (obj.chainPosition(:,2,:));
            my = my(:);
            mz = obj.chainPosition(:,3,:);
            mz = mz(:);
            if min(mz)==max(mz)
                zlim = [0 0.00001];
            end
            set(obj.handles.mainAxes,'XLim',[min(mx), max(mx)],'YLim',[min(my), max(my)],'ZLim',[zlim(1) zlim(2)]);
            obj.handles.frameSlider  = uicontrol('style','slider','Units','norm',...
                                                 'Parent',obj.handles.controlPanel,'Position',[0.3 0, 0.5 0.9],...
                                                 'Callback',@obj.SliderMotion,'Max',size(obj.chainPosition,3),'Min',1,'Value',1,...
                                                 'SliderStep',[10 10].*(1/(size(obj.chainPosition,3)-1)));
           % initialize the chain graphics
           obj.handles.chain           = line('XData',obj.chainPosition(:,1,1),'YData',obj.chainPosition(:,2,1),...
               'ZData', obj.chainPosition(:,3,1),'Marker','o','Parent',obj.handles.mainAxes,...
               'MarkerEdgeColor','b');
           % plot extra connectors
            for cIdx = 1:size(obj.connectedMonomers,1)
                obj.handles.connectors(cIdx) = line(...
                    'XData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),1,1),obj.chainPosition(obj.connectedMonomers(cIdx,2),1,1)],...
                    'YData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),2,1),obj.chainPosition(obj.connectedMonomers(cIdx,2),2,1)],...
                    'ZData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),3,1),obj.chainPosition(obj.connectedMonomers(cIdx,2),3,1)],...
                    'Parent',obj.handles.mainAxes,'Color','y');
            end
            
           obj.handles.damagedMonomers = line('XData',obj.chainPosition(obj.affectedMonomers,1,1),...
                                              'YData',obj.chainPosition(obj.affectedMonomers,2,1),...
                                              'ZData',obj.chainPosition(obj.affectedMonomers,3,1),...
                                              'Marker','o','MarkerFaceColor','r','lineStyle','none',...
                                              'Parent',obj.handles.mainAxes,'Visible','off');
            % expansion radius 
            t   = linspace(0,2*pi,20);
            cmz = mean(obj.chainPosition(:,3,1));
            obj.handles.radiusOfExpansion.damaged       = line('XData',cos(t),'YData',sin(t),'ZData',ones(1,numel(t)).*cmz,'Color','r','Visible','off','Parent',obj.handles.mainAxes,'LineWidth',2);
            obj.handles.radiusOfExpansion.nonDamaged    = line('XData',cos(t),'YData',sin(t),'ZData',ones(1,numel(t)).*cmz,'Color','b','Visible','off','Parent',obj.handles.mainAxes,'LineWidth',2);
            obj.handles.roi                             = line('XData',obj.roi*cos(t),'YData',obj.roi*sin(t),'ZData',ones(1,numel(t)).*cmz,'Color','m','Parent',obj.handles.mainAxes,'Visible','off','LineWidth',3);
            obj.handles.expansion.affected              = line('XData',1:numel(obj.radiusOfExpansion.affected),'YData',obj.radiusOfExpansion.affected,'Parent',obj.handles.expansionAxes,'Color','r','DisplayName','damaged');
            obj.handles.expansion.indicator.affected    = line('XData',1,'YData',obj.radiusOfExpansion.affected(1),'marker','o','markerFacecolor','r','Parent',obj.handles.expansionAxes,'HandleVisibility','off');
            obj.handles.expansion.nonAffected           = line('XData',1:numel(obj.radiusOfExpansion.nonAffected),'YData',obj.radiusOfExpansion.nonAffected,'Parent',obj.handles.expansionAxes,'Color','b','Displayname','non damaged');
            obj.handles.expansion.indicator.nonAffected = line('XData',1,'YData',obj.radiusOfExpansion.nonAffected(1),'marker','o','markerFacecolor','b','Parent',obj.handles.expansionAxes,'HandleVisibility','off');
            obj.handles.numMonomersInROI                = line('XData',1:numel(obj.numMonomersInROI),'YData',obj.numMonomersInROI,'Parent',obj.handles.densityAxes);
            obj.handles.numMonomersIndicator            = line('XData',1,'YData',obj.numMonomersInROI(1),'Parent',obj.handles.densityAxes,'marker','o','markerFaceColor','k');
            obj.handles.concentricNumMonomersInROI      = line('XData',1:numel(obj.concentricDensity(1,:)),'YData',obj.concentricDensity(1,:),'Parent',obj.handles.concentricDensityAxes);            
            legend(obj.handles.expansionAxes,'show')
            obj.AddTimeLineToAxes(obj.handles.expansionAxes);
            obj.AddTimeLineToAxes(obj.handles.densityAxes);
            axes(obj.handles.mainAxes)% give focus to main axes
        end
        
        function SliderMotion(obj, sliderHandle,varargin)
            frameNum = round(get(sliderHandle,'Value')); 
            % update chain position
            set(obj.handles.chain,'XData',obj.chainPosition(:,1,frameNum),...
                                  'YData',obj.chainPosition(:,2,frameNum),...
                                  'ZData',obj.chainPosition(:,3,frameNum));
            % update affected monomes                              
            set(obj.handles.damagedMonomers,'XData',obj.chainPosition(obj.affectedMonomers,1,frameNum),...
                                            'YData',obj.chainPosition(obj.affectedMonomers,2,frameNum),...
                                            'ZData',obj.chainPosition(obj.affectedMonomers,3,frameNum));
                                        
            set(obj.handles.expansion.indicator.affected,'XData',frameNum,'YData',obj.radiusOfExpansion.affected(frameNum));
            set(obj.handles.expansion.indicator.nonAffected,'XData',frameNum,'YData',obj.radiusOfExpansion.nonAffected(frameNum));
            set(obj.handles.numMonomersIndicator,'XData',frameNum,'YData',obj.numMonomersInROI(frameNum));
            set(obj.handles.concentricNumMonomersInROI,'XData',1:numel(obj.concentricDensity(frameNum,:)),'YData',obj.concentricDensity(frameNum,:));
            
            
            if frameNum>obj.params.numRecordingSteps && frameNum<=(obj.params.numRecordingSteps+obj.params.numBeamSteps) % beam
                set(obj.handles.damagedMonomers,'Visible','on','markerFaceColor','r')
                set(obj.handles.frameSlider,'BackgroundColor','r')
                cm = mean(obj.chainPosition(obj.affectedMonomers,:,frameNum),1);
                % plot radius of expansion 
                t= linspace(0,2*pi,20);
                set(obj.handles.radiusOfExpansion.damaged,'XData',cm(1)+obj.radiusOfExpansion.affected(frameNum)*cos(t),...
                    'YData',cm(2)+obj.radiusOfExpansion.affected(frameNum)*sin(t),'Visible','on');
                set(obj.handles.radiusOfExpansion.nonDamaged,'XData',cm(1)+obj.radiusOfExpansion.nonAffected(frameNum)*cos(t),...
                    'YData',cm(2)+obj.radiusOfExpansion.nonAffected(frameNum)*sin(t),'Visible','on');
                
                set(obj.handles.roi,'XData', cm(1)+obj.roi*cos(t),'Ydata',cm(2)+obj.roi*sin(t),'Visible','on');
                
             % plot extra connectors
            for cIdx = 1:size(obj.connectedMonomers,1)
                if obj.connectedMonomersAfterBeam(cIdx);
                vState = 'on';
                else
                    vState = 'off';
                end
                set(obj.handles.connectors(cIdx),...
                    'XData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),1,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),1,frameNum)],...
                    'YData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),2,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),2,frameNum)],...
                    'ZData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),3,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),3,frameNum)],...
                    'Visible',vState);
            end
                
                % Repair
            elseif frameNum>(obj.params.numRecordingSteps+obj.params.numBeamSteps)&& frameNum<=(obj.params.numRecordingSteps+obj.params.numBeamSteps+obj.params.numRepairSteps)% repair
                set(obj.handles.frameSlider,'BackgroundColor','g')
                set(obj.handles.damagedMonomers,'MarkerFacecolor','g')
                             % plot extra connectors
            for cIdx = 1:size(obj.connectedMonomers,1)
                if obj.connectedMonomersAfterRepair(cIdx);
                    vState = 'on';
                else
                    vState = 'off';
                end
                set(obj.handles.connectors(cIdx),...
                    'XData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),1,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),1,frameNum)],...
                    'YData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),2,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),2,frameNum)],...
                    'ZData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),3,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),3,frameNum)],...
                    'Visible',vState);
            end
            
            else % recording
                set(obj.handles.damagedMonomers,'Visible','off')
                set(obj.handles.frameSlider,'BackgroundColor',[0.9400 0.9400 0.9400]);
                             % plot extra connectors
            for cIdx = 1:size(obj.connectedMonomers,1)                
                set(obj.handles.connectors(cIdx),...
                    'XData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),1,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),1,frameNum)],...
                    'YData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),2,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),2,frameNum)],...
                    'ZData',[obj.chainPosition(obj.connectedMonomers(cIdx,1),3,frameNum),obj.chainPosition(obj.connectedMonomers(cIdx,2),3,frameNum)],...
                    'Visible','on');
            end
            end
        end
        
        function AddTimeLineToAxes(obj,axesHandle)
            % plot timeline 
            ylim      = get(axesHandle,'YLim');
            fontSize  = 12;
            textYpos  = ylim(1)+(ylim(2)-ylim(1))*0.05;
            text(obj.params.numRecordingSteps,textYpos,'Beam shot','Color','r','LineStyle','-.','FontSize',fontSize,'HandleVisibility','off','Parent',axesHandle)
            line('XData',obj.params.numRecordingSteps*[1 1], 'YData',ylim,'Color','r','Parent',axesHandle,'HandleVisibility','off');% beam start
            line('XData',(obj.params.numRecordingSteps+obj.params.numBeamSteps)*[1 1],'YData',ylim,'Color','g','Parent',axesHandle,'HandleVisibility','off')
            text(obj.params.numRecordingSteps+obj.params.numBeamSteps,textYpos,'Repair','Color','g','LineStyle','-.','FontSize',fontSize,'HandleVisibility','off','Parent',axesHandle)
            
        end
        
        function ExportGIF(obj)
            % set the slider at step(1)
            res = 15;
            set(obj.handles.frameSlider,'Value',1);
            obj.SliderMotion(obj.handles.frameSlider)
            f = getframe;
            [im,map] = rgb2ind(f.cdata,256,'nodither');                        
            im(1,1,1,numel(1:res:get(obj.handles.frameSlider,'Max'))) = 0;
            for k = 1:res:get(obj.handles.frameSlider,'Max')
                set(obj.handles.frameSlider,'Value',k);
                obj.SliderMotion(obj.handles.frameSlider)
                 f = getframe;
                 im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
            end
           imwrite(im,map,'simulation.gif','DelayTime',0,'LoopCount',inf) %g443800     
        end
            
        end
 end
    
