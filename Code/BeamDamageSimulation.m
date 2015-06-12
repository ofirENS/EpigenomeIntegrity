classdef BeamDamageSimulation<handle %[UNFINISHED]
    properties
        params
        handles
        particlePosition
    end
    
    methods
        
        function obj = BeamDamageSimulation(params)
            obj.params = params;
            if ~obj.params.loadConfiguration
                obj.InitializeChain
                obj.Initialize
            else
                obj.RunRelaxationSteps;
                obj.RunRecordingSteps;
                obj.RunBeamSteps;
            end
        end
        
        function Run(obj)
            
        end
        
        function InitializeChain(obj)
            % Initial position for the chain 
            obj.particlePosition = cumsum(sqrt(obj.params.forceParams.diffusionConst*obj.params.dt)*randn(obj.params.numParticles,obj.params.dimension));
        end
        
        function Initialize(obj)
                                              
            % -Relaxation time for the Rouse chain defined by the longest relaxation time--
            % relaxation times - 300 beads ~= 1000 steps
            %                    400 beads ~= 2000 steps
            %                    500 beads ~= 3000 steps
            
            % (numBeads*b)^2 / (3*D*pi^2) %[Doi &  Edwards p.96 eq. 4.37)
            % using: b  = sqrt(3)~=1.7
            %        dt = 0.01
            %        D  = 1; diffusion const.
            % (500*sqrt(3))^2 /(3*pi^2 * 1)
            
            
        end
        
        function RunRelaxationSteps(obj)
            for sIdx = 1:obj.params.numRelaxationSteps
                % get particle distance
             particleDist  = ForceParams.GetParticleDistance(obj.particlePosition);
             % get spring force
             springForce  = ForceManager.GetSpringForce(obj.params.forceParams.springForce,obj.particlePosition,...
                 particleDist,obj.params.forceParams.springConst,connectivityMap,obj.params.forceParams.minParticleEqDistance,...
                 obj.params.fixedParticleNum);
             % get diffusion force
             diffusionForce = ForceManager.GetDiffusionForce(obj.params.forceParams.diffusionForce,obj.particlePosition,...
                 obj.params.forceParams.diffusionConst,obj.params.dt,obj.params.fixedParticleNum);
             % set new pariclePosition
             obj.particlePosition = obj.particlePosition + springForce*obj.params.dt +diffusionForce;
            end
        end
        
        function RunBeamSteps(obj)
            % start beam
            chainPos = r.objectManager.GetPosition(1);
            chainPos = chainPos{1};
            % move the beam to the chain's center of mass
            UpdateBeamPosition(chainPos,r,2);
            if r.params.simulator.showSimulation
                UpdateBeamGraphics(r,2)
            end
            
            % Create DNA damages in the ray
            inBeam        = r.handles.classes.domain.InDomain(chainPos,2);
            inBeamInds    = find(inBeam);
            % choose a fraction of inBeam to exclude
            frInBeam      = randperm(sum(inBeam));
            fracToExclude = 0; % fraction of damages to induce on the edges in the ray
            inBeam(inBeamInds(frInBeam(1:round(numel(frInBeam)*fracToExclude))))= false;
            
            r.runSimulation = true;
            
            connectivityMat        = r.objectManager.GetConnectivityMapAsOne(1);
            bendingElasticityConst = 1/(r.params.simulator.dt);
            
            disp ('Beam on')
            
            while all([r.simulationData.step<(r.simulationData.step+numBeamSteps),r.runSimulation])
                % Advance one simulation step
                [r,h,chainPos] = Step(r,h);
                [chainPos]     = ApplyDamageEffect(chainPos,inBeam,connectivityMat,bendingElasticityConst,r.params.simulator.dt);
                r.objectManager.DealCurrentPosition(1,chainPos)
                
                % Update projectionPlane position according to the cm of the chain
                [rectX,rectY] = UpdateProjectionPlanePositionByCM(chainPos,rectWidth,rectHeight);
                
                % Calculate histone and DNA densities in the ROI
                [histoneDensity, dnaDensity] = CalculateDensitiesInROI(chainPos,h.curPos,rectX,rectY,rectWidth,rectHeight,histoneParams.numHistones);
                if r.params.simulator.showSimulation
                    %         gyrationsphereData   = r.simulationGraphics.handles.graphical.domain(3);% update gyration sphere center to the chain's center of mass
                    %         UpdateGyrationSphere(r,3,gyrationsphereData,chainPos);
                    UpdateGraphics(r.simulationData.step,r.params.simulator.dt,h.curPos,chainPos,histoneDensity,dnaDensity,...
                        rectX,rectY,rectWidth, rectHeight,dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
                end
            end
        end
        
        function RunRecordingSteps(obj)
            % Start recording densities with no beam effect
            for sIdx = 1:obj.params.numRecordingSteps
                % advance one step
                [r,h,chainPos] = obj.Step(r,h);
                % Update projectionPlane position according to the cm of the chain
                [rectX,rectY] = obj.UpdateProjectionPlanePositionByCM(chainPos,rectWidth,rectHeight);
                
                % Calculate histone and DNA densities in the ROI
                [histoneDensity, dnaDensity] = obj.CalculateDensitiesInROI(chainPos,h.curPos,rectX,rectY,rectWidth,rectHeight,histoneParams.numHistones);
                
                % update graphics
                if r.params.simulator.showSimulation
                    UpdateGraphics(r.simulationData.step,r.params.simulator.dt,h.curPos,chainPos,histoneDensity,dnaDensity,rectX,rectY,rectWidth, rectHeight,...
                        dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
                end
            end
        end
        
        function [r,h,chainPos] = Step(r,h)
            r.Step;
            [~,chainPos] = r.objectManager.GetMembersPosition(1);
            chainPos     = chainPos{1};
            
            % Move the histones
            h.Step(chainPos,r.params.simulator.dt);% update current position
        end
        
        function [chainPos] = ApplyDamageEffect(chainPos,inBeam,connectivityMat,bendingElasticityConst,dt)
            % Apply forces on the chain falling in the beam
            bForce             = ForceManager.GetBendingElasticityForce(true,chainPos,connectivityMat,bendingElasticityConst ,[]);
            bForce(~inBeam,:)  = 0;
            % Update affected edges of the chain
            chainPos(inBeam,:) = chainPos(inBeam,:) + bForce(inBeam,:)*dt;
        end
        
        function InitializeGraphics(obj)
            %     mAxes = r.simulationGraphics.handles.graphical.mainAxes;
            
            % Create graphics
            obj.handles.mainFig       = figure('Units','norm');
            obj.handles.densityFigure = figure('Units','norm');
            obj.handles.projFig       = figure('Units','norm');
            obj.handles.mainAxes      = axes('Parent',obj.handles.mainFig,'Units','norm','Color','k','XLim',.5*obj.params.gyrationRadius.*[-1 1],...
                'YLim',.5*obj.params.gyrationRadius.*[-1 1],'ZLim',.5*obj.params.gyrationRadius.*[-1 1],'FontSize',25);
            obj.handles.densityAxes   = axes('Parent',obj.handles.densityFigure,'Units','norm','FontSize',25,'YLim',[0 obj.params.numParticles]);
            obj.handles.projAxes      = axes('Parent',obj.handles.projFig,'Units','norm','FontSize',25,'XLim',get(obj.handles.mainAxes,'XLim'),'YLim',get(obj.handles.mainAxes,'YLim'));
            
            xlabel(obj.handles.densityAxes,'Time'); ylabel(obj.handles.densityAxes,'Num. Beads In ROI')
            
            cameratoolbar(obj.handles.mainFig);
            daspect(obj.handles.mainAxes,[1 1 1]);
            daspect(obj.handles.projAxes,[1 1 1]);
            
            obj.handles.particleHandle = line('XData',obj.handles.particlePosition(:,1),...
                'Ydata',obj.handles.particlePosition(:,2),...
                'Zdata',obj.handles.particlePosition(:,3),...
                'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
                'LineStyle','-','Color','w','LineWidth',1,'Parent',obj.handles.mainAxes);
            
            obj.handles.projBeadHandle = line('XData',obj.particlePosition(:,1),...
                'Ydata',obj.particlePosition(:,2),...
                'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',6,...
                'LineStyle','-','Color','k','LineWidth',1,'Parent',projAxes);
            
            
            obj.handles.affectedBeadsHandle = line('XData',NaN,...
                'Ydata',NaN,...
                'Zdata',NaN,...
                'Marker','o','MarkerFaceColor','g','MarkerEdgeColor','g',...
                'LineStyle','none','Color','w','LineWidth',1,'Parent',obj.handles.mainAxes);
            obj.handles.affectedBeadsProjHandle   = line('XData',obj.particlePosition(affectedBeads,1),...
                'Ydata',obj.particlePosition(affectedBeads,2),...
                'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r',...
                'LineStyle','none','Color','r','LineWidth',1,'Parent',obj.handles.projAxes,'Visible','off');
            
            obj.handles.numBeadsInHandle    = line('XData',NaN,'YData',NaN,'Color','b','Parent',obj.handles.densityAxes,'LineWidth',4);
            
            cm    = mean(obj.particlePosition,1); % center of mass
            
            % create the beam
            numCirclesInBeam = 40;
            circlezPos       = linspace(-2*obj.params.gyrationRadius,2*obj.params.gyrationRadius,numCirclesInBeam);
            theta    = 0:0.1:2*pi;
            cosTheta = cos(theta);
            sinTheta = sin(theta);
            for cIdx = 1:numCirclesInBeam;
                obj.handles.beamCircleHandle(cIdx)= line('XData',obj.params.beamRad*cosTheta+cm(1),...
                                                         'YData',obj.params.beamRad*sinTheta+cm(2),...
                                                         'ZData',ones(1,numel(theta))*circlezPos(cIdx),...
                    'Parent',obj.handles.mainAxes,'Color','y','Visible','off');
            end
            
            
            % projection plane in 3D
            obj.handles.roiHandle           = patch([cm(1)-obj.params.roiRadius,cm(1)+obj.params.roiRadius,cm(1)+obj.params.roiRadius,cm(1)-obj.params.roiRadius],...
                [cm(2)-obj.params.roiRadius,cm(2)-obj.params.roiRadius,cm(2)+obj.params.roiRadius,cm(2)+obj.params.roiRadius],'y',...
                'Parent',obj.handles.mainAxes,'FaceAlpha',0.5);
            obj.handles.projRoiHandle       = patch([cm(1)-obj.params.roiRadius,cm(1)+obj.params.roiRadius,cm(1)+obj.params.roiRadius,cm(1)-obj.params.roiRadius],...
                [cm(2)-obj.params.roiRadius,cm(2)-obj.params.roiRadius,cm(2)+obj.params.roiRadius,cm(2)+obj.params.roiRadius],'y',...
                'Parent',obj.handles.projAxes,'FaceAlpha',0.5);
        end
        
        function [rectX,rectY]= UpdateProjectionPlanePositionByCM(chainPos,rectWidth, rectHeight)
            % set x and y to be at the cm
            
            rectX  = mean(chainPos(:,1));
            rectY  = mean(chainPos(:,2));
            rectX  = rectX- rectWidth/2;
            rectY  = rectY- rectHeight/2;
        end
        
        function UpdateGraphics(step,dt,histPosition,chainPos,histoneDensity,dnaDensity,rectX,rectY,rectWidth, rectHeight,dnaDensityHandle,histoneDensityHandle,histHandle,projPlane2D,projPlane3D,pHistHandle,pPolyHandle)
            
            UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,dnaDensityHandle,histoneDensityHandle)
            UpdateHistoneGraphics(histHandle,histPosition)
            UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
            UpdateProjectedHistoneGraphics(pHistHandle,histPosition);
            UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
            drawnow
        end
        
        function UpdateHistoneGraphics(histHandle,curPos)
            set(histHandle,'XData',curPos(:,1),'YData',curPos(:,2),'ZData',curPos(:,3));
        end
        
        function UpdateProjectedHistoneGraphics(pHistHandle,curPos)
            set(pHistHandle,'XData',curPos(:,1),'YData',curPos(:,2));
        end
        
        function UpdateProjectedPolymerGraphics(pPolyHandle,chainPos)
            set(pPolyHandle,'XData',chainPos(:,1),'YData', chainPos(:,2))
        end
        
        function UpdateDensityGraphics(step,dt,histoneDensity,dnaDensity,dnaDensityHandle,histoneDensityHandle)
            dnaDensityDataX     = get(dnaDensityHandle,'XData');
            dnaDensityDataY     = get(dnaDensityHandle,'YData');
            histoneDensityDataX = get(histoneDensityHandle,'XData');
            histoneDensityDataY = get(histoneDensityHandle,'YData');
            
            set(dnaDensityHandle,'XData',[dnaDensityDataX, step*dt], 'YData',[dnaDensityDataY,dnaDensity]);
            set(histoneDensityHandle,'XData',[histoneDensityDataX, step*dt], 'YData',[histoneDensityDataY,histoneDensity]);
            %        line('XData',step*dt,'YData',histoneDensity,'Marker','o','markerFaceColor','y','Parent',dAxes)
            %        line('XData',step*dt,'YData',dnaDensity,'Marker','o','Parent',dAxes,'markerFaceColor','b')
        end
        
        function UpdateProjectionPlaneGraphics(rectX,rectY, rectWidth, rectHeight,projPlane2D, projPlane3D)
            % projPlane3d,projPlane2D are graphical handles
            set(projPlane3D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
            set(projPlane2D,'XData',[rectX, (rectX+rectWidth), (rectX+rectWidth), rectX],...
                'YData',[rectY, rectY, (rectY+rectHeight), (rectY+rectHeight)]);
        end
        
        function UpdateGyrationSphere(simulationFrameworkHandle,domainNumber, gyrationSphereData,chainPos)
            % update the sphere graphics
            ccm  = mean(chainPos,1);% chain center of mass
            xPos = gyrationSphereData.points.x;
            yPos = gyrationSphereData.points.y;
            zPos = gyrationSphereData.points.z;
            
            xPos = xPos -xPos(1,1)+ccm(1);
            yPos = yPos -yPos(1,1)+ccm(2);
            zPos = zPos -zPos((size(zPos,1)-1)/2,1)+ccm(3);
            
            set(gyrationSphereData.mesh,'XData',xPos,...
                'YData',yPos,...
                'ZData',zPos);
            
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.x = xPos;
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.y = yPos;
            simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).points.z = zPos;
        end
        
        function UpdateBeamPosition(chainPos,simulationFrameworkHandle,domainNumber)
            % Change the parametr such that the inDomainfunction will work properly
            cm = mean(chainPos,1);
            simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainCenter = cm;
        end
        
        function UpdateBeamGraphics(simulationFrameworkHandle,domainNumber)
            domainCenter = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainCenter;
            %     domainRad    = simulationFrameworkHandle.handles.classes.domain.params(domainNumber).domainWidth;
            beamHandle   = simulationFrameworkHandle.simulationGraphics.handles.graphical.domain(domainNumber).mesh;
            xPos         = get(beamHandle,'XData');
            yPos         = get(beamHandle,'YData');
            set(beamHandle,'XData',xPos+domainCenter(1),...
                'YData',yPos+domainCenter(2));
            
        end
        
        function [histoneDensity, dnaDensity] = CalculateDensitiesInROI(chainPos,histonePos,rectX,rectY, rectWidth, rectHeight,numHistones)
            % calculate histone density
            histoneDensity =  sum((histonePos(:,1)<=(rectX+rectWidth) & histonePos(:,1)>=rectX &...
                histonePos(:,2)<=(rectY+rectHeight) & histonePos(:,2)>=rectY))/numHistones;
            
            % calculate DNA density
            [dnaLengthIn,totalDNALength] = PolygonLengthInRoi(chainPos(:,1:2),rectX,rectY,rectWidth, rectHeight);
            dnaDensity = dnaLengthIn./totalDNALength;
            
            %         dnaDensity    = sum((chainPos(:,1)<=(rectX+rectWidth) & chainPos(:,1)>=rectX &...
            %                              chainPos(:,2)<=(rectY+rectHeight) & chainPos(:,2)>=rectY));
        end
        
        function simFramework = LoadConfiguration(loadRelaxationConfiguration)
            
            if loadRelaxationConfiguration
                [fName,fPath] = uigetfile('*.mat');
                load(fullfile(fPath,fName)); % simulation framework, saved as r
                f = whos;
                cInd = strcmpi({f.class},'RouseSimulatorFramework');
                if ~any(cInd)% class indicator
                    error('this is not a file of class rouseSimulatorFramework')
                end
                simFramework = eval(f(cInd).name);
            end
        end
        
        function SaveConfiguration(r,saveConfiguration)
            % save RouseSimulatorFramework class
            if saveConfiguration
                uisave('r')
            end
            
        end
    end
end