function [p]=FETIpost(p)
    NodalPos = [];
    fsPost = [];
    xlim = [inf -inf];
    ylim = [inf -inf];
    %PlotMargin = p.elHeight;
    
    
    
    %% compute initial nodal positions
    NodalPos0 = cell(1,p.Ns);
    if p.nonconforming==1
        for s=1:p.Ns
            for r = 1:(p.Nely(s)+1) % r: row (substructure-local)
                for c = 1:(p.Nelx(s)+1) % c: column (substructure-local)
                    i = (r-1)*(p.Nelx(s)+1)+c; % i: node nr
                    NodalPos0{s}(2*i-1:2*i,1)=[p.poss(1,s)+p.elHeight(s)*(c-1);p.poss(2,s)+p.elHeight(s)*(r-1)];
                end
            end
        end
    else
        for iNsy = 1:p.Nsy
            for iNsx = 1:p.Nsx
                s = (iNsy-1)*p.Nsx + iNsx;
                SubBaseX = (iNsx-1)*(p.Nelx(s)*p.elHeight);
                SubBaseY = (iNsy-1)*(p.Nely(s)*p.elHeight);
                for r = 1:(p.Nely(s)+1) % r: row (substructure-local)
                    for c = 1:(p.Nelx(s)+1) % c: column (substructure-local)
                        i = (r-1)*(p.Nelx(s)+1)+c; % i: node nr
                        % NodalPos0{s}(1,2*nodenr-1) = x0
                        NodalPos0{s}(2*i-1:2*i,1) = [SubBaseX + p.elHeight*c; SubBaseY + p.elHeight*r];
                    end
                end
            end
        end
    end
    
    
    %% Dynamic
    if strcmp(p.mode,'dynamic')
        
        tPlot = 0:p.Plotstep:p.End;
        PlotSteps = length(tPlot);        

        tiPlot(PlotSteps) = 0;
        
        if strcmp(p.Solver,'FULLsolver')
            u_direct_loc=p.L_man*p.uFull;
            %disp('u_direct_loc:')
            %disp(size(u_direct_loc))
            fsPost_glob = p.L_man*p.Full.f;
        end
        fsPost = cell(1,PlotSteps);       
            
        NodalPos = cell(1,PlotSteps);
        tracking = zeros(2,PlotSteps);
        
        i = 0;
        for timePlotSet = tPlot
            if strcmp(p.Solver,'FULLsolver')
                end_u=0;
            end
            i = i+1;
            tiPlot(i) = max(1,floor(timePlotSet/p.Step));
            for s = 1:p.Ns
                if p.nonconforming==1
                    PlotMargin=p.elHeight(s);%p.sizes(2,s)/p.Nely(s);
                else
                    PlotMargin=p.elHeight;
                end
                if strcmp(p.Solver,'FULLsolver')
                    start_u=1+end_u;
                    end_u=end_u+size(p.Bs{s},2);
                    uDeformation = u_direct_loc(start_u:end_u,tiPlot(i));
                    %disp('uDeformation:')
                    %disp(uDeformation)
                    %disp(p.L{s}'*uDeformation)
                    fsPost{i}{s} = fsPost_glob(start_u:end_u,i);
                    NodalPos{i}{s} = NodalPos0{s} + p.L{s}'*uDeformation;
                else                
                    % blow up external forces to full coordinates
                    dofEnd = sum(p.Ndof_s(1:s));
                    dofStart = dofEnd - p.Ndof_s(s) + 1;
                    fsPost{i}{s} = p.L{s}'*p.fs{s}(:,tiPlot(i));
                    NodalPos{i}{s} = NodalPos0{s} + p.L{s}'*p.u(dofStart:dofEnd,tiPlot(i));
                end
                ylim(1) = min(min(NodalPos{i}{s}(2:2:end)) - PlotMargin,ylim(1));
                ylim(2) = max(max(NodalPos{i}{s}(2:2:end)) + PlotMargin,ylim(2));
                xlim(1) = min(min(NodalPos{i}{s}(1:2:end)) - PlotMargin,xlim(1));
                xlim(2) = max(max(NodalPos{i}{s}(1:2:end)) + PlotMargin,xlim(2));
            end
            tracking(1,i)=NodalPos{i}{s}(end-1);
            tracking(2,i)=NodalPos{i}{s}(end);
        end
        disp('NodalPos')
        disp(NodalPos{end}{2})
        figure(854756);
        hold on
        %plot([tPlot,tracking(1,:)]);
        plot(tPlot',tracking(2,:));
        hold off
        
        % plot frames
        if p.PlotDeform
            figure(p.Figh.structure);
            set(p.Figh.structure,'units','normalized','outerposition',[0 0 1 1]);
            axis off;
            
            if p.WriteVideo
                Frames = struct('cdata', cell(1,length(PlotSteps)), 'colormap', cell(1,length(PlotSteps)));
            end
            
            SubsToPlot = 1:p.Ns;
            i = 0;
            for ti = tiPlot
                i = i + 1;
                timePlot = p.t(ti);
                description = ['time = ' num2str(timePlot) 's'];
                clf;
                FETIplot( p, p.Figh.structure, NodalPos{i}, SubsToPlot, fsPost{i}, description, xlim, ylim);
                set(gcf,'color','w');
                set(gca,'xtick',[],'ytick',[]);
                axis off;
                drawnow
                %plot2svg([p.Experiment '\PLOT_n' num2str(i) '.svg']);
                if p.WriteVideo
                    Frames(i) = getframe;
                end
            end
        end
        
        %plot2svg([p.Experiment '\PLOT_n' num2str(i) '.svg']);
        
        %% render video file
        if p.WriteVideo && ~isempty(Frames)
            writerObj = VideoWriter([p.Experiment '\Animation.mp4'],'MPEG-4');
            writerObj.FrameRate = 1/timestepPlot;
            writerObj.Quality = 100;
            open(writerObj);
            for frame = Frames
                writeVideo(writerObj,frame);
            end
            close(writerObj);
        end
        
        %% Te
        % plot iterations per time step
        if isfield(p,'PlotTe') && p.PlotTe
            figure(2) 
            subplot(2,1,1)
            h = title('Number of iterations per time step','FontSize',14);
            hold on
            grid on
            plot(p.it_ti);
            xlabel('timesteps','FontSize',14);
            ylabel('iterations','FontSize',14);
            legend(['total iterations:' num2str(sum(p.it_ti))]);
            axis([0 length(p.it_ti) 0 max(p.it_ti)+5])
            hold off
            % plot calculation time per time step
            subplot(2,1,2)
            title('Calculation time per time step','FontSize',14);
            hold on
            grid on
            plot(p.ti_ti);
            xlabel('timesteps','FontSize',14);
            ylabel('calculation time [s]','FontSize',14);
            legend(['total time:' num2str(sum(p.ti_ti)) 's']);
            %axis([0 length(p.ti_ti) 0 max(p.it_ti)+5])
            axis([0 length(p.ti_ti) 0 0.06])
            hold off
        end



        %%
        
        %% plot geneo modes
        if p.PlotGeneo
            s = p.PlotGeneoSubs;
            SelectedEigenvalues{s} = [1:size(p.GeneoEigenValues{s},1)];
            
            if (size(p.GeneoEigenValues{s},1) > p.PlotGeneoMax)
                SelectedEigenvalues{s} = 1:p.PlotGeneoMax;
            end
            SelectedEigenvalues{s} = [SelectedEigenvalues{s}; ...
                            p.GeneoEigenValues{s}(SelectedEigenvalues{s})'];
            
            % SelectedEigenvalues{s}: [Eigenvalue Numbers (asc. order);
            %                          Corresponding Eigenvalues         ]
            
            if isfield(p,'Odd') && p.NoPattern && ~isempty(p.Odd)
                figure('name',['Substrukur: ' num2str(p.Odd)])
            else
                figure(p.Figh.geneomodes);
                set(p.Figh.geneomodes,'units','normalized','outerposition',[0 0 1 1]);                
            end
            hold on;
            PlotCols = 4;
            PlotRows = ceil(size(SelectedEigenvalues{s},2)/PlotCols);
            i = 0;
            for EigenValueNr = SelectedEigenvalues{s}(1,:)
                if p.nonconforming==1
                    PlotMargin=p.elHeight(s);%p.sizes(2,s)/p.Nely(s);
                else
                    PlotMargin=p.elHeight;
                end
                i = i + 1;
                buMode = p.VqGen{s}(:,EigenValueNr);
                uMode{i}(p.bDOF{s},1) = buMode;
                uMode{i}(p.iDOF{s},1) = - inv(p.Ksii{s})*p.Ksib{s}*buMode;
                uMode{i} = p.L{s}'*uMode{i};
                uMode{i} = NodalPos0{s} + uMode{i}./(norm(uMode{i}));
                NodalPos{s} = 0.3*uMode{i};

                ylim(1) = min(NodalPos{s}(2:2:end)) - PlotMargin;
                ylim(2) = max(NodalPos{s}(2:2:end)) + PlotMargin;
                xlim(1) = min(NodalPos{s}(1:2:end)) - PlotMargin;
                xlim(2) = max(NodalPos{s}(1:2:end)) + PlotMargin;

                subplot(PlotRows,PlotCols,i);
                description = ['\lambda_{' num2str(SelectedEigenvalues{s}(1,i)) '} = ' num2str(SelectedEigenvalues{s}(2,i))];
                %description = [];
                SubsToPlot = s;
                %% Te - tausche x und y-Achse
                %g = 1:2:length(NodalPos{SubsToPlot});
                %h = 2:2:length(NodalPos{SubsToPlot});
                %he = NodalPos{SubsToPlot};
                %NodalPos{SubsToPlot}(g) = he(h);    
                %NodalPos{SubsToPlot}(h) = he(g);
                %he = xlim;
                %xlim = ylim;
                %ylim = he;
                %%
                %FETIplot( p, NodalPos, SubsToPlot, fsPost{1}, description, xlim, ylim);
                FETIplot( p, p.Figh, NodalPos, SubsToPlot, fsPost{1}, description, xlim, ylim) %FigureHandle
                set(gcf,'color','w');
                set(gca,'xtick',[],'ytick',[]);
            end

            hold off;
            if ~p.NoWrite
                %plot2svg(['case_' num2str(p.Case) '_geneo.svg']);
            end
        end
        
        
    %% Static
    elseif strcmp(p.mode,'static')
        for i = length(p.PlotIterations):-1:1
            if p.PlotIterations(i) > p.StaticIterations
                p.PlotIterations(i) = [];
            end
        end

        if strcmp(p.Solver,'FULLsolver')
            p.PlotIterations = 1;
        else
            if p.PlotLastIteration && isempty(p.PlotIterations)
                p.PlotIterations(1) = p.StaticIterations;
            elseif p.PlotLastIteration && p.PlotIterations(end) ~= p.StaticIterations
                p.PlotIterations(end+1) = p.StaticIterations;
            end
        end

        
        n = 0;
        for iteration = p.PlotIterations
            % for every iterationto plot the following plotting procedure
            % is carried out:
            
            n = n + 1;
            descriptions{1} = iteration;
            descriptions{2} = { ...
                'Actual State', ...
                'Actual State Without RBM', ...
                'Actual State Exploded', ...
                'Undeformed State', ...
                'Undeformed State Exploded', ...
                'Actual State Without RBM Exploded'};

            Nrbm_buf = 0;
            alphas = cell(1,p.Ns);
            NumPlots = 6;
            flms = cell(NumPlots,p.Ns);
            
            if strcmp(p.Solver,'FULLsolver')
                u_direct_loc=p.L_man*p.uFull;
                fsPost_glob = p.L_man*p.Full.f;
                end_u=0;
                %disp('u_direct_loc:')
                %disp(size(u_direct_loc))
            end
            for s=1:p.Ns
                if p.nonconforming ==1
                    gapx = 0.06*s;
                    gapy = 0.06*s;
                else
                    iNsy=int8(s/p.Nsy)+1;
                    iNsx=mod(s,p.Nsx);
            %for iNsy = 1:p.Nsy
             %   for iNsx = 1:p.Nsx
                    %s = (iNsy-1)*p.Nsx + iNsx;
                    gapx = 0.06*iNsx;
                    gapy = 0.06*iNsy;
                end
                uAdd = [];
                uAdd(1:2:size(p.L{s}',1),1) = gapx;
                uAdd(2:2:size(p.L{s}',1),1) = gapy;
                
                if strcmp(p.Solver,'FULLsolver')
                    start_u=1+end_u;
                    end_u=end_u+size(p.Bs{s},2);
                    %disp(['end_u: ' num2str(end_u)])
                    uDeformation = p.L{s}'*u_direct_loc(start_u:end_u);

                    fsPost{s} = fsPost_glob(start_u:end_u,1);
                    uRBM = zeros(size(uDeformation));
                    for i=1:NumPlots
                        flms{i}{s}=0;
                    end
                else
                   
                    alphas{s} = p.alpha(Nrbm_buf+1:Nrbm_buf+p.Nrbm_s(s),iteration);
                    uDeformation = p.L{s}'*(p.Ksp{s}*(p.fs{s}-p.Bs{s}'*p.lambda(:,iteration)));
                    uRBM = p.L{s}'*(p.RsK{s}*alphas{s});
                    
                    % visualize lambdas:

                    % actual state exploded
                    flms{1}{s} = -p.L{s}'*p.Bs{s}'*p.lambda(:,iteration);
                    flms{2}{s} = flms{1}{s};
                    flms{3}{s} = flms{1}{s};
                    flms{4}{s} = flms{1}{s};

                    % undeformed state exploded
                    demoLambda(1:p.Nlm,1) = 0;
                    demoLambda(1:2:p.Nlm,1) = 1;
                    flms{5}{s} = -p.L{s}'*p.Bs{s}'*demoLambda;

                    % actual state without rbm exploded
                    flms{6}{s} = flms{1}{s};   
                    
                    Nrbm_buf = Nrbm_buf + p.Nrbm_s(s);
                end
                us{1}{s} = uDeformation - uRBM;
                us{2}{s} = uDeformation;
                us{3}{s} = uDeformation - uRBM + uAdd;
                us{4}{s}(1:size(p.L{s}',1),1) = 0;
                us{5}{s}(1:size(p.L{s}',1),1) = uAdd;
                us{6}{s} = uDeformation + uAdd;

            end
            
            
            
            %% identify reasonable values for xlim and ylim
            NodalPosMult = cell(NumPlots,p.Ns);
            for s = 1:p.Ns
                dofEnd = sum(p.Ndof_s(1:s));
                dofStart = dofEnd - p.Ndof_s(s) + 1;
                if strcmp(p.Solver,'FULLsolver')==0
                    fsPost{s} = p.L{s}'*p.fs{s}(:,1);
                end

                %if p.nonconforming==1
                %    PlotMargin=p.sizes(2,s)/p.Nely(s);
                %else
                if p.nonconforming==1
                    PlotMargin=p.elHeight(s);
                else
                    PlotMargin=p.elHeight;
                end
                %end
                for i = p.StaticPlots
                    fAll{s} = fsPost{s} + flms{i}{s};
                    NodalPosMult{i}{s} = NodalPos0{s} + us{i}{s};
                    ylim(1) = min(min(NodalPosMult{i}{s}(2:2:end)) - PlotMargin,ylim(1));
                    ylim(2) = max(max(NodalPosMult{i}{s}(2:2:end)) + PlotMargin,ylim(2));
                    xlim(1) = min(min(NodalPosMult{i}{s}(1:2:end)) - PlotMargin,xlim(1));
                    xlim(2) = max(max(NodalPosMult{i}{s}(1:2:end)) + PlotMargin,xlim(2));
                end
            end
            p.tracking=us{1}{p.tracking_dof(1)}(p.tracking_dof(2));
            disp('Deformation of rightmost and uppermost point:')
            disp(p.tracking)
            
            
            
            %% plot different displacements
            if p.PlotDeform 
                if ~exist('DisplacementPlots','var')
                    DisplacementPlots = [];
                end
                if p.PlotDeformClosePrev
                    close(DisplacementPlots);
                end
                
                for iPlot = p.StaticPlots % array
                    description = [descriptions{2}{iPlot} '_it' num2str(iteration)];
                    PlotDescription = [descriptions{2}{iPlot} ' at iteration ' num2str(iteration)];

                    FigH = p.Figh.structure*10000+1000*iPlot+iteration;
                    figure(FigH);
                    set(FigH,'units','normalized','outerposition',[0 0 1 1]);
                    hold on;
                    axis off;
                    SubsToPlot = 1:p.Ns;
                    %FETIplot( p, FigH, NodalPosMult{iPlot}, SubsToPlot, fsPost, PlotDescription, xlim, ylim);
                    FETIplot( p, FigH, NodalPosMult{iPlot}, SubsToPlot, fAll, PlotDescription, xlim, ylim);
                    set(gcf,'color','w');
                    set(gca,'xtick',[],'ytick',[]);
                    hold off;
                    if ~p.NoWrite
                        plot2svg([p.Experiment '\case' num2str(p.Case) '-' regexprep(description,'[^\w'']','') '.svg']);
                    end

                    DisplacementPlots(end+1) = FigH;
                end
            end
        end
        
        
        
        %% plot geneo modes
        if p.PlotGeneo
            s = p.PlotGeneoSubs;
            SelectedEigenvalues{s} = [1:size(p.GeneoEigenValues{s},1)];
            
            if (size(p.GeneoEigenValues{s},1) > p.PlotGeneoMax)
                SelectedEigenvalues{s} = 1:p.PlotGeneoMax;
            end
            SelectedEigenvalues{s} = [SelectedEigenvalues{s}; ...
                            p.GeneoEigenValues{s}(SelectedEigenvalues{s})'];
            
            % SelectedEigenvalues{s}: [Eigenvalue Numbers (asc. order);
            %                          Corresponding Eigenvalues         ]

            figure(p.Figh.geneomodes);
            set(p.Figh.geneomodes,'units','normalized','outerposition',[0 0 1 1]);
            hold on;
            PlotCols = 4;
            PlotRows = ceil(size(SelectedEigenvalues{s},2)/PlotCols);
            i = 0;
            for EigenValueNr = SelectedEigenvalues{s}(1,:)
                if p.nonconforming==1
                    PlotMargin=p.elHeight(s);%p.sizes(2,s)/p.Nely(s);
                else
                    PlotMargin=p.elHeight;
                end
                i = i + 1;
                buMode = p.VqGen{s}(:,EigenValueNr);
                uMode{i}(p.bDOF{s},1) = buMode;
                uMode{i}(p.iDOF{s},1) = - inv(p.Ksii{s})*p.Ksib{s}*buMode;
                uMode{i} = p.L{s}'*uMode{i};
                uMode{i} = NodalPos0{s} + uMode{i}./(norm(uMode{i}));
                NodalPos{s} = uMode{i};

                ylim(1) = min(NodalPos{s}(2:2:end)) - PlotMargin;
                ylim(2) = max(NodalPos{s}(2:2:end)) + PlotMargin;
                xlim(1) = min(NodalPos{s}(1:2:end)) - PlotMargin;
                xlim(2) = max(NodalPos{s}(1:2:end)) + PlotMargin;

                subplot(PlotRows,PlotCols,i);
                description = ['\lambda_{' num2str(SelectedEigenvalues{s}(1,i)) '} = ' num2str(SelectedEigenvalues{s}(2,i))];
                %description = [];
                SubsToPlot = s;
                FETIplot( p, p.Figh.geneomodes, NodalPos, SubsToPlot, fsPost, description, xlim, ylim);

                set(gcf,'color','w');
                set(gca,'xtick',[],'ytick',[]);
            end

            hold off;
            if ~p.NoWrite
                %plot2svg(['case_' num2str(p.Case) '_geneo.svg']);
            end
        end
    end
end