function Frames = FULLpost(Full, p, uFull)
    NodalPos = [];
    fsPost = [];
    xlim = [inf -inf];
    ylim = [inf -inf];
    PlotMargin = p.elHeight;
    
    % compute initial nodal positions
    for r = 1:(p.Nely*p.Nsy+1) % r: row (substructure-local)
        for c = 1:(p.Nelx*p.Nsx+1) % c: column (substructure-local)
            i = (r-1)*(p.Nelx*p.Nsx+1)+c; % i: node nr
            NodalPos0{1}(2*i-1:2*i,1) = [p.elHeight*c; p.elHeight*r];
        end
    end
    
    tPlot = 0:p.Plotstep:p.End;
    
    n = 0;
    for time = p.t
        n = n + 1;
        % blow up external forces to full coordinates
        fsPost{n}{1} = p.Lfull'*Full.f(:,n);
        NodalPos{n}{1} = NodalPos0{1} + p.Lfull'*uFull(:,n);

        ylim(1) = min(min(NodalPos{n}{1}(2:2:end)) - PlotMargin,ylim(1));
        ylim(2) = max(max(NodalPos{n}{1}(2:2:end)) + PlotMargin,ylim(2));
        xlim(1) = min(min(NodalPos{n}{1}(1:2:end)) - PlotMargin,xlim(1));
        xlim(2) = max(max(NodalPos{n}{1}(1:2:end)) + PlotMargin,xlim(2));
    end

    figure(p.Figh.structure);
    %set(p.Figh.structure,'units','pixels','outerposition',[0 0 1920 1080]);
    set(p.Figh.structure,'units','normalized','outerposition',[0 0 1 1]);
    
    % plot frames
    n = 0;
    for timePlotSet = tPlot
        n = n + 1;
        i = max(1,floor(timePlotSet/p.Step));
        timePlot = p.t(i);
        
        Fullgeometry = feti;
        Fullgeometry.Nsx = 1;
        Fullgeometry.Nsy = 1;
        Fullgeometry.Nelx = p.Nelx*p.Nsx;
        Fullgeometry.Nely = p.Nely*p.Nsy;
        
        FETIplot( Fullgeometry, feti, NodalPos{i}, fsPost{i}, timePlot, xlim, ylim);
        set(gcf,'color','w');
        set(gca,'xtick',[],'ytick',[]);
        
        Frames(n) = getframe;
    end
end