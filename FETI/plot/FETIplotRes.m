function FETIplotRes(p)

    figure(p.Figh.its);
    set(p.Figh.its,'units','normalized','outerposition',[0 0 0.3 0.6]);
    %resPlot = sqrt(sum((w).^2,1));
    %semilogy(resPlot);
    semilogy(p.Residual);
    set(gca,'YGrid','on')
    set(gca,'YMinorGrid','off');
    %set(gca,'YTick',[1e-10 1e-8 1e-6 1e-4 1e-2 1 1e2 1e4]);
    %set(gca,'XTick',[0 20 40 60 80 100 120 140 160]);
    xlim([0 p.StaticIterations+1]);
    disp('ylim:')
    disp(p.Residual)
    ylim([min(p.Residual)/7 max(p.Residual)*7]);
    set(gcf,'color','w');
    if ~p.NoWrite
        %plot2svg([Experiment '\Ns' num2str(p.Nsx) '_It_n' num2str(n??) '_h' num2str(Step??) '_' CoarseGrid '.svg']);
    end
    
end