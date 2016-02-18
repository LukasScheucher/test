function FETIplotMACcoarse(p)
MAConlyForSub = 3;
%% mac plot between geneo and rbms
CorrCols = p.Nsx;
CorrRows = ceil(p.nFloating/CorrCols);
close all;
figure(p.Figh.mac);
if MAConlyForSub > 0
    set(p.Figh.mac,'units','normalized','outerposition',[0.1 0.1 0.3 0.5]);
else
    set(p.Figh.mac,'units','normalized','outerposition',[0 0 1 1]);
end
set(gcf,'color','w');

ploti = 0;
for iNsy = p.Nsy:-1:1
    for iNsx = 1:p.Nsx
        s = (iNsy-1)*p.Nsx + iNsx;
        ploti = ploti+1;
        if size(p.Rs{s},2) == 0
            continue;
        end

        if MAConlyForSub > 0
            if s ~= MAConlyForSub
                continue;
            end
        else
            subplot(CorrRows,CorrCols,ploti);
        end

        GeneoNorms = sqrt(sum((p.VqGen{s}(:,1:3)).^2,1));
        RNorms = sqrt(sum((p.Rs{s}(:,1:3)).^2,1));
        Norming = inv(diag([GeneoNorms, RNorms]));

        Corr = [p.VqGen{s}(:,1:3), p.Rs{s}(:,1:3)]*Norming;
        CorrVals{s} = sort(eig(Corr'*Corr));
        semilogy(CorrVals{s},'*');
        %title(['Subs. ', num2str(s)]);
        title([num2str(s)]);
        %set(gca,'YGrid','on')
        %set(gca,'YMinorGrid','off');
        %set(gca,'XTick',[1 2 3 4 5 6]);
        xlim([0 7]);
        %set(gca,'YTick',[1e-10 1e-8 1e-6 1e-4 1e-2 0.1 1]);
        %ylim([min(CorrVals{s})/22 max(CorrVals{s})*5]);
        ylim([1e-14 1e2]);
        %set(gca,'YTickLabel',{'';1e-3;'';'$10^{-1}$';'$1$';});
    end
end
if ~p.NoWrite
    plot2svg([p.Experiment '\Ns' num2str(p.Nsx) '_MAC_h' num2str(Integration.Step) '_' CoarseGrid '.svg']);
end
end