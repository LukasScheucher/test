function FETIplotOperatorEig(p)
%% plot operator eigenvalues

if strcmp(p.mode,'dynamic')
    Op = p.FItrans;
    OpPre = p.FItransPre;
    P = p.Ptrans;
elseif strcmp(p.mode,'static')
    Op = p.FI;
    OpPre = p.FIPre;
    P = p.P;
end

EigenvaluesOp = [];
EigenvaluesPrecondOp = [];
EigenvaluesProjPrecondOp = [];

EigenvaluesOp = sort(eig(Op));
if p.Coarse
    EigenvaluesPrecondOp = sort(eig(OpPre*Op));
    EigenvaluesProjPrecondOp = sort(eig(OpPre*P'*Op*P));
end

abs(EigenvaluesPrecondOp)
figh = p.Figh.eig + p.Case;
figure(figh);
set(figh,'units','normalized','outerposition',[0 0 1 1]);

numCrossPoints = (p.Nsx-1)*(p.Nsy-1)*6;

subplot(3,1,1);
semilogy(abs(EigenvaluesOp(numCrossPoints+1:end)),'*');
xlim([0 (p.Nlm-numCrossPoints)]);
title(['Before Preconditioning (Case ' num2str(p.Case) ')']);
if p.Coarse
    subplot(3,1,2);
    semilogy(abs(EigenvaluesPrecondOp(numCrossPoints+1:end)),'*');
    xlim([0 (p.Nlm-numCrossPoints)]);
    title(['After Preconditioning (Case ' num2str(p.Case) ')']);
    subplot(3,1,3);
    semilogy(abs(EigenvaluesProjPrecondOp(numCrossPoints+1:end)),'*');
    xlim([0 (p.Nlm-numCrossPoints)]);
    title(['With Coarse Grid (Case ' num2str(p.Case) ')']);
end


% figure(8);
% set(8,'units','normalized','outerposition',[0 0 1 1]);
% 
% binranges = [1e-8 0.1 0.8 1.2 10 100 1000 10000];
% 
% hold on;
% 
% subplot(3,1,1);
% [bincounts] = histc(Eigenvalues1(numCrossPoints+1:end),binranges);
% bar(binranges,bincounts,'histc');
% set(gca, 'Xscale', 'log')
% 
% if p.Coarse
%     subplot(3,1,2);
%     [bincounts] = histc(Eigenvalues3(numCrossPoints+1:end),binranges);
%     bar(binranges,bincounts,'histc');
%     set(gca, 'Xscale', 'log')
% 
%     subplot(3,1,3);
%     [bincounts] = histc(Eigenvalues2(numCrossPoints+1:end),binranges);
%     bar(binranges,bincounts,'histc');
%     set(gca, 'Xscale', 'log')
% end
% 
% hold off;


set(gcf,'color','w');
if ~p.NoWrite
%plot2svg([p.Experiment '\Eig_h' num2str(p.Step) '_' p.Coarse '.svg']);
end
end