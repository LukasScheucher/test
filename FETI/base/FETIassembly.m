function [p] = FETIassembly(p)
tic
% compute element stiffness matrix

% set last argument to 1 for plane stress (elements are thin and have free
% surfaces at their upside/downside, can change thickness
if p.nonconforming==1
    s_end=p.Ns; % Assembly of matrices has to be done separately for every substructure in case of nonconforming meshes
else
    s_end=1;
end
p.Ls_dof=0;
tol=p.geom_tol;
for s=1:s_end
    %p.Ls_ass{s} = null(p.Bs{s});
    %disp(['Dimension Ls{' num2str(s) '}:'])
    %disp(size(p.Ls_ass{s}))

    %disp(p.Bs{s}*p.Ls_ass{s})
    %if p.nonconforming==1
    %    p.Ls_dof=p.Ls_dof+size(p.Ls_ass{s},2);
    %else
    %    p.Ls_dof=size(p.Ls_ass{s},2)*p.Ns;
    %end
    %if s==1
    %    Full.K_multi=zeros(size(p.Ls_ass{s},2));
    %end
    kSt{s} = BilinearQuadElementStiffness(p.ESt,p.nuSt,p.elThick, 0,0, p.elHeight(s),0, p.elHeight(s),p.elHeight(s), 0,p.elHeight(s), p.plain);
    kRub{s} = BilinearQuadElementStiffness(p.ERub,p.nuRub,p.elThick, 0,0, p.elHeight(s),0, p.elHeight(s),p.elHeight(s), 0,p.elHeight(s), p.plain);

    % [E] = Einheit E-Modul
    % [s] = Einheit Elementdicke, Elementkoordinaten
    % [Kelement] = [s]*[E]
    %
    % z.B. [E] = kN/m^2 (z.B. ESt = 210e6 kN/m^2)
    %      [s] = m
    % -->  [Kelement] = kN/m
    % mit  [u] = m
    % -->  [f] = kN !!!

    % compute element mass matrix for bilinear quadrilateral element
    A = p.elHeight(s)^2;
    mPre = [4 2 1 2; ...
            2 4 2 1; ...
            1 2 4 2; ...
            2 1 2 4];
    mPre = [mPre zeros(4); ... % x direction
            zeros(4) mPre];    % y direction
    LmassEl = [1 0 0 0 0 0 0 0; ...
               0 0 1 0 0 0 0 0; ...
               0 0 0 0 1 0 0 0; ...
               0 0 0 0 0 0 1 0; ...
               0 1 0 0 0 0 0 0; ...
               0 0 0 1 0 0 0 0; ...
               0 0 0 0 0 1 0 0; ...
               0 0 0 0 0 0 0 1];
    MelementWithoutRHO = A*p.elThick/36.*LmassEl'*mPre*LmassEl;

    % [rho] = Einheit Dichte
    % [s] = Einheit Dicke, Elementhöhe etc.
    % [Melement] = [rho]*[s]^3
    %
    % z.B. [rho] = g/cm^3 = kN/((m/s^2)*m^3) (z.B. rhoSt = 7,85 g/cm^3)
    %      [s] = m
    % -->  [Melement] = kN/(m/s^2)
    % mit  [u] = m
    % -->  [f] = kN !!!

    % assemble global stiffness matrix K and global mass matrix M of one substructure
    KprotoOdd = zeros(p.Ndof_s_Free(s));
    MprotoOdd = zeros(p.Ndof_s_Free(s));
    KprotoEven = zeros(p.Ndof_s_Free(s));
    MprotoEven = zeros(p.Ndof_s_Free(s));


    for r = 1:p.Nely(s) % r: row
        for c = 1:p.Nelx(s) % c: column
            i = (r-1)*p.Nelx(s)+c; % i: element number
            Lass = zeros(8,p.Ndof_s_Free(s));
            for node2ass = 1:4
                Lass(node2ass*2,2*p.ec{s}(i,node2ass)-1) = 1;
                Lass(node2ass*2-1,2*p.ec{s}(i,node2ass)) = 1;
            end

            if (~isempty(find(r == p.SteelRowNrsOddNsx,1)))
                KprotoOdd = BilinearQuadAssemble(KprotoOdd,kSt{s},p.ec{s}(i,1),p.ec{s}(i,2),p.ec{s}(i,3),p.ec{s}(i,4));
                MprotoOdd = MprotoOdd + Lass'*(p.rhoSt.*MelementWithoutRHO)*Lass;
                %display(['row ' num2str(r) ': steel']);
            else
                KprotoOdd = BilinearQuadAssemble(KprotoOdd,kRub{s},p.ec{s}(i,1),p.ec{s}(i,2),p.ec{s}(i,3),p.ec{s}(i,4));
                MprotoOdd = MprotoOdd + Lass'*(p.rhoRub.*MelementWithoutRHO)*Lass;
                %display(['row ' num2str(r) ': rubber']);
            end

            if (~isempty(find(r == p.SteelRowNrsEvenNsx,1)))
                KprotoEven = BilinearQuadAssemble(KprotoEven,kSt{s},p.ec{s}(i,1),p.ec{s}(i,2),p.ec{s}(i,3),p.ec{s}(i,4));
                MprotoEven = MprotoEven + Lass'*(p.rhoSt.*MelementWithoutRHO)*Lass;
                %display(['row ' num2str(r) ': steel']);
            else
                KprotoEven = BilinearQuadAssemble(KprotoEven,kRub{s},p.ec{s}(i,1),p.ec{s}(i,2),p.ec{s}(i,3),p.ec{s}(i,4));
                MprotoEven = MprotoEven + Lass'*(p.rhoRub.*MelementWithoutRHO)*Lass;
                %display(['row ' num2str(r) ': rubber']);
            end
        end
    end



    % distribute the prototype of the substructure stiffness and mass matrix
    if (p.NoPattern == 0)
        if p.nonconforming==1
            %for s=1:p.Ns
            if mod(s,2) || ~p.ChangeForEvenNsy
                NsxTest = s;
            else
                NsxTest = s+1;
            end

            if mod(NsxTest,2)
               p.Ks{s} = KprotoOdd;
               p.Ms{s} = MprotoOdd;
            else
               p.Ks{s} = KprotoEven;
               p.Ms{s} = MprotoEven;
            end    
            %end
        else
            for iNsy = 1:p.Nsy
                for iNsx = 1:p.Nsx
                    i = (iNsy-1)*p.Nsx + iNsx;

                    if mod(iNsy,2) || ~p.ChangeForEvenNsy
                       NsxTest = iNsx;
                    else
                       NsxTest = iNsx+1;
                    end

                    if mod(NsxTest,2)
                       p.Ks{i} = KprotoOdd;
                       p.Ms{i} = MprotoOdd;
                    else
                       p.Ks{i} = KprotoEven;
                       p.Ms{i} = MprotoEven;
                    end
                end
            end
        end

    else
        %% TE   
        % Unterscheide nicht zwischen odd und even sondern wähle speziell aus
        % welche Substruktur wie sein soll (immernoch nur Reihenweise Unterschiede
        % pro Substruktur)
        if p.nonconforming==1
            if (~isempty(find(s == p.Odd,1)))
                p.Ks{s} = KprotoOdd; % Odd bedeutet hier: die "anderen" Substrukturen
                p.Ms{s} = MprotoOdd;
            else
               p.Ks{s} = KprotoEven;
               p.Ms{s} = MprotoEven;
            end          
        else
            for i = 1:p.Ns
                if (~isempty(find(i == p.Odd,1)))
                    p.Ks{i} = KprotoOdd; % Odd bedeutet hier: die "anderen" Substrukturen
                    p.Ms{i} = MprotoOdd;
                else
                   p.Ks{i} = KprotoEven;
                   p.Ms{i} = MprotoEven;
                end      
            end
        end
    end
    
    %Full.K_multi=Full.K_multi+p.Ls_ass{s}'*p.Ks{s}*p.Ls_ass{s};
end
% #Zeitmessung l2-178
time1=toc;
time_rel1=time1/(178-2);
%%
%% assemble matrices for direct solver
tic
if p.nonconforming ~=1
    Full.M = zeros(p.NdofFull);
    Full.K = Full.M;
    Full.f = zeros(p.NdofFull,length(p.t));

    Full.K = sparse(Full.K);
    Full.M = sparse(Full.M);


    for iNsy = 1:p.Nsy
        for iNsx = 1:p.Nsx
            s = (iNsy-1)*p.Nsx + iNsx;
            for r = 0:p.Nely(s)
                base = (p.Nely(s)*(iNsy-1)*(p.Nelx(s)*p.Nsx+1)+(p.Nelx(s)*p.Nsx+1)*r+(iNsx-1)*p.Nelx(s)); 
                firstNode = r*(p.Nelx(s)+1)+1;
                lastNode = r*(p.Nelx(s)+1)+(p.Nelx(s)+1);
                ecFull{s}(firstNode*2-1:lastNode*2) = ((base+1)*2-1):((base+p.Nelx(s)+1)*2);
            end

            Full.M(ecFull{s},ecFull{s}) = Full.M(ecFull{s},ecFull{s}) + p.Ms{s};
            Full.K(ecFull{s},ecFull{s}) = Full.K(ecFull{s},ecFull{s}) + p.Ks{s};
            Full.f(ecFull{s},:) = Full.f(ecFull{s},:) + p.fs{s};
        end
    end

    % FULL: delete DBC rows and cols
    Full.K(p.zerodofsFull,:) = [];
    Full.M(p.zerodofsFull,:) = [];
    Full.K(:,p.zerodofsFull) = [];
    Full.M(:,p.zerodofsFull) = [];
    Full.f(p.zerodofsFull,:) = [];
end
% #Zeitmessung l182-216
time2=toc;
time_rel2=time2/(216-182);

% feti: delete DBC rows and cols
tic
for s = 1:p.Ns
    p.Ks{s}(p.zerodofs{s},:) = [];
    p.Ks{s}(:,p.zerodofs{s}) = [];
    p.Ms{s}(p.zerodofs{s},:) = [];
    p.Ms{s}(:,p.zerodofs{s}) = [];
    p.fs{s}(p.zerodofs{s},:) = [];

    p.Bs{s}(p.zerolms,:) = [];
    p.Bs{s}(:,p.zerodofs{s}) = [];
end
% #Zeitmessung l220-245
time3=toc;
time_rel3=time3/(245-220);
tic
p.B = p.Bs{1};
for s = 2:p.Ns
    p.B = [p.B p.Bs{s}];
end
p.B = horMatrix(p.Bs,p.Nlm,p.Ndof);
p.B2=p.B;

%{
for l=1:size(p.B2,1)
    if sum(abs(p.B2(l,:)))>tol
        crossnode=0;
        %{
        for n=1:size(p.B2,2)
            if abs(abs(p.B2(l,n))-1)<=tol
                crossnode=1;
            end
        end
        %}
        if crossnode==0
            LMnorm=sqrt(sum(p.B2(l,:).^2));
            p.B2(l,:)=p.B2(l,:)/LMnorm;
            for s=1:p.Ns
                p.Bs{s}=p.Bs{s}/LMnorm;
            end
        end
    end
end
%}
disp('B2:')
disp(p.B2)

disp(size(p.B2))
disp(rank(p.B2))

%p.L_ass=null(p.B2);

%disp('L:')
%disp(p.L_ass)
%disp(size(p.L_ass))

if p.globalassembly==1
    disp('L manuell:')

    [p] = LMatrix(p);

    disp(p.L_man)
    disp(size(p.L_man))
    %disp('B*L:')
    %disp(sum(abs(p.B2*p.L_man)))

    %disp('Kern(L):')
    %disp(null(p.L_man')')
    %disp(['Rang(L): ' num2str(rank(p.L_man))])

    %p.B=null(p.L_man')';zeros(size(p.B2,1)-rank(p.B2),size(p.B2,2))];

    %disp('Kern L:')
    %disp(size(p.B))
    %disp(rank(p.B))

    % #Zeitmessung l247-313
    time4=toc;
    time_rel4=time4/(313-247);
    tic
    p.K = diagMatrix(p.Ks,p.Ndof,p.Ndof);
    p.M = diagMatrix(p.Ms,p.Ndof,p.Ndof);
    Full.K2=p.L_man'*p.K*p.L_man;
    Full.M2=p.L_man'*p.M*p.L_man;

    if p.nonconforming~=1
        Full.K=full(Full.K);
        Full.M=full(Full.M);
    end
    Full.K=Full.K2;
    Full.K=sparse(Full.K);
    Full.M=Full.M2;
    Full.M=sparse(Full.M);
end
if p.nonconforming~=1
    p.CornerLMs(p.zerolms,:) = [];
    p.NdofFull = size(Full.K,1);
end
% #Zeitmessung l315-376
time5=toc;
time_rel5=time5/(376-315);
tic
%% eigenanalysis
if p.globalassembly==1
    if exist([p.Experiment '\Eigenanalysis.mat'], 'file')
        load([p.Experiment '\Eigenanalysis.mat']);

        display('ATTENTION: Loading saved Eigenanalysis!');
    else
        if p.DisplayEigen
            display('Recalc. Eigv.');
        end
        [Vfull, MuFullDiag] = eigs(Full.K,Full.M,10,'sm');
        %[Vfull, MuFullDiag] = eig(Full.K,Full.M);
        % (K - mu*M)*v = 0
        % Vfull = (v1, v2, ... )
        % MuFullDiag = diag(mu1, mu2, ...)
        if ~p.NoWrite
            save([p.Experiment '\Eigenanalysis.mat'],'Vfull','MuFullDiag');
        end
    end

    [Vfull, MuFullDiag] = sortem(Vfull,MuFullDiag);

    % (K - mu*M)*v = 0
    % (K - omega^2*M)*v = 0
    % --> omega = sqrt(mu)
    OmegaFull = sqrt(diag(MuFullDiag));

    % omega = 2*pi*f
    FreqFull = OmegaFull./(2*pi);
    PeriodsFull = 1./FreqFull;

    if p.DisplayEigen
        Frequencies = [1 2 3 4 5 10];
        for Freq = Frequencies
            if Freq > length(FreqFull)
                break;
            end
            display(['Frequency of the ' num2str(Freq) '. Eigenfrequency: ' num2str(FreqFull(Freq)) 'Hz  (Period:' num2str(PeriodsFull(Freq)) ' s)']);
        end
        %MidsizeDeltaT = SelectedPeriod50/10
    end
    % #Zeitmessung l378-421
    time6=toc;
    time_rel6=time6/(421-378);
    tic
    %% damping
    DampingRatioSteel = 0.001;
    DampingRatioRubber = 0.025;
    StiffnessCoefficient = 2*DampingRatioRubber/(OmegaFull(1)+OmegaFull(2));
    MassCoefficient = OmegaFull(1)*OmegaFull(2)*StiffnessCoefficient;
    if p.DisplayDamping
        display('Todo display damping');
    end

    if p.nonconforming==1
        for s=1:p.Ns
            p.Cs{s} = MassCoefficient.*p.Ms{s} + StiffnessCoefficient.*p.Ks{s};
        end
        p.C=diagMatrix(p.Cs,p.Ndof,p.Ndof);

        Full.C2=p.L_man'*p.C*p.L_man;
        %disp('Full.C2:')
        %disp(Full.C2)    

    else
        %#########
        for s=1:p.Ns
            p.Cs{s} = MassCoefficient.*p.Ms{s} + StiffnessCoefficient.*p.Ks{s};
        end
        p.C=diagMatrix(p.Cs,p.Ndof,p.Ndof);

        Full.C2=p.L_man'*p.C*p.L_man;
        %disp('Full.C2:')
        %disp(Full.C2)  
        %#########
        s=1;
        CprotoOdd = MassCoefficient.*MprotoOdd + StiffnessCoefficient.*KprotoOdd;
        CprotoEven = MassCoefficient.*MprotoEven + StiffnessCoefficient.*KprotoEven;
        Full.C = zeros((p.Nsx*p.Nelx(s)+1)*(p.Nsy*p.Nely(s)+1)*2);
        %}
        if (p.NoPattern == 0)        
            for iNsy = 1:p.Nsy
                for iNsx = 1:p.Nsx
                    s = (iNsy-1)*p.Nsx + iNsx;

                    if mod(iNsy,2) || ~p.ChangeForEvenNsy
                        NsxTest = iNsx;
                    else
                        NsxTest = iNsx+1;
                    end

                    if mod(NsxTest,2)
                        p.Cs{s} = CprotoOdd;
                        %disp(['size p.Cs{' num2str(s) '}:'])
                        %disp(size(p.Cs{s}))
                        if (length(p.zerodofs)>=s && ~isempty(p.zerodofs{s}))
                            p.Cs{s} = p.L{s}*p.Cs{s}*p.L{s}';
                        end
                        Full.C(ecFull{s},ecFull{s}) = Full.C(ecFull{s},ecFull{s}) + CprotoOdd;
                    else
                        p.Cs{s} = CprotoEven;
                        if (length(p.zerodofs)>=s && ~isempty(p.zerodofs{s}))
                            p.Cs{s} = p.L{s}*p.Cs{s}*p.L{s}';
                        end
                        Full.C(ecFull{s},ecFull{s}) = Full.C(ecFull{s},ecFull{s}) + CprotoEven;
                    end
                end
            end
            %}
        else 

            %% TE   
            % Unterscheide nicht zwischen odd und even sondern wähle speziell aus
            % welche Substruktur wie sein soll (immernoch nur Reihenweise Unterschiede
            % pro Substruktur)
            for s = 1:p.Nsx*p.Nsy
                if (~isempty(find(s == p.Odd,1)))
                    p.Cs{s} = CprotoOdd;
                    if (length(p.zerodofs)>=s && ~isempty(p.zerodofs{s}))
                        p.Cs{s} = p.L{s}*p.Cs{s}*p.L{s}';
                    end
                    Full.C(ecFull{s},ecFull{s}) = Full.C(ecFull{s},ecFull{s}) + CprotoOdd;
                else
                    p.Cs{s} = CprotoEven;
                    if (length(p.zerodofs)>=s && ~isempty(p.zerodofs{s}))
                        p.Cs{s} = p.L{s}*p.Cs{s}*p.L{s}';
                    end
                    Full.C(ecFull{s},ecFull{s}) = Full.C(ecFull{s},ecFull{s}) + CprotoEven;
                end       
            end
        end
        if p.nonconforming~=1
            Full.C(p.zerodofsFull,:) = [];
            Full.C(:,p.zerodofsFull) = [];
            % Full.C = p.Lfull*Full.C*p.Lfull';
            % Full.C = MassCoefficient.*Full.M + StiffnessCoefficient.*Full.K;
        end
        %disp('Full C:')
        Full.C_diag=zeros(1,size(Full.C2,2));
        for i=1:size(Full.C2,2)
            Full.C2_diag(i)=Full.C2(i,i);
        end
        %disp(Full.C)
    end

    Full.C=Full.C2;
end
% #Zeitmessung l423-526
time7=toc;
time_rel7=time7/(526-423);
tic
for s=1:p.Ns
    if strcmp(p.mode,'dynamic')

        % build up D matrix (generalized alpha)
        p.Ds{s} = p.Ms{s} + p.Step*p.GAgamma*(1-p.GAalphaf)/(1-p.GAalpham).*p.Cs{s} + p.Step^2*p.GAbeta*(1-p.GAalphaf)/(1-p.GAalpham).*p.Ks{s};

        % compute (b)oundary(=interface!) DOF parts and (i)nternal DOF parts
        % (b)oundary = interface, NOT boundaries without neighbours
        Dsii{s} = p.Ds{s}(p.iDOF{s},p.iDOF{s});
        Dsib{s} = p.Ds{s}(p.iDOF{s},p.bDOF{s});
        Dsbi{s} = p.Ds{s}(p.bDOF{s},p.iDOF{s});
        p.Dsbb{s} = p.Ds{s}(p.bDOF{s},p.bDOF{s});

        % Schur Complement with D:
        p.SDs{s} = (p.Dsbb{s} - Dsbi{s}*(Dsii{s}\Dsib{s}));
    end

    % compute (b)oundary DOF part of B matrices
    
    p.Bbs{s} = p.Bs{s}(:,p.bDOF{s});

    % compute (b)oundary DOF parts and (i)nternal DOF parts
    % (b)oundary = interface, NOT boundaries without neighbours
    p.Ksii{s} = p.Ks{s}(p.iDOF{s},p.iDOF{s});
    p.Ksib{s} = p.Ks{s}(p.iDOF{s},p.bDOF{s});
    p.Ksbi{s} = p.Ks{s}(p.bDOF{s},p.iDOF{s});
    p.Ksbb{s} = p.Ks{s}(p.bDOF{s},p.bDOF{s});

    % pseudoinverse (--> check if still needed, do not use)
    p.Ksp{s} = pinv(p.Ks{s});

    % Schur Complement with K:
    p.SKs{s} = (p.Ksbb{s} - p.Ksbi{s}*(p.Ksii{s}\p.Ksib{s}));
end

%% build global feti matrices
%{ 
Moved some lines up for calculation of L-matrix 
p.B = p.Bs{1};
for s = 2:p.Ns
    p.B = [p.B p.Bs{s}];
end
%}
p.f = [];
for i = 1:p.Ns
    p.f = [p.f;p.fs{i}];
    disp(['p.fs{' num2str(i) '}'])
    disp(p.fs{i})
end
assert(size(p.f,1)==p.Ndof,['Size of vector p.f is wrong! (size is ' num2str(size(p.f,1)) ' and should be ' num2str(p.Ndof) ')']);
if p.globalassembly==1
    Full.f2=p.L_man\p.f;
    Full.f=Full.f2;
    %disp('Full.f:')
    %disp(Full.f)
end
%p.K = diagMatrix(p.Ks,p.Ndof,p.Ndof);
if strcmp(p.mode,'dynamic')
    p.M = diagMatrix(p.Ms,p.Ndof,p.Ndof);
    p.C = diagMatrix(p.Cs,p.Ndof,p.Ndof);
    p.D = diagMatrix(p.Ds,p.Ndof,p.Ndof);

    p.Dbb = diagMatrix(p.Dsbb,p.Nbdof,p.Nbdof);
    p.SD = diagMatrix(p.SDs,p.Nbdof,p.Nbdof);
end
p.Kp = diagMatrix(p.Ksp,p.Ndof,p.Ndof);
p.Kbb = diagMatrix(p.Ksbb,p.Nbdof,p.Nbdof);
p.SK = diagMatrix(p.SKs,p.Nbdof,p.Nbdof);
%p.Bdiag = diagMatrix(p.Bs,p.Nlm*p.Ns,p.Ndof);

%p.B = horMatrix(p.Bs,p.Nlm,p.Ndof);
p.Bb = horMatrix(p.Bbs,p.Nlm,p.Nbdof);
%GAbetaScale = horMatrix(GAbetaScales,p.Nlm,p.Nlm*p.Ns);

% build up p.R matrix
p.R = [];
p.RK = [];
p.nFloating = 0;
p.Nrbm = 0;
for i = 1:p.Ns

%     p.Rs{i} = null(p.SKs{i});
%     p.RsK{i} = null(p.Ks{i});

    p.RsK{i} = null(p.Ks{i});
    p.Rs{i} = p.RsK{i}(p.bDOF{i},:);
    
    p.Nrbm_s(i) = size(p.Rs{i},2);
    p.Nrbm = p.Nrbm + p.Nrbm_s(i);
    if p.Nrbm_s(i)>0
        p.nFloating = p.nFloating + 1;
        p.R = [p.R zeros(size(p.R,1),size(p.Rs{i},2)); ...
            zeros(size(p.Rs{i},1),size(p.R,2)) p.Rs{i}];
        p.RK = [p.RK zeros(size(p.RK,1),size(p.RsK{i},2)); ...
            zeros(size(p.RsK{i},1),size(p.RK,2)) p.RsK{i}];
    else
        p.R = [p.R; ...
            zeros(size(p.SKs{i},1),size(p.R,2))];
        p.RK = [p.RK; ...
            zeros(size(p.Ks{i},1),size(p.RK,2))];
    end
end

if strcmp(p.mode,'dynamic')
    p.FItrans = p.B*(p.D\p.B'); % sum over substruct
elseif strcmp(p.mode,'static')
    p.FI = p.B*p.Kp*p.B';
    p.d = p.B*p.Kp*p.f;
    p.e = p.RK'*p.f;
end
% #Zeitmessung l528-639
time8=toc;
time_rel8=time8/(639-528);
%{
disp('Zeit l:2-178')
disp(time1)
disp(time_rel1)
disp('Zeit l:182-216')
disp(time2)
disp(time_rel2)
disp('Zeit l:220-245')
disp(time3)
disp(time_rel3)
disp('Zeit l:247-313')
disp(time4)
disp(time_rel4)
disp('Zeit l:315-376')
disp(time5)
disp(time_rel5)
disp('Zeit l:378-421')
disp(time6)
disp(time_rel6)
disp('Zeit l:423-526')
disp(time7)
disp(time_rel7)
disp('Zeit l:528-639')
disp(time8)
disp(time_rel8)
%}
if p.globalassembly==1
    p.Full = Full;
end
end