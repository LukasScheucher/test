function [p]=LMatrix(p)
%{
disp('Ndof_s:')
disp(p.Ndof_s)
disp(['NdofFull: ' num2str(p.NdofFull)])
disp(['Rank B: ' num2str(rank(p.B2))])
%}
p.L_man=zeros(size(p.B2,2),size(p.B2,2)-rank(p.B2));%sum(p.Ndof_s),p.NdofFull-rank(p.B));
%disp('L_man size:')
%disp(size(p.L_man))                   
tol=p.geom_tol;
%{
disp('B2-Matrix:')
disp(p.B2)
disp(size(p.B2))
%}
B2_help=p.B2;
masterdofs_subs=[];
slavedofs=[];
m=1;
n_s=1;
lm_conf=[];
n_set_conf=[];
lmc=1;
sign_change=1; 
if strcmp(p.mesh_method,'Mortar')
    sign_change=-1; % change algebraic sign for Mortar, as Master-values become negative instead of Slave-values for NTS
end
for l=1:size(p.B2,1)
    n_set=zeros(1,size(p.B2,2));
    i=1;
    for n=1:size(p.B2,2)
        if sign_change*p.B2(l,n)<-tol % find slave entries
            slavedofs(n_s)=n;
            n_s=n_s+1;
            n_set(i)=n;
            i=i+1;
        else
            if sign_change*p.B2(l,n)>tol % find master entries
                n_set(i)=n;
                i=i+1;
                b=find(n==masterdofs_subs);
                %c=find(n==slavedofs);
                if isempty(b) %&& isempty(c) 
                    masterdofs_subs(m)=n;
                    m=m+1;
                end
            end
        end
    end
    n=1;
    while n<=size(n_set,2) % delete zero-entries in n_set as they were never edited
        if n_set(1,n)==0
            n_set(:,n)=[];
            n=1;
        else
            n=n+1;
        end
    end
    %disp('n_set:')
    %disp(n_set)
    % Note conforming dofs in n_set_conf and corresponding LMs in
    % lm_conf
    if size(n_set,2)==2 && abs(abs(p.B2(l,n_set(1,1)))-abs(p.B2(l,n_set(1,2))))<=tol % condition for conforming dofs
        lm_conf(1,lmc)=l;
        for i=1:2
            n_set_conf(i,lmc)=n_set(1,i);
        end
        lmc=lmc+1;
    end            
end
%% Search for redundant LMs

disp('lm_conf:')
disp(lm_conf)

slavedofs=[1:size(B2_help,2)];

m=1;
zerocol=[];
for n=1:size(B2_help,2)
    if sum(abs(B2_help(:,n)))<=tol
        zerocol(m)=slavedofs(n);
        m=m+1;
    end
end

%{
disp('B_s:')
disp(B_s)
disp(size(B_s))
disp('slavedofs:')
disp(slavedofs)
disp('masterdofs:')
disp(masterdofs)
%}

masterdofs_subs=sort(masterdofs_subs);
%disp('masterdofs_subs')
%disp(masterdofs_subs)
masterdofs=[];
ma=1;
%% Filling the L-matrix
% Define global dofs for conforming dofs, which are part of redundant LMs
disp('defined conforming dofs:')
g_dof=1;

%mastdofs_subs2=masterdofs_subs;
%masterdofs_subs=[33,34,35,36,37,38,39,40];
disp('masterdofs_subs:')
disp(masterdofs_subs)
disp(size(masterdofs_subs))
masterdofs=masterdofs_subs;%(1,1:size(p.L_man,2)-size(zerocol,2));

%masterdofs=[73    74    77    78    79    80    83    84    85    86    87    88    89    90   151   152   153   154   155   156   169   170   173   174];
slavedofs=[1:size(B2_help,2)];
slavedofs([masterdofs,zerocol])=[]; 
%mdof_free=[size(p.L_man,2)-size(zerocol,2)+1:size(masterdofs_subs,2)];
%disp('mdof_free:')
%disp(mdof_free)
%disp('masterdofs')
%disp(masterdofs)
%disp(['Rang(B2_help): ' num2str(rank(B2_help))])
%disp(['Rang(B_s): ' num2str(rank(B2_help(:,slavedofs)))])

% Choose global dofs in "masterdofs"
if rank(B2_help)>rank(B2_help(:,slavedofs))
    masterdofs=masterdofs_subs;             % masterdofs_subs = all potential master-dofs with negative (Mortar) or positive (NTS) entries
    slavedofs=[1:size(B2_help,2)];
    slavedofs([masterdofs,zerocol])=[];     % zerocol = all dofs with zero entries in B-Matrix only
    B_s=B2_help(:,slavedofs);
    n_collision=[];                         % n_collision = all dofs to delete from "masterdofs"
    i=1;
    for l=1:size(B2_help,1)         % Search rows, which are zero, if all potential global dofs are deleted and save their slave-dofs in "n_collision"
        if sum(abs(B_s(l,:)))<=tol && sum(abs(B2_help(l,:)))>tol
            n_slave_col=find(sign_change*B2_help(l,:)<-tol);
            for j=1:size(n_slave_col,2)
                searchvar=find(n_slave_col(j)==n_collision);
                if isempty(searchvar)
                    n_collision=[n_collision,n_slave_col(j)];
                end
            end
            i=i+1;
        end
    end
    %disp('n_collision')
    %disp(n_collision)
    disp(size(n_collision))
    disp(size(masterdofs))
    disp(size(p.L_man,2)-size(zerocol,2))
    n_collision2=n_collision;
    for n=1:size(n_collision,2) % save dofs from "n_collision", if they are in "masterdofs" as well and replace "n_collision" with them
        i=find(n_collision(n)==masterdofs);
        n_collision2(1,n)=i;
    end
    n_collision=n_collision2;
    % (Note Andreas Seibold 18.02.2016): This worked only for the special
    % cases I've calculated so far. For a more general approach masterdofs
    % are deleted together with LMs below
    
    %{
    if (size(masterdofs,2)-size(n_collision,2))<size(p.L_man,2)-size(zerocol,2) % If "n_collision" is still to long, delete waste dofs
        n_collision=n_collision(1:size(masterdofs,2)-(size(p.L_man,2)-size(zerocol,2)));
    end
    %}
    masterdofs(n_collision)=[]; % delete dofs from "n_collision" from "masterdofs" to achieve global dofs in "masterdofs" only
    slavedofs=[1:size(B2_help,2)];
    slavedofs([masterdofs,zerocol])=[];     
end

disp('masterdofs:')
disp(masterdofs)
disp(['Rang(B2_help): ' num2str(rank(B2_help))])
disp(['Rang(B_s): ' num2str(rank(B2_help(:,slavedofs)))])
disp(size(B2_help(:,slavedofs)))
slavedofs=[1:size(B2_help,2)];
slavedofs([masterdofs,zerocol])=[];
disp('slavedofs:')
disp(slavedofs)
disp(size(slavedofs))

% Search for redundant LMs
lm_delete=[];   % Vektor of redundant or zero LMs, which have to be deleted
lm_del=1;
lm_zero=[];
z=1;
lm_set=[1:size(B2_help,1)]; % Set of all LMs
l=1;
N_check=2;
mdof_tolarge=0;     % "mdof_tolarge": 1="masterdofs" is to large; 2=an entry in "masterdofs" has just been deleted; 0="masterdofs" is equal to global dofs

while size(lm_set,2)>rank(B2_help(:,slavedofs)) && l<=size(lm_set,2) && N_check<=3 % There are linear dependent equations in "B2_help"
    if mdof_tolarge==2  % Reset the set of LMs befor starting with deleting LMs after waste masterdofs were deleted
        lm_set=[1:size(B2_help,1)]; % Set of all LMs
        lm_set(lm_zero)=[];
        l=1;
        N_check=2;
    end
    lm_redundant=0; 
    
    if size(masterdofs,2)>size(p.L_man,2)-size(zerocol,2) % Check if "masterdofs" is to large
        slavedofs=[1:size(B2_help,2)];
        slavedofs([masterdofs,zerocol])=[];
        mdof_tolarge=1;   
    else
        mdof_tolarge=0;
    end
    if sum(abs(B2_help(lm_set(l),slavedofs)))>tol   % If the current equation (here: compare-LM) is no zero LM
        for lm=1:size(lm_set,2)                         % Loop over all other LMs
            if lm_set(lm)~=lm_set(l) && sum(abs(B2_help(lm_set(lm),slavedofs)))>tol % The looped LM is not equal to the compare-LM and the equation is not zero
                if N_check==2   % Two equations are compared
                    lm_check=[lm_set(l);lm_set(lm)];    
                    B_check=B2_help(lm_check,slavedofs);    % Matrix of two equations with compare-LM and looped LM
                    if rank(B_check)<N_check    % If the rank of the Check-Matrix is lower than the number of compared equations, the equations are linear dependent and one LM is written in "lm_delete"
                        if mdof_tolarge==1 % If there are to many masterdofs left, delete one masterdof first, before you continue with the LMs
                            mdof=find(abs(B2_help(lm_set(l),:))>tol);
                            for m=1:size(mdof,2)
                                i=find(mdof(m)==masterdofs);
                                if isempty(i)==0
                                    %disp(['deleted masterdof: ' num2str(masterdofs(i))])
                                    masterdofs(i)=[];
                                    mdof_tolarge=2;
                                    break
                                end
                            end
                        else
                            %disp('B_check')
                            %disp(B_check)
                            %disp('lm_check')
                            %disp(lm_check)
                            lm_delete(lm_del)=lm_set(l);
                            lm_del=lm_del+1;
                            lm_set(l)=[];
                            lm_redundant=1;
                        end
                        break
                    end
                else    % Three equations are compared
                    for i=1:size(lm_set,2)  % Loop again over all other LMs
                        if lm_set(i)~=lm_set(l) && lm_set(i)~=lm_set(lm) && sum(abs(B2_help(lm_set(i),slavedofs)))>tol
                            lm_check=[lm_set(l);lm_set(lm);lm_set(i)];  % Analog to comparison of two equations, three equations are checked 
                            B_check=B2_help(lm_check,slavedofs);
                            if rank(B_check)<N_check
                                if mdof_tolarge==1
                                    mdof=find(abs(B2_help(lm_set(l),:))>tol);
                                    for m=1:size(mdof,2)
                                        j=find(mdof(m)==masterdofs);
                                        if isempty(j)==0
                                            %disp(['deleted masterdof: ' num2str(masterdofs(j))])
                                            masterdofs(j)=[];
                                            mdof_tolarge=2;
                                            break
                                        end
                                    end
                                else
                                    %disp('B_check')
                                    %disp(B_check)
                                    %disp('lm_check')
                                    %disp(lm_check)
                                    lm_delete(lm_del)=lm_set(l);
                                    lm_del=lm_del+1;
                                    lm_set(l)=[];
                                    lm_redundant=1;
                                end
                                break
                            end
                        end
                        if lm_redundant>0 || mdof_tolarge==2
                            break
                        end
                    end
                end
            end
            if lm_redundant>0 || mdof_tolarge==2
                break
            end
        end
    else    % For zero-rows:
        lm_delete(lm_del)=lm_set(l);
        lm_del=lm_del+1;
        lm_zero(z)=lm_set(l); % In lm_zero are all zero-rows saved, as lm_set is reset after masterdofs are deleted and must not be considered in the following loops
        z=z+1;
        lm_set(l)=[];
        lm_redundant=1;
    end
    if lm_redundant==1 || mdof_tolarge==2
        l=1;
    else
        l=l+1;  % If there are no redundant lms, or to many global dofs, go to further LMs
    end
    if l>size(lm_set,2)
        l=1;
        N_check=N_check+1; % If all LMs are checked, do the loop again for three equations
    end
end
disp('masterdofs:')
disp(masterdofs)
disp('lm_delete:')
disp(lm_delete)
%masterdofs=[33, 34, 35, 36, 37, 38, 39, 40, 65, 66];
slavedofs=[1:size(B2_help,2)];
slavedofs([masterdofs,zerocol])=[];
%disp(lm_delete)
%lm_delete=[9, 10, 27, 28];
%lm_delete=[5,6,7,8,13,14,15,16]; %NTS2x2
disp('abs Zeilensumme B2_help:')
rowsum=zeros(size(B2_help,1),1);
for l=1:size(B2_help,1)
    rowsum(l)=sum(abs(B2_help(l,:)));
end
disp(rowsum)

B2_help(lm_delete,:)=[];

for l=1:size(B2_help,1)
    B2_help(l,:)=B2_help(l,:)/norm(B2_help(l,:));
end

disp(['globale Interface-dofs: ' num2str(size(p.L_man,2)-size(zerocol,2))])
disp(['globale Nicht-Interface dofs: ' num2str(size(zerocol,2))])
disp(['Rang(B2_help): ' num2str(rank(B2_help))])
disp(['Rang(B2_help)^T: ' num2str(rank(B2_help'))])
disp(['Rang(p.B2): ' num2str(rank(p.B2))])
disp(['Rang(B_s): ' num2str(rank(B2_help(:,slavedofs)))])
disp(size(B2_help(:,slavedofs)))
disp(['det(B_s): ' num2str(det(B2_help(:,slavedofs)))])
p.L_man=zeros(size(p.B2,2),size(p.B2,2)-rank(B2_help));

disp('#############')
disp('slavedofs:')
disp(slavedofs)
disp('masterdofs:')
disp(masterdofs)

masterdofs_subs=masterdofs;
disp('Beschriebene dofs:')
for m=1:size(masterdofs_subs,2)
    disp(masterdofs_subs(m))
    p.L_man(masterdofs_subs(m),g_dof)=1;
    g_dof=g_dof+1;
end
disp(['g_dof: ' num2str(g_dof-1)])

%% Split B-matrix in B_s and B_m for local and global dofs


disp('slavedofs:')
disp(slavedofs)

B_m=B2_help(:,masterdofs);
B_s=B2_help(:,slavedofs);

disp('B_s:')
disp(B_s)
disp(size(B_s))

disp('Determinant of B_s:')
disp(det(B_s))
L=-inv(B_s)*B_m;
disp('B_s invers:')
disp(inv(B_s))
disp('Hilfsmatrix L:')
disp(L)
disp(size(L))


for m=1:size(L,2) % go through interface-slave dofs = number of rows of L
    g=find(p.L_man(masterdofs(m),:)==1);
    if isempty(g)
        g=find(abs(p.L_man(masterdofs(m),:))>tol);
    end
    for n=1:size(L,1) % amount of columns of L = amount of master dofs
        if abs(L(n,m))>tol
            p.L_man(slavedofs(n),g)=L(n,m);%*p.L_man(masterdofs(m),g);
        end
    end
end
counter=0;
for i=1:size(p.L_man,1)
    if sum(p.L_man(i,:))==0
        counter=counter+1;
    end
end
disp(['Nullzeilen: ' num2str(counter)])
disp(size(p.L_man))
disp(['uebrige globale dofs: ' num2str(size(p.L_man,2)-g_dof+1)])
%%Translate remaining local dofs into global dofs
for n=1:size(p.L_man,1)
    if sum(abs(p.L_man(n,:)))<=tol && g_dof<=size(p.L_man,2)
        p.L_man(n,g_dof)=1;
        g_dof=g_dof+1;
    end
end

disp('L_man:')
disp(p.L_man)
disp('Zeilensumme von L:')
zeilensum=zeros(1,size(p.L_man,1));
for i=1:size(p.L_man,1)
    zeilensum(i)=sum(p.L_man(i,:));
end
disp(zeilensum)
end