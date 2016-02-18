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
masterdofs=masterdofs_subs(1,1:size(p.L_man,2)-size(zerocol,2));

%masterdofs=[73    74    77    78    79    80    83    84    85    86    87    88    89    90   151   152   153   154   155   156   169   170   173   174];
slavedofs=[1:size(B2_help,2)];
slavedofs([masterdofs,zerocol])=[]; 
m=1;
mdof_free=[size(p.L_man,2)-size(zerocol,2)+1:size(masterdofs_subs,2)];
disp('mdof_free:')
disp(mdof_free)
disp('masterdofs')
disp(masterdofs)
disp(['Rang(B2_help): ' num2str(rank(B2_help))])
disp(['Rang(B_s): ' num2str(rank(B2_help(:,slavedofs)))])
if rank(B2_help)>rank(B2_help(:,slavedofs)) && m<=size(masterdofs,2)/2
    masterdofs=masterdofs_subs;
    slavedofs=[1:size(B2_help,2)];
    slavedofs([masterdofs,zerocol])=[]; 
    B_s=B2_help(:,slavedofs);
    lm_collision=[];
    n_collision=[];
    i=1;
    for l=1:size(B2_help,1)
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
    disp('n_collision')
    disp(n_collision)
    disp(size(n_collision))
    disp(size(masterdofs))
    disp(size(p.L_man,2)-size(zerocol,2))
    n_collision2=n_collision;
    for n=1:size(n_collision,2)
        i=find(n_collision(n)==masterdofs);
        n_collision2(1,n)=i;
    end
    n_collision=n_collision2;
    if (size(masterdofs,2)-size(n_collision,2))<size(p.L_man,2)-size(zerocol,2)
        n_collision=n_collision(1:size(masterdofs,2)-(size(p.L_man,2)-size(zerocol,2)));
    end
    masterdofs(n_collision)=[];
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

lm_delete=[];
lm_del=1;
lm_set=[1:size(B2_help,1)];
l=1;
N_check=2;
while size(lm_set,2)>rank(B2_help(:,slavedofs)) && l<=size(lm_set,2) && N_check<=3
    lm_redundant=0;
    if sum(abs(B2_help(lm_set(l),slavedofs)))>tol
        for lm=1:size(lm_set,2)
            if lm_set(lm)~=lm_set(l) && sum(abs(B2_help(lm_set(lm),slavedofs)))>tol
                if N_check==2
                    lm_check=[lm_set(l);lm_set(lm)];
                    B_check=B2_help(lm_check,slavedofs);
                    if rank(B_check)<N_check
                        disp('B_check')
                        disp(B_check)
                        disp('lm_check')
                        disp(lm_check)
                        lm_delete(lm_del)=lm_set(l);
                        lm_del=lm_del+1;
                        lm_set(l)=[];
                        lm_redundant=1;
                        break
                    end
                else
                    for i=1:size(lm_set,2)
                        if lm_set(i)~=lm_set(l) && lm_set(i)~=lm_set(lm) && sum(abs(B2_help(lm_set(i),slavedofs)))>tol
                            lm_check=[lm_set(l);lm_set(lm);lm_set(i)];
                            B_check=B2_help(lm_check,slavedofs);
                            if rank(B_check)<N_check
                                disp('B_check')
                                disp(B_check)
                                disp('lm_check')
                                disp(lm_check)
                                lm_delete(lm_del)=lm_set(l);
                                lm_del=lm_del+1;
                                lm_set(l)=[];
                                lm_redundant=1;
                                break
                            end
                        end
                        if lm_redundant>0
                            break
                        end
                    end
                end
            end
            if lm_redundant>0
                break
            end
        end
    else
        lm_delete(lm_del)=lm_set(l);
        lm_del=lm_del+1;
        lm_set(l)=[];
        lm_redundant=1;
    end
    if lm_redundant==1
        l=1;
    else
        l=l+1;
    end
    if l>size(lm_set,2)
        l=1;
        N_check=N_check+1;
    end
end

disp('lm_delete:')
disp(lm_delete)

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
disp('Rang:')
disp(rank(B_s))
disp(rank(B2_help))
disp(rank(p.B2))
%disp('B_m:')
%disp(B_m)
%disp(size(B_m))
disp('Determinante of B_s:')
disp(det(B_s))
L=-inv(B_s)*B_m;
disp('Hilfsmatrix L:')
disp(L)
disp(size(L))
disp(size(masterdofs))
disp(size(slavedofs))
disp(['g_dof: ' num2str(g_dof-1)])

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