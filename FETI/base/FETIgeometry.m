function [p] = FETIgeometry(p)
% there are two orderings: ordering of the dofs u, and ordering of the
% lagrange multipliers lambda (the two fields of this hybrid two-field
% method)
p.Nscorner=zeros(1,4); % Matrix for indices of substructures on domain's corners in mathematical positiv direction beginning at bottom left corner

if p.nonconforming==1
    p.Ns = size(p.sizes,2);
    
    %% Control of geometry input parameters 
    if size(p.sizes,2)~=size(p.elcount,2)
        error('Warning: p.elcount must have the same number of columns as p.sizes')
    end
    
    
    %% Position vectors at left lower corners of substructures
    len_struc=0; % current length of assembled structure
    length_top=0;
    height=0;
    sub=1;
    p.poss=zeros(2,p.Ns); % Position of bottom left corner for every substructure
    p.cornerset=[];
    i_below=1;
    tol=p.geom_tol;
    p.Nsy=1; %### For similar shapes of substructures
    p.Nscorner(1)=1;
    p.Nscorner(3)=p.Ns;
    % Write bottom left corner positions in matrix p.poss
    p.right_side=[];
    i=1;
    for s = 1:p.Ns % Build undermost substructures 
        p.poss(1,s)=len_struc;
        p.poss(2,s)=height;
        % Current top corners of substructures to fit the following substructures on top
        p.cornerset(1,2*s-1)=len_struc; % x-coordinate of left top corner 
        p.cornerset(2,2*s-1)=p.poss(2,s)+p.sizes(2,s); % y-coordinate
        p.cornerset(1,2*s)=p.poss(1,s)+p.sizes(1,s); % x-coordinate of right top corner
        p.cornerset(2,2*s)=p.poss(2,s)+p.sizes(2,s);
        
        len_struc=len_struc+p.sizes(1,s);
        disp(['Bottom: ' num2str(s) ' Length: ' num2str(len_struc)])
        if abs(len_struc-p.Length)<=tol || (len_struc-p.Length)>tol
            if (len_struc-p.Length)>tol
                p.sizes(1,s)=p.sizes(1,s)+(p.Length-len_struc);
                display(['Warning: Substructure ' num2str(s) ' exceeded Length of Cantilever with tollerance ' num2str(tol) ' and has been shortened to ' num2str(p.sizes(1,s)) '!'])
            end
            p.Nscorner(2)=s;
            len_struc=0;
            p.right_side(i)=s;
            p.Nsx=s; %### For similar shapes of substructures
            break; 
        end
    end
    for k =1:p.Ns % Build topmost substructures to get the substructures with highest index at the top of the cantilever
        inv_counter=p.Ns+1-k;
        if p.poss(1,inv_counter)==0 && p.poss(2,inv_counter)==0 && inv_counter~=1
            p.poss(1,inv_counter)=p.Length-p.sizes(1,inv_counter)-length_top;
            p.poss(2,inv_counter)=p.Height-p.sizes(2,inv_counter);
            length_top_backup=length_top;
            length_top=length_top+p.sizes(1,inv_counter);
            disp(['Top: ' num2str(inv_counter) ' Length: ' num2str(length_top)])
            if abs(length_top-p.Length)<=tol || (length_top-p.Length)>tol
                if (length_top-p.Length)>tol
                    p.sizes(1,inv_counter)=p.sizes(1,inv_counter)+(p.Length-length_top);
                    p.poss(1,inv_counter)=p.sizes(1,inv_counter)+length_top_backup;
                    display(['Warning: Substructure ' num2str(inv_counter) ' exceeded Length of Cantilever with tollerance ' tol ' and has been shortened to ' num2str(p.sizes(1,inv_counter)) '!'])
                end
                p.Nscorner(4)=inv_counter;
                break;
            end
        end
    end
    disp('p.poss mit oberster und unterster Reihe:')
    disp(p.poss)
    disp(p.right_side)
    disp(inv_counter)
    k=1;
    for i =p.right_side(1)+1:inv_counter-1 % Build rest of substructures from left to right and from bottom to top
        if p.poss(1,i)==0 && p.poss(2,i)==0 && i~=1
            sub=sub+1;
            [i_below2,i_below3]=min(p.cornerset(2,:)); % minimal y-coordinate
            height=p.cornerset(2,i_below3);
            disp('p.cornerset')
            disp(p.cornerset)
            len_struc=p.cornerset(1,i_below3);
            p.poss(1,i)=len_struc;
            p.poss(2,i)=height;
            
            len_struc=len_struc+p.sizes(1,i);
            disp(['Rest: ' num2str(i) ' Length: ' num2str(len_struc)])
            
            if i_below+2<size(p.cornerset,2) % Check, that the last substructure is not the corresponding one
                if p.cornerset(2,i_below+1)~=p.cornerset(2,i_below+2) && abs(len_struc-p.cornerset(1,i_below+1))<=tol
                    i_below=i_below+2;
                end
            end
            
            if abs(len_struc-p.Length)<=tol || (len_struc-p.Length)>tol % Substructure at the cantilever's end reached
                if (len_struc-p.Length)>tol
                    p.sizes(1,i)=p.sizes(1,i)+(p.Length-len_struc);
                    display(['Warning: Substructure ' num2str(i) ' exceeded Length of Cantilever with tollerance ' tol ' and has been shortened to ' num2str(p.sizes(1,i)) '!'])
                end
                len_struc=0;
                sub=1;
                height=height+p.sizes(2,i_below);
                p.Nsy=2+p.Nsy; % For similar shapes of substructures
                k=k+1;
                p.right_side(k)=i;
            end

            % Current top corners of substructures to fit the following substructures on top
            p.cornerset(1,2*sub-1)=p.poss(1,i); % x-coordinate of left top corner 
            p.cornerset(2,2*sub-1)=p.poss(2,i)+p.sizes(2,i); % y-coordinate
            p.cornerset(1,2*sub)=p.poss(1,i)+p.sizes(1,i);
            p.cornerset(2,2*sub)=p.poss(2,i)+p.sizes(2,i);
        end
    end
    p.right_side(k+1)=p.Ns;
    disp('Substructure positions: ')
    disp(p.poss)
    disp(['Geometric tolerance: ' num2str(tol)])
    if size(p.elcount,2)~=p.Ns
        i=size(p.elcount,2);
        while size(p.elcount,2)~=p.Ns
            p.elcount=horzcat(p.elcount,p.elcount(i));
            i=i+1;
        end
        disp('Count of elementsizes expanded with last value to match count of substructures.')
    end
    
    % Element and node generation
    Ndof_n=2;
    p.Nely=p.elcount;
    p.Nelx=p.Nely;
    p.elHeight=p.elcount;
    p.Nn_s_nc=p.elHeight; % Number of nodes in one substructure
    p.Nn_s_ncFull=0;
    for s=1:p.Ns
        p.elHeight(s)=p.sizes(2,s)/p.Nely(s);
        p.Nelx(s)=int32(p.sizes(1,s)/(p.sizes(2,s)/p.Nely(s)));
        if abs(p.Nelx(s)-p.sizes(1,s)/p.elHeight(s))<=tol % Check, if element-size matches the susbtructure's x-length
            p.Nn_s_nc(s)=(p.Nelx(s)+1)*(p.Nely(s)+1); % numer of nodes in each substructure for quadratic elements
            p.Nn_s_ncFull=p.Nn_s_ncFull+p.Nn_s_nc(s);
        else
            disp(['Substructure ' num2str(i) 's element-size does not match its length. Process canceled!'])
            break
        end 
    end
    disp('p.sizes(2): ') 
    disp(p.sizes(2,:))
    disp('p.Nely(s): ')
    disp(p.Nely)
    disp('p.elHeight: ') 
    disp(p.elHeight)
    p.Ndof_s = p.Nn_s_nc*Ndof_n; % dof per substructure
    p.Nn_s=p.Nn_s_nc;
    
    % If there's no relative movement between substructures: determine starting positions of interface nodes to detect conforming
    % nodes => those nodes can be calculated as in FETI with conforming meshes
    %n_int_full=0;
    %N_int=0;
    p.Nlm=0;
    int_full=p.Ns*p.Ns;
    int=1;
    p.int=zeros(4,int_full);
    for s = 1:p.Ns
        z=s;
        while z<=p.Ns
            if abs(p.poss(1,s)+p.sizes(1,s)-p.poss(1,z))<=tol && p.poss(2,s)+p.sizes(2,s)-p.poss(2,z)>0 && p.poss(2,z)+p.sizes(2,z)-p.poss(2,s)>0 % Condition for interface on eastern side
                p.int(1,int)=s;
                p.int(2,int)=z;
                p.int(3,int)=1;         % Horizontal (2) or vertical (1) interface or cross connection (3)
                if p.poss(2,z)>=p.poss(2,s)
                    if p.poss(2,z)+p.sizes(2,z)>=p.poss(2,s)+p.sizes(2,s)
                        int_length=p.poss(2,s)+p.sizes(2,s)-p.poss(2,z);
                    else
                        int_length=p.sizes(2,z);
                    end
                else
                    if p.poss(2,s)+p.sizes(2,s)>=p.poss(2,z)+p.sizes(2,z)
                        int_length=p.poss(2,z)+p.sizes(2,z)-p.poss(2,s);
                    else
                        int_length=p.sizes(2,s);
                    end
                end
                disp(['Interface Length: ' num2str(int_length)])
                disp(ceil(int_length/p.elHeight(s)))
                if p.elHeight(s)<p.elHeight(z) % definition of master(1), slave(2), conforming interface (3) and cross connection (4)
                    p.int(4,int)=2;
                    p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(s)))+1;%+p.Nely(s)+1;
                else
                    if p.elHeight(s)>p.elHeight(z)
                        p.int(4,int)=1;
                        p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(z)))+1;%+p.Nely(z)+1;
                    else                    % element sizes are equal
                        if abs(modfl(p.poss(2,s)-p.poss(2,z),p.elHeight(s),tol))<=tol % if the offset between interface-meshes is equal to a multiple of the elementheight
                            p.int(4,int)=3;
                        else
                            p.int(4,int)=2;
                        end
                        p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(s)))+1;%+p.Nely(s)+1;
                    end
                end
                disp(['p.Nlm nach int=' num2str(int) '(int_dir: ' num2str(p.int(3,int)) '): ' num2str(p.Nlm)])
                int=int+1;
            else
                if abs(p.poss(2,s)+p.sizes(2,s)-p.poss(2,z))<=tol && p.poss(1,s)+p.sizes(1,s)-p.poss(1,z)>0 && p.poss(1,z)+p.sizes(1,z)-p.poss(1,s)>0 % Condition for interface on northern side
                    p.int(1,int)=s;
                    p.int(2,int)=z;
                    p.int(3,int)=2;
                    if p.poss(1,z)>=p.poss(1,s)
                        if p.poss(1,z)+p.sizes(1,z)>=p.poss(1,s)+p.sizes(1,s)
                            int_length=p.poss(1,s)+p.sizes(1,s)-p.poss(1,z);
                        else
                            int_length=p.sizes(1,z);
                        end
                    else
                        if p.poss(1,s)+p.sizes(1,s)>=p.poss(1,z)+p.sizes(1,z)
                            int_length=p.poss(1,z)+p.sizes(1,z)-p.poss(1,s);
                        else
                            int_length=p.sizes(1,s);
                        end
                    end
                    disp(['Interface Length: ' num2str(int_length)])
                    disp(['s: ' num2str(ceil(int_length/p.elHeight(s)))])
                    disp(['z: ' num2str(ceil(int_length/p.elHeight(z)))])
                    if p.elHeight(s)<p.elHeight(z)
                        p.int(4,int)=2;
                        p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(s)))+1;%+p.Nelx(s)+1;
                    else
                        if p.elHeight(s)>p.elHeight(z)
                            p.int(4,int)=1;
                            p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(z)))+1;%+p.Nelx(z)+1;
                        else
                            if abs(modfl(p.poss(1,s)-p.poss(1,z),p.elHeight(s),tol))<=tol
                                p.int(4,int)=3;
                            else
                                p.int(4,int)=2;
                            end
                            p.Nlm=p.Nlm+int32(ceil(int_length/p.elHeight(s)))+1;%+p.Nelx(s)+1;
                        end
                    end
                    disp(['p.Nlm nach int=' num2str(int) '(int_dir: ' num2str(p.int(3,int)) '): ' num2str(p.Nlm)])
                    int=int+1;
                else
                    if abs(p.poss(1,s)+p.sizes(1,s)-p.poss(1,z))<=tol && abs(p.poss(2,s)+p.sizes(2,s)-p.poss(2,z))<=tol % Condition for cross connection
                        p.int(1,int)=s;
                        p.int(2,int)=z;
                        p.int(3,int)=3;
                        p.int(4,int)=4;
                        p.Nlm=p.Nlm+2;
                        int=int+1;
                    else
                        p.int(:,int)=[];
                    end
                end
            end
            z=z+1;
        end
    end
    
    int_full=int-1;
    %disp(['int_full: ' num2str(int_full)])
    %disp(['Nlm: ' num2str(p.Nlm)])
    
    %p.posn_int=zeros(5,p.Nlm*2);
    N_s_self=zeros(1,int_full);
    
    for int = 1:int_full % Nodes related to interfaces with node positions
        p.posn_int{int}=zeros(4,p.Nlm*2);
        s=p.int(1,int);
        z=p.int(2,int);
        if p.int(3,int)==1 % vertical interface
            N_s_self(int)=p.Nely(s)+1;
            for n = 1:p.Nely(s)+p.Nely(z)+2
                if n<=N_s_self(int)
                    p.posn_int{int}(1,n)=p.poss(1,s)+p.sizes(1,s);
                    p.posn_int{int}(2,n)=p.poss(2,s)+p.elHeight(s)*(n-1);
                    p.posn_int{int}(3,n)=(p.Nelx(s)+1)*n;  % local node number
                    p.posn_int{int}(4,n)=s;
                else
                    p.posn_int{int}(1,n)=p.poss(1,z);
                    p.posn_int{int}(2,n)=p.poss(2,z)+p.elHeight(z)*(n-N_s_self(int)-1);
                    p.posn_int{int}(3,n)=(p.Nelx(z)+1)*(n-N_s_self(int)-1)+1;% local node number
                    p.posn_int{int}(4,n)=z;
                end
            end
        else
            if p.int(3,int)==2 % horizontal interface
                N_s_self(int)=p.Nelx(s)+1;
                for n = 1:p.Nelx(s)+p.Nelx(z)+2
                    if n<=N_s_self(int)
                        p.posn_int{int}(1,n)=p.poss(1,s)+p.elHeight(s)*(n-1);
                        p.posn_int{int}(2,n)=p.poss(2,s)+p.sizes(2,s);
                        p.posn_int{int}(3,n)=(p.Nelx(s)+1)*p.Nely(s)+n;
                        p.posn_int{int}(4,n)=s;
                    else
                        p.posn_int{int}(1,n)=p.poss(1,z)+p.elHeight(z)*(n-N_s_self(int)-1);
                        p.posn_int{int}(2,n)=p.poss(2,z);
                        p.posn_int{int}(3,n)=n-N_s_self(int);
                        p.posn_int{int}(4,n)=z;
                    end
                end
            else % cross connection
                s_cross=s;
                z_cross=z;
                for n=1:2
                    s=s_cross+(n-1);
                    z=z_cross-(2-n);
                    p.posn_int{int}(1,n)=p.poss(1,s)+p.sizes(1,s)*(2-n);
                    p.posn_int{int}(2,n)=p.poss(2,s)+p.sizes(2,s);
                    p.posn_int{int}(3,n)=(p.Nelx(s)+1)*(p.Nely(s)+(2-n))+(n-1);
                    p.posn_int{int}(4,n)=s;
                    p.posn_int{int}(1,n+2)=p.poss(1,z)+p.sizes(1,z)*(2-n);
                    p.posn_int{int}(2,n+2)=p.poss(2,z);
                    p.posn_int{int}(3,n+2)=(p.Nelx(z)*(2-n)+1);
                    p.posn_int{int}(4,n+2)=z;
                end
            end
        end 
        for n=1:2*p.Nlm % Delete unnecessary columns
            if p.posn_int{int}(3,2*p.Nlm-n+1)==0
                p.posn_int{int}(:,2*p.Nlm-n+1)=[];
            end
        end
    end
    p.Nlm=Ndof_n*(p.Nlm+1); % 2 dof per node
    disp(['Anzahl Lagrange-Multiplikatoren: ' num2str(p.Nlm)])
    for s = 1:p.Ns
        % set NBC
        p.interfaceDOF{s}(1:p.Ndof_s(s),1) = 0;
        interfaceDOFlm{s}(1:p.Ndof_s(s),1) = 0;
        valenceNodes{s}(1:p.Ndof_s(s)) = 2;
    end
else
    p.Ns = p.Nsx*p.Nsy;
    p.Nn_s = (p.Nelx+1)*(p.Nely+1); % number of nodes in one element-assembled substructure
    Ndof_n = 2; % Ndof_n: dof per node, x-/y-transl.
    p.Ndof_s(1:p.Ns) = p.Nn_s*Ndof_n; % dof per substructure
    p.NdofFull = (p.Nsx*p.Nelx+1)*(p.Nsy*p.Nely+1)*2;
    
    disp('Ndof_s:')
    disp(p.Ndof_s)
    disp(p.NdofFull)
    
    p.Nelx_backup=p.Nelx;
    p.Nely_backup=p.Nely;
    p.Nelx=p.Ndof_s;
    p.Nely=p.Nelx;
    
    for s=1:p.Ns
        p.Nely(s) = p.Nely_backup;
        p.Nelx(s) = p.Nelx_backup;
    end   
    % global node ordering: all substructures sequentially

    % p.Nlm = (Nr. of horizontal connections + Nr. of vertical connections + Nr.
    % of cross connections)*2dofpercon
    p.Nlm = 2*(...
        (p.Nsx-1)*(p.Nely(1)+1)*p.Nsy ...
        + (p.Nsy-1)*(p.Nelx(1)+1)*p.Nsx ...
        + (p.Nsx-1)*(p.Nsy-1)*2 ...
        );

    % compute the element connectivity matrix p.ec
    s=1;
end
if p.nonconforming==1
    s_end=p.Ns;
else
    s_end=1;
end

for s=1:p.Ns
    p.ec{s} = zeros(p.Nelx(s)*p.Nely(s),4); % 4 columns for 4 nodes of the bilinear quadrilateral elements
    for r = 1:p.Nely(s) % r: row
        for c = 1:p.Nelx(s) % c: column
            i = (r-1)*p.Nelx(s)+c; % i: element number
            % p.ec(element,localnode) = globalnode
            p.ec{s}(i,1) = (r-1)*(p.Nelx(s)+1)+c;
            p.ec{s}(i,2) = (r-1)*(p.Nelx(s)+1)+c+1;
            p.ec{s}(i,3) = r*(p.Nelx(s)+1)+c+1;
            p.ec{s}(i,4) = r*(p.Nelx(s)+1)+c;
        end
    end
end

if p.nonconforming==1
    for s=1:p.Ns
        p.Bs{s} = zeros(p.Nlm,2*((p.Nelx(s)+1)*(p.Nely(s)+1)));
    end
else
    for s = 1:p.Ns
        % set NBC
        p.interfaceDOF{s}(1:p.Ndof_s(s),1) = 0;
        interfaceDOFlm{s}(1:p.Ndof_s(s),1) = 0;
        valenceNodes{s}(1:p.Ndof_s(s)) = 2;
        p.Bs{s} = zeros(p.Nlm,p.Ndof_s(s));
        betaScales{s}(1:p.Nlm,1) = 0;
    end
end

%% set up B-matrices

% identify node pairings = lagrange multipliers between a pair of nodes
% travel through all substructures row-wise from left to right (x++)
% and from the bottom row to the uppermost row (y++)
if p.nonconforming==0
    ipairings = 0;
    CornerPairings = [];

    for iNsy = 1:p.Nsy
        for iNsx = 1:p.Nsx
            s=iNsx+(iNsy-1)*p.Nsx;
            % calculate the index s of the substructure in row iNsy and col
            % iNsx
            iNs_self = (iNsy-1)*p.Nsx + iNsx;

            % calculate the indices of the adjacent substructures
            iNs_right = iNs_self + 1;
            iNs_above = iNs_self + p.Nsx;
            iNs_below = iNs_self - p.Nsx;
            iNs_rightabove = iNs_self + p.Nsx + 1;
            iNs_rightbelow = iNs_self - p.Nsx + 1;

            if iNsx < p.Nsx
                % this (self) is not a rightmost substructure, so setup
                % constraints with the substructure to the right
                for nodeRow = 1:(p.Nely(s)+1)
                    node_self = nodeRow*(p.Nelx(s)+1);
                    node_right = (nodeRow-1)*(p.Nelx(s)+1)+1;
                    ipairings = ipairings + 1;
                    nodepairings(ipairings,1:4) = [iNs_self iNs_right node_self node_right];
                    betaScalePre{iNs_self}{node_self} = [iNs_right node_right];
                    betaScalePre{iNs_right}{node_right} = [iNs_self node_self];

                    if nodeRow == 1 || nodeRow == (p.Nely(s)+1)
                        CornerPairings(end+1) = ipairings;
                    end
                end
            end
            if iNsy < p.Nsy
                % this (self) is not an uppermost substructure, so setup
                % constraints with the substructure above
                for nodeCol = 1:(p.Nelx(s)+1)
                    node_self = (p.Nelx(s)+1)*p.Nely(s) + nodeCol;
                    node_above = nodeCol;
                    ipairings = ipairings + 1;
                    nodepairings(ipairings,1:4) = [iNs_self iNs_above node_self node_above];
                    betaScalePre{iNs_self}{node_self} = [iNs_above node_above];
                    betaScalePre{iNs_above}{node_above} = [iNs_self node_self];

                    if nodeCol == 1 || nodeCol == (p.Nelx(s)+1)
                        CornerPairings(end+1) = ipairings;
                    end
                end
            end
            if iNsy < p.Nsy && iNsx < p.Nsx
                % this (self) is neither a rightmost nor an uppermost
                % substructure so we need a cross LM to the substructure right
                % above

                node_self = (p.Nelx(s)+1)*(p.Nely(s)+1);
                node_right = p.Nely(s)*(p.Nelx(s)+1)+1;
                node_above = p.Nelx(s)+1;
                node_rightabove = 1;
                ipairings = ipairings + 1;
                nodepairings(ipairings,1:4) = [iNs_self iNs_rightabove node_self node_rightabove];
                % the valence of the node "node_self" of "iNs_self" is 4
                % the valence of the node "node_rightabove" of "iNs_rightabove" is 4
                valenceNodes{iNs_self}(node_self) = 4;
                valenceNodes{iNs_rightabove}(node_rightabove) = 4;
                betaScalePre{iNs_self}{node_self} = [iNs_above node_above; iNs_right node_right; iNs_rightabove node_rightabove];
                betaScalePre{iNs_above}{node_above} = [iNs_self node_self; iNs_right node_right; iNs_rightabove node_rightabove];
                betaScalePre{iNs_right}{node_right} = [iNs_above node_above; iNs_self node_self; iNs_rightabove node_rightabove];
                betaScalePre{iNs_rightabove}{node_rightabove} = [iNs_above node_above; iNs_self node_self; iNs_right node_right];

                CornerPairings(end+1) = ipairings;
            end
            if iNsy > 1 && iNsx < p.Nsx
                % this (self) is neither a rightmost nor a lowermost
                % substructure so we need a cross LM to the substructure right
                % below

                node_self = p.Nelx(s)+1;
                node_right = 1;
                node_below = (p.Nely(s)+1)*(p.Nelx(s)+1);
                node_rightbelow = (p.Nelx(s)+1)*p.Nely(s) + 1;
                ipairings = ipairings + 1;
                nodepairings(ipairings,1:4) = [iNs_self iNs_rightbelow node_self node_rightbelow];
                % the valence of the node "node_self" of "iNs_self" is 4
                % the valence of the node "node_rightabove" of "iNs_rightabove" is 4
                valenceNodes{iNs_self}(node_self) = 4;
                valenceNodes{iNs_rightbelow}(node_rightbelow) = 4;
                betaScalePre{iNs_self}{node_self} = [iNs_below node_below; iNs_right node_right; iNs_rightbelow node_rightbelow];
                betaScalePre{iNs_below}{node_below} = [iNs_self node_self; iNs_right node_right; iNs_rightbelow node_rightbelow];
                betaScalePre{iNs_right}{node_right} = [iNs_below node_below; iNs_self node_self; iNs_rightbelow node_rightbelow];
                betaScalePre{iNs_rightbelow}{node_rightbelow} = [iNs_below node_below; iNs_self node_self; iNs_right node_right];

                CornerPairings(end+1) = ipairings;
            end
        end
    end
    p.CornerLMs = zeros(p.Nlm,1);
end
if p.nonconforming==1
    lm=1;
    if strcmp(p.mesh_method,'Mortar')==1
        slave_substructures=zeros(1,int_full);
        lm_backup=zeros(int_full,4);
    end
    for int = 1:int_full
        disp(['#Interface: ' num2str(int)])
        if p.int(3,int)==3 % cross nodes
            subs1=p.posn_int{int}(4,1);
            subs2=p.posn_int{int}(4,2);
            subs3=p.posn_int{int}(4,3);
            subs4=p.posn_int{int}(4,4);

            p.Bs{subs1}(lm+1,p.posn_int{int}(3,1)*2)=1;
            p.Bs{subs1}(lm,p.posn_int{int}(3,1)*2-1)=1;
            p.Bs{subs4}(lm+1,p.posn_int{int}(3,4)*2)=-1;
            p.Bs{subs4}(lm,p.posn_int{int}(3,4)*2-1)=-1;
            lm=lm+2;
            p.interfaceDOF{subs1}(p.posn_int{int}(3,1)*2,1)=1;
            p.interfaceDOF{subs1}(p.posn_int{int}(3,1)*2-1,1)=1;
            p.interfaceDOF{subs4}(p.posn_int{int}(3,4)*2,1)=1;
            p.interfaceDOF{subs4}(p.posn_int{int}(3,4)*2-1,1)=1;

            p.Bs{subs2}(lm+1,p.posn_int{int}(3,2)*2)=1;
            p.Bs{subs2}(lm,p.posn_int{int}(3,2)*2-1)=1;
            p.Bs{subs3}(lm+1,p.posn_int{int}(3,3)*2)=-1;
            p.Bs{subs3}(lm,p.posn_int{int}(3,3)*2-1)=-1;         
            lm=lm+2;
            p.interfaceDOF{subs2}(p.posn_int{int}(3,2)*2,1)=1;
            p.interfaceDOF{subs2}(p.posn_int{int}(3,2)*2-1,1)=1;
            p.interfaceDOF{subs3}(p.posn_int{int}(3,3)*2,1)=1;
            p.interfaceDOF{subs3}(p.posn_int{int}(3,3)*2-1,1)=1;   

            p.Nn_s_ncFull=p.Nn_s_ncFull+1; %Compensate the nodes detected as conforming
        else
            n_convert=0;
            N_int=size(p.posn_int{int},2);
            if p.int(4,int)==1
                n_convert=1; % switch from master-side to slave-side
            end

            switch p.mesh_method
                case 'NTS-LM'
                    disp('Slave-Nodes:')
                    disp([num2str(1+n_convert*N_s_self(int)) ' bis ' num2str(N_s_self(int)+n_convert*(N_int-N_s_self(int)))])
                    disp('Master-Nodes:')
                    disp([num2str(1+(1-n_convert)*N_s_self(int)) ' bis ' num2str(N_int-n_convert*(N_int-N_s_self(int)))])
                    for n=1+n_convert*N_s_self(int):N_s_self(int)+n_convert*(N_int-N_s_self(int)) % jump through slave-nodes
                        n_conf=0;
                        p.connection=zeros(2,N_int);
                        for m=1+(1-n_convert)*N_s_self(int):N_int-n_convert*(N_int-N_s_self(int)) % jump through master-nodes
                            if abs(p.posn_int{int}(1,n)-p.posn_int{int}(1,m))<=tol && abs(p.posn_int{int}(2,n)-p.posn_int{int}(2,m))<=tol
                                n_conf=m; % conform
                            else
                                p.connection(1,m)=m;
                                p.connection(2,m)=sqrt((p.posn_int{int}(1,n)-p.posn_int{int}(1,m))^2+(p.posn_int{int}(2,n)-p.posn_int{int}(2,m))^2);
                            end
                        end
                        if n_conf==0
                            i=1;
                            while i<=size(p.connection,2)
                                if p.connection(1,i)==0
                                    p.connection(:,i)=[];   % Not overwritten columns are deleted, to avoid zero-columns being detected as minimum
                                    i=i-1;
                                end
                                i=i+1;
                            end
                            [ratio,m_C]=min(p.connection(2,:));
                            m_C=int32(p.connection(1,m_C));

                            if p.int(3,int)==1 % vertical interface
                                if p.posn_int{int}(2,m_C)<p.posn_int{int}(2,n) % slave-node n is above master-node m_C
                                    m_R=m_C+1;
                                else
                                    m_R=m_C-1;
                                end
                                if m_R>N_int-n_convert*(N_int-N_s_self(int)) % Slave-Node is not between outermost master nodes
                                    if p.addNTSLMs==1
                                        % Slave-node is on the right side of the rightmost master-node
                                        if p.posn_int{int}(2,n)>p.posn_int{int}(2,m_C) && p.posn_int{int}(2,n)<p.posn_int{int}(2,m_C)+p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(2,n-1)>p.posn_int{int}(2,m_C)-p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(2,n-1)<p.posn_int{int}(2,m_C) 
                                            disp('vertikal')
                                            disp(['m_R ' num2str(m_R)])
                                            m_R=m_C-1;
                                            disp(m_R)
                                            ratio=ratio/abs(p.posn_int{int}(2,m_C)-p.posn_int{int}(2,m_R));
                                        else
                                            m_R=0;
                                        end
                                    else
                                        m_R=0;
                                    end
                                else
                                    if m_R<1+(1-n_convert)*N_s_self(int)
                                        if p.addNTSLMs==1
                                            if p.posn_int{int}(2,n)<p.posn_int{int}(2,m_C) && p.posn_int{int}(2,n)>p.posn_int{int}(2,m_C)-p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(2,n+1)>p.posn_int{int}(2,m_C) && p.posn_int{int}(2,n+1)<p.posn_int{int}(2,m_C)+p.elHeight(p.posn_int{int}(4,m_C))
                                                disp('vertikal')
                                                disp(['m_R ' num2str(m_R)])
                                                m_R=m_C+1;
                                                disp(m_R)
                                                ratio=-ratio/abs(p.posn_int{int}(2,m_C)-p.posn_int{int}(2,m_R));
                                            else
                                                m_R=0;
                                            end
                                        else
                                            m_R=0;
                                        end
                                    else
                                        ratio=ratio/abs(p.posn_int{int}(2,m_C)-p.posn_int{int}(2,m_R));
                                    end
                                end
                            else % horizontal interface
                                if p.posn_int{int}(1,m_C)<p.posn_int{int}(1,n) % slave-node n on the right hand of master-node m_C
                                    m_R=m_C+1;
                                else
                                    m_R=m_C-1;
                                end
                                if m_R>N_int-n_convert*(N_int-N_s_self(int)) % Slave-Node is not between outermost master nodes
                                    if p.addNTSLMs==1
                                        % Slave-node is on the right side of the rightmost master-node
                                        if p.posn_int{int}(1,n)>p.posn_int{int}(1,m_C) && p.posn_int{int}(1,n)<p.posn_int{int}(1,m_C)+p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(1,n-1)>p.posn_int{int}(1,m_C)-p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(1,n-1)<p.posn_int{int}(1,m_C) 
                                            disp('horizontal')
                                            disp(['m_R ' num2str(m_R)])
                                            m_R=m_C-1;
                                            disp(m_R)
                                            ratio=ratio/abs(p.posn_int{int}(1,m_C)-p.posn_int{int}(1,m_R));
                                        else
                                            m_R=0;
                                        end
                                    else
                                        m_R=0;
                                    end
                                else
                                    if m_R<1+(1-n_convert)*N_s_self(int)
                                        if p.addNTSLMs==1
                                            if p.posn_int{int}(1,n)<p.posn_int{int}(1,m_C) && p.posn_int{int}(1,n)>p.posn_int{int}(1,m_C)-p.elHeight(p.posn_int{int}(4,m_C)) && p.posn_int{int}(1,n+1)>p.posn_int{int}(1,m_C) && p.posn_int{int}(1,n+1)<p.posn_int{int}(1,m_C)+p.elHeight(p.posn_int{int}(4,m_C))
                                                disp('horizontal')
                                                disp(['m_R ' num2str(m_R)])
                                                m_R=m_C+1;
                                                disp(m_R)
                                                ratio=-ratio/abs(p.posn_int{int}(1,m_C)-p.posn_int{int}(1,m_R));
                                            else
                                                m_R=0;
                                            end
                                        else
                                            m_R=0;
                                        end
                                    else
                                        ratio=ratio/abs(p.posn_int{int}(1,m_C)-p.posn_int{int}(1,m_R));
                                    end
                                end
                            end
                            
                            if m_R>0
                                disp(['lm ' num2str(lm) ' und '  num2str(lm+1)])
                                p.Bs{p.posn_int{int}(4,n)}(lm+1,p.posn_int{int}(3,n)*2)=-1;
                                p.Bs{p.posn_int{int}(4,n)}(lm,p.posn_int{int}(3,n)*2-1)=-1;
                                p.Bs{p.posn_int{int}(4,m_C)}(lm+1,p.posn_int{int}(3,m_C)*2)=1-ratio;
                                p.Bs{p.posn_int{int}(4,m_C)}(lm,p.posn_int{int}(3,m_C)*2-1)=1-ratio;
                                p.Bs{p.posn_int{int}(4,m_R)}(lm+1,p.posn_int{int}(3,m_R)*2)=ratio;
                                p.Bs{p.posn_int{int}(4,m_R)}(lm,p.posn_int{int}(3,m_R)*2-1)=ratio;
                                lm=lm+2;
                                p.interfaceDOF{p.posn_int{int}(4,n)}(p.posn_int{int}(3,n)*2,1)=1;
                                p.interfaceDOF{p.posn_int{int}(4,n)}(p.posn_int{int}(3,n)*2-1,1)=1;
                                p.interfaceDOF{p.posn_int{int}(4,m_C)}(p.posn_int{int}(3,m_C)*2,1)=1;
                                p.interfaceDOF{p.posn_int{int}(4,m_C)}(p.posn_int{int}(3,m_C)*2-1,1)=1;  
                                p.interfaceDOF{p.posn_int{int}(4,m_R)}(p.posn_int{int}(3,m_R)*2,1)=1;
                                p.interfaceDOF{p.posn_int{int}(4,m_R)}(p.posn_int{int}(3,m_R)*2-1,1)=1; 

                                p.Nn_s_ncFull=p.Nn_s_ncFull-1;
                            end
                        else
                            disp(['conf lm ' num2str(lm) ' und '  num2str(lm+1)])
                            p.Bs{p.posn_int{int}(4,n)}(lm+1,p.posn_int{int}(3,n)*2)=1;
                            p.Bs{p.posn_int{int}(4,n)}(lm,p.posn_int{int}(3,n)*2-1)=1;
                            p.Bs{p.posn_int{int}(4,n_conf)}(lm+1,p.posn_int{int}(3,n_conf)*2)=-1;
                            p.Bs{p.posn_int{int}(4,n_conf)}(lm,p.posn_int{int}(3,n_conf)*2-1)=-1;
                            lm=lm+2;
                            p.interfaceDOF{p.posn_int{int}(4,n)}(p.posn_int{int}(3,n)*2,1)=1;
                            p.interfaceDOF{p.posn_int{int}(4,n)}(p.posn_int{int}(3,n)*2-1,1)=1;
                            p.interfaceDOF{p.posn_int{int}(4,n_conf)}(p.posn_int{int}(3,n_conf)*2,1)=1;
                            p.interfaceDOF{p.posn_int{int}(4,n_conf)}(p.posn_int{int}(3,n_conf)*2-1,1)=1;

                            p.Nn_s_ncFull=p.Nn_s_ncFull-1;
                        end
                    end
                    
                case 'Mortar'
                    % build up set of slave nodes and projected
                    % master-nodes
                    s=p.int(1+n_convert,int);   %Slave substructure of interface
                    z=p.int(2-n_convert,int);   %Master substructure of interface
                    
                    switch p.int(3,int)
                        case 1 % vertical interface
                            int_dir=2;
                            active_set=zeros(2,p.Nely(s)+1+modfl(p.Nely(s),p.Nely(z),tol));%size(p.posn_int{int},2)); % first row: slave side nodes, second row: master-side nodes
                        case 2 % horizontal interface
                            int_dir=1;
                            active_set=zeros(2,p.Nelx(s)+1+modfl(p.Nelx(s),p.Nelx(z),tol));%size(p.posn_int{int},2)); % first row: slave side nodes, second row: master-side nodes
                    end
                    
                    i=1;
                    disp(['int_dir: ' num2str(int_dir)])
                    disp('size active_set:')
                    disp(size(active_set))
                    disp('Slave Nodes:')
                    disp(1+n_convert*N_s_self(int):N_s_self(int)+n_convert*(N_int-N_s_self(int)))
                    disp(p.posn_int{int}(3,1+n_convert*N_s_self(int):N_s_self(int)+n_convert*(N_int-N_s_self(int))))
                    disp('Master Nodes:')
                    disp(1+(1-n_convert)*N_s_self(int):N_int-n_convert*(N_int-N_s_self(int)))
                    for n=1+n_convert*N_s_self(int):N_s_self(int)+n_convert*(N_int-N_s_self(int)) % jump through slave-nodes
                        %disp(['n: ' num2str(n)])
                        j=1;
                        p.connection=[];
                        for m=1+(1-n_convert)*N_s_self(int):N_int-n_convert*(N_int-N_s_self(int)) % jump through master-nodes
                            if int==4 && m_C==7
                                disp(['m: ' num2str(m)])
                                disp(['p.posn_int{int}(int_dir,m)-p.posn_int{int}(int_dir,n): ' num2str(p.posn_int{int}(int_dir,m)-p.posn_int{int}(int_dir,n))])
                            end
                            %if p.posn_int{int}(int_dir,m)-p.posn_int{int}(int_dir,n)>-tol % only save difference to master nodes on top or equal to current slave node
                            p.connection(1,j)=m;
                            p.connection(2,j)=p.posn_int{int}(int_dir,m)-p.posn_int{int}(int_dir,n);
                            j=j+1;
                            %end
                        end
                        [dist,m_C]=min(abs(p.connection(2,:)));
                        m_C=int32(p.connection(1,m_C));
                        disp(['n: ' num2str(n)])
                        disp(['m_C: ' num2str(m_C)])
                        disp(['dist: ' num2str(dist)])
                        if abs(dist)<=tol
                            active_set(1,i)=n;
                            active_set(2,i)=m_C;
                            i=i+1;
                        else
                            %k= p.connection(1,:)==m_C;                               
                            
                            disp('p.connection: ')
                            disp(p.connection)
                            k=1;
                            while k<=size(p.connection,2) % Search for closest master-node
                                %disp(['i<size(active_set,2): ' num2str(i<size(active_set,2))])
                                %disp(['p.connection(2,k)-p.elHeight(s)<-tol: ' num2str(p.connection(2,k)-p.elHeight(s)<-tol) ])                                    
                                if p.connection(2,k)<0 && i==1
                                    if p.connection(2,k+1)>0
                                        active_set(1,i)=int32(p.connection(1,k));
                                        active_set(2,i)=int32(p.connection(1,k));
                                        i=i+1;
                                    end
                                else
                                    if abs(p.connection(2,k))-p.elHeight(s)<-tol%&& n<N_s_self(int)+n_convert*(N_int-N_s_self(int))
                                        active_set(1,i)=n;
                                        active_set(2,i)=n;
                                        i=i+1;
                                        l=k;
                                        while l<=size(p.connection,2)
                                            if abs(p.connection(2,l))-p.elHeight(s)<-tol && p.connection(2,l)>0
                                                active_set(1,i)=int32(p.connection(1,l));
                                                active_set(2,i)=int32(p.connection(1,l));
                                                i=i+1;
                                            end
                                            l=l+1;
                                        end
                                        break                                            
                                    end
                                    if k==size(p.connection,2) && p.connection(2,k)>0 && i>1
                                        active_set(1,i)=n;
                                        active_set(2,i)=n;
                                        i=i+1;
                                    end
                                end
                                k=k+1;
                            end
                        end
                    end
                    i=1;
                    while i<=size(active_set,2)
                        if active_set(1,i)==0
                            active_set(:,i)=[];
                            i=1;
                        else
                            i=i+1;
                        end
                    end
                    % Check, if slave nodes have been slave nodes at
                    % another interface e.g. substructure 3 has a
                    % slave-contact with subs 1 and subs 2 and subs1 and
                    % subs2 are contacting each other as well
                    slaveint=find(s==slave_substructures(1,:));
                    disp('p.int:')
                    disp(p.int)
                    disp(['int: ' num2str(int) ', slaveint: ' num2str(slaveint) ', s: ' num2str(s)])
                    disp(slave_substructures)
                    disp(size(slaveint))
                    contact_case=0;
                    int_counter=1;
                    if isempty(slaveint)==0
                        for k=1:size(slaveint,2)
                            if p.int(3,int)==p.int(3,slaveint(k))
                                if s==p.int(1,slaveint(k)) % note the substructure in subscontact_old, which was in contact with s at interface slaveint
                                    subscontact_old=p.int(2,slaveint(k));
                                else
                                    subscontact_old=p.int(1,slaveint(k));
                                end
                                int_counter=int-1;
                                disp(['subscontact_old: ' num2str(subscontact_old) ', z: ' num2str(z)])
                                while int_counter>0
                                    disp(['int_counter: ' num2str(int_counter)])
                                    %disp(p.int)
                                    % Check, if there was already a contact between
                                    % the two master-substructures for s
                                    if (p.int(1,int_counter)==z && p.int(2,int_counter)==subscontact_old) || (p.int(2,int_counter)==z && p.int(1,int_counter)==subscontact_old)

                                        if active_set(1,1)==active_set(2,1)
                                            contact_case=1; % Last node of interface was nonconforming => all four last LMs needed
                                        else
                                            contact_case=2; % Last node of interface was conforming => only two last LMs needed
                                        end
                                        break

                                    end
                                    int_counter=int_counter-1;
                                end
                            end
                        end
                    end
                    disp(['contact_case: ' num2str(contact_case) ', int_counter: ' num2str(int_counter)])
                    lag_mult=zeros(1,(N_s_self(int)+n_convert*(N_int-N_s_self(int))-(1+n_convert*N_s_self(int))+1)*2);
                    disp('lm_backup: ')
                    disp(lm_backup)
                    i_add=0;
                    for i=1:size(lag_mult,2)
                        if i<=4 && contact_case==1
                            lag_mult(i)=lm_backup(slaveint,i);%int_counter
                            i_add=i_add+1;
                        else
                            if i<=2 && contact_case==2
                                lag_mult(i)=lm_backup(slaveint,i+2);
                                i_add=i_add+1;
                            else
                                lag_mult(i)=lm+i-i_add-1;
                            end
                        end
                    end
                    
                    disp('lag_mult:')
                    disp(lag_mult)
                    %{
                    disp(['p.elHeight: ' num2str(p.elHeight(s))])
                    %}
                    disp(['Active Set (Interface ' num2str(int) '):'])
                    disp(active_set)
                    disp(size(active_set))
                    
                    % calculate mortar-integrals
                    subs=[s,z];
                    %disp(['subs(2): ' num2str(z)])
                    lm_i=1; % index for lagrange multipliers in lag_mult
                    interface_boundaries=0;
                    n_expect=active_set(1,1);
                    N_counter=0;
                    n_ma=0;
                    for n=1:size(active_set,2)
                        if active_set(1,n)==n_expect
                            N_counter=N_counter+1;
                            n_expect=n_expect+1;
                            if active_set(1,1)~=active_set(2,1)
                                n_ma=active_set(2,1);
                            end
                        else
                            if n_ma==0
                                n_ma=active_set(2,n);
                            end
                        end
                    end
                    if N_counter>size(active_set,2)/2
                        n_sl=active_set(1,1);
                        if n_ma==0
                            n_ma=active_set(2,1);
                        end
                    else
                        n_sl=n_ma;
                        n_ma=active_set(1,1);
                    end
                    
                    
                    disp(['Lower boundary: ' num2str(1+n_convert*N_s_self(int))])
                    activate_sl=0;
                    activate_ma=0;
                    for n=1:size(active_set,2)
                        if active_set(1,n)==n_sl
                            activate_sl=1;
                        end
                        if active_set(2,n)==n_ma
                            activate_ma=1;
                        end
                        if activate_sl==1 && activate_ma==1%active_set(1,n)==1+n_convert*N_s_self(int) || active_set(2,n)==1+(1-n_convert)*N_s_self(int)
                            interface_boundaries=1;
                            n_start=n;
                            activate_sl=0;
                            activate_ma=0;
                        %{
                        else
                            if n==size(active_set,2) || active_set(2,n)==N_int-n_convert*(N_int-N_s_self(int)) || active_set(1,n)==N_s_self(int)+n_convert*(N_int-N_s_self(int))
                                interface_boundaries=0;
                            end
                            %}
                        end
                        if n==size(active_set,2) || active_set(2,n)==N_int-n_convert*(N_int-N_s_self(int)) || active_set(1,n)==N_s_self(int)+n_convert*(N_int-N_s_self(int))
                            interface_boundaries=0;
                        end
                        if interface_boundaries==1
                            disp(['n: ' num2str(n)])
                            %disp(['n: ' num2str(n)])
                            % calculate element-coordinates of the
                            % mortar-segment's nodes
                            xi_a=zeros(1,2);                            xi_b=xi_a;
                            for i=1:2
                                disp(['Subs(' num2str(i) '): ' num2str(subs(i))])
                                disp(['elHeight(' num2str(subs(i)) '): ' num2str(p.elHeight(subs(i)))])
                                disp(['p.posn_int{int=' num2str(int) '}(int_dir=' num2str(int_dir) ',active_set(i,n)=' num2str(active_set(i,n)) ')-p.poss(int_dir=' num2str(int_dir) ',subs(i)=' num2str(subs(i)) '): ' num2str(p.posn_int{int}(int_dir,active_set(i,n))-p.poss(int_dir,subs(i)))])
                                disp(['Modulo(' num2str(i) '): ' num2str(modfl(p.posn_int{int}(int_dir,active_set(i,n))-p.poss(int_dir,subs(i)),p.elHeight(subs(i)),tol))])
                                xi_a(i)=-1+2*modfl(p.posn_int{int}(int_dir,active_set(i,n))-p.poss(int_dir,subs(i)),p.elHeight(subs(i)),tol)/p.elHeight(subs(i));
                                xi_b(i)=xi_a(i)+2*(p.posn_int{int}(int_dir,active_set(i,n+1))-p.posn_int{int}(int_dir,active_set(i,n)))/p.elHeight(subs(i));
                            end
                            
                            disp('xi_a: ')
                            disp(xi_a)
                            disp('xi_b: ')
                            disp(xi_b)
                            
                            J_seg=p.elHeight(subs(1))*(xi_b(1)-xi_a(1))/4;
                            %{
                            for l=1:2
                                for k=1:2
                                    M_seg{n}(l,k)=(shapefunc(l,xi_a(1))*shapefunc(k,xi_a(1))+shapefunc(l,xi_b(1))*shapefunc(k,xi_b(1)))*J_seg;
                                    N_seg{n}(l,k)=(shapefunc(l,xi_a(1))*shapefunc(k,xi_a(2))+shapefunc(l,xi_b(1))*shapefunc(k,xi_b(2)))*J_seg;
                                end
                            end
                            %}
                            %{
                            n_sl=n;
                            n_ma=n;
                            while active_set(1+n_convert,n_sl)>=active_set(2-n_convert,1) && n_sl>1 %active_set(1+n_convert,size(active_set,2)) && n_sl>1
                                n_sl=n_sl-1;
                            end
                            while active_set(2-n_convert,n_ma)<active_set(2-n_convert,1) && n_ma>1 %active_set(1+n_convert,size(active_set,2)) && n_ma>1
                                n_ma=n_ma-1;
                            end
                            disp(['n_convert: ' num2str(n_convert)])
                            
                            if n_convert==1
                                n_sl_backup=n_sl; % Switch n_sl and n_ma, because in p.posn_int nodes are numbered from 1 to last interface-node beginning with nodes of current substructure independent of master and slave sides
                                n_sl=n_ma;
                                n_ma=n_sl_backup;
                            end
                            %}
                            
                            if n>n_start && (xi_a(1)+1)<=tol%n==n_sl && n>n_start
                                lm_i=lm_i+2;
                                n_sl=n_sl+1;
                            end                                
                            if n>n_start && (xi_a(2)+1)<=tol%n==n_sl && n>n_start
                                n_ma=n_ma+1;
                            end  
                            disp(['n_sl: ' num2str(n_sl)])
                            disp(['n_ma: ' num2str(n_ma)])
                            disp(['xi_a(1)+1: ' num2str(xi_a(1)+1)])
                            disp(['lm_i: ' num2str(lm_i)])
                            
                            % Gauss-Points:
                            eta1=-sqrt(1/3);
                            eta2=sqrt(1/3);  
                            
                            for l=1:2
                                %{
                                disp(['l: ' num2str(l)])
                                disp(['lm_i: ' num2str(lm_i)])
                                disp(['lag_mult(lm_i+(l-1)*2): ' num2str(lag_mult(lm_i+(l-1)*2))])
                                %}
                                for k=1:2
                                    %disp(['k: ' num2str(k)])
                                    n_slave=n_sl+k-1;%active_set(1,n_sl)+k-1;
                                    n_master=n_ma+k-1;%active_set(2,n_ma)+k-1;
                                    
                                    disp(['n_slave: ' num2str(n_slave)])
                                    disp(['n_master: ' num2str(n_master)])
                                    
                                    disp(['Bs(slave)-Eintrag: ' num2str((shapefunc(l,0,xi_a(1))*shapefunc(k,0,xi_a(1))+shapefunc(l,0,xi_b(1))*shapefunc(k,0,xi_b(1)))*J_seg)])
                                    disp(['Bs(master)-Eintrag: ' num2str(-(shapefunc(l,0,xi_a(1))*shapefunc(k,0,xi_a(2))+shapefunc(l,0,xi_b(1))*shapefunc(k,0,xi_b(2)))*J_seg)])
                                                                   
                                    subs1=p.posn_int{int}(4,n_slave);
                                    subs2=p.posn_int{int}(4,n_master);
                                    disp(['subs1: ' num2str(subs1)])
                                    disp(['subs2: ' num2str(subs2)])
                                    disp(['lm1: ' num2str(lag_mult(lm_i+(l-1)*2))])
                                    disp(['lm2: ' num2str(lag_mult(lm_i+1+(l-1)*2))])
                                    
                                    %p.Bs{subs1}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_slave)*2-1)=p.Bs{subs1}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_slave)*2-1)+(shapefunc(l,0,xi_a(1))*shapefunc(k,0,xi_a(1))+shapefunc(l,0,xi_b(1))*shapefunc(k,0,xi_b(1)))*J_seg;
                                    %p.Bs{subs2}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_master)*2-1)=p.Bs{subs2}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_master)*2-1)-(shapefunc(l,0,xi_a(1))*shapefunc(k,0,xi_a(2))+shapefunc(l,0,xi_b(1))*shapefunc(k,0,xi_b(2)))*J_seg;
                                    xig1_1=0.5*(1-eta1)*xi_a(1)+0.5*(1+eta1)*xi_b(1);
                                    xig1_2=0.5*(1-eta2)*xi_a(1)+0.5*(1+eta2)*xi_b(1);
                                    p.Bs{subs1}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_slave)*2-1)=p.Bs{subs1}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_slave)*2-1)+(shapefunc(l,0,xig1_1)*shapefunc(k,0,xig1_1)+shapefunc(l,0,xig1_2)*shapefunc(k,0,xig1_2))*J_seg;
                                    xig2_1=0.5*(1-eta1)*xi_a(2)+0.5*(1+eta1)*xi_b(2);
                                    xig2_2=0.5*(1-eta2)*xi_a(2)+0.5*(1+eta2)*xi_b(2);
                                    p.Bs{subs2}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_master)*2-1)=p.Bs{subs2}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_master)*2-1)-(shapefunc(l,0,xig1_1)*shapefunc(k,0,xig2_1)+shapefunc(l,0,xig1_2)*shapefunc(k,0,xig2_2))*J_seg;
                                    p.Bs{subs1}(lag_mult(lm_i+1+(l-1)*2),p.posn_int{int}(3,n_slave)*2)=p.Bs{subs1}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_slave)*2-1);
                                    p.Bs{subs2}(lag_mult(lm_i+1+(l-1)*2),p.posn_int{int}(3,n_master)*2)=p.Bs{subs2}(lag_mult(lm_i+(l-1)*2),p.posn_int{int}(3,n_master)*2-1);
                                    
                                    p.interfaceDOF{subs1}(p.posn_int{int}(3,n_slave)*2-1,1)=1;
                                    p.interfaceDOF{subs1}(p.posn_int{int}(3,n_slave)*2,1)=1;
                                    p.interfaceDOF{subs2}(p.posn_int{int}(3,n_master)*2-1,1)=1;
                                    p.interfaceDOF{subs2}(p.posn_int{int}(3,n_master)*2,1)=1;
                                end
                                %lm=lm+2;
                            end
                            %{
                            if n<size(active_set,2)-1
                                lm=lm-2;
                            end
                            %}
                            disp('--------')
                            %{
                            if active_set(1,n+1)==active_set(2,n+1)
                                lm=lag_mult(lm_i);
                            else
                                lm=lag_mult(lm_i+1+(l-1)*2)+1;
                            end
                            %}
                        end
                    end
                    slave_substructures(1,int)=subs1;
                    %{
                    if int<int_full
                        if p.int(4,int)==p.int(4,int+1) && subs1==p.int(1+n_convert,int+1) && p.int(3,int)==p.int(3,int+1)
                            lm=lag_mult(lm_i);
                        else
                            lm=lag_mult(lm_i+1+(l-1)*2)+1;
                        end
                    else
                    %}
                    lm=lag_mult(lm_i+1+(l-1)*2)+1;
                    %end
                    for i=1:4
                        lm_backup(int,5-i)=lm-i;
                    end
                    %disp(['p.Nn_s_ncFull vorher: ' num2str(p.Nn_s_ncFull)])
                    p.Nn_s_ncFull=p.Nn_s_ncFull-int32(size(lag_mult,2)/2);
                    %lm=lag_mult(lm_i+1+(l-1)*2)+1;
                    disp(['lm: ' num2str(lm)])
                    %{
                    for n=1:size(active_set,2)-1
                        disp(['M_seg{' num2str(n) '}: '])
                        disp(M_seg{n}) 
                        disp(['N_seg{' num2str(n) '}: '])
                        disp(-N_seg{n})
                    end
                    %}
            end
        end
    end
    %{
    disp('Interface dof (1):')
    for i=1:size(p.interfaceDOF{1},1)
        if p.interfaceDOF{1}(i,1)==1
            disp(i)
        end
    end

    disp(['p.Nn_s_ncFull: ' num2str(p.Nn_s_ncFull)])
    disp(['lm: ' num2str(lm)])
    %}
    
    while lm-1<size(p.Bs{1},1) % lm-1, because last lm must be a conforming one
        for s=1:p.Ns
            p.Bs{s}(p.Nlm,:)=[];
        end
        p.Nlm=p.Nlm-1;
    end
    for s=1:p.Ns
        disp(['Bs(' num2str(s) '):'])
        disp(p.Bs{s})
        disp(size(p.Bs{s}))
    end
    disp(['p.Nlm: ' num2str(p.Nlm)])
    p.NdofFull = p.Nn_s_ncFull*Ndof_n;
else
    Nnodepairings = ipairings;
    lm = 0;
    p.W = zeros(p.Nlm);
    %p.CornerLMs = zeros(p.Nlm,1);
    for ipairings = 1:Nnodepairings
        iNs_self = nodepairings(ipairings,1);
        iNs_other = nodepairings(ipairings,2);
        node_self = nodepairings(ipairings,3);
        node_other = nodepairings(ipairings,4);
        valence = valenceNodes{iNs_self}(node_self);
        assert(valence == valenceNodes{iNs_other}(node_other), 'valence inconsistent!');
        xdof_self = node_self*2-1;
        xdof_other = node_other*2-1;
        ydof_self = node_self*2;
        ydof_other = node_other*2;

    %     stiffnessSelf_x = Ks{iNs_self}(xdof_self,xdof_self);
    %     stiffnessOther_x = Ks{iNs_other}(xdof_other,xdof_other);
    %     betaScaleDivisor_x = stiffnessSelf_x;
    %     
    %     stiffnessSelf_y = Ks{iNs_self}(ydof_self,ydof_self);
    %     stiffnessOther_y = Ks{iNs_other}(ydof_other,ydof_other);
    %     betaScaleDivisor_y = stiffnessSelf_y;
    %     
    %     for i = 1:size(betaScalePre{iNs_self}{node_self},1)
    %         otherSub = betaScalePre{iNs_self}{node_self}(i,1);
    %         otherNode = betaScalePre{iNs_self}{node_self}(i,2);
    %         betaScaleDivisor_x = betaScaleDivisor_x + Ks{otherSub}(otherNode*2-1,otherNode*2-1);
    %         betaScaleDivisor_y = betaScaleDivisor_y + Ks{otherSub}(otherNode*2,otherNode*2);
    %     end

        % setup a LM for the xdof
        lm = lm + 1;
        p.W(lm,lm) = valence;
        %betaScales{iNs_self}(lm,1) = stiffnessOther_x/betaScaleDivisor_x;
        %betaScales{iNs_other}(lm,1) = stiffnessSelf_x/betaScaleDivisor_x;
        p.Bs{iNs_self}(lm,xdof_self) = 1;
        p.Bs{iNs_other}(lm,xdof_other) = -1;

        % setup a LM for the ydof
        lm = lm + 1;
        p.W(lm,lm) = valence;
        %betaScales{iNs_self}(lm,1) = stiffnessOther_y/betaScaleDivisor_y;
        %betaScales{iNs_other}(lm,1) = stiffnessSelf_y/betaScaleDivisor_y;
        p.Bs{iNs_self}(lm,ydof_self) = 1;
        p.Bs{iNs_other}(lm,ydof_other) = -1;

        if ~isempty(find(CornerPairings == ipairings,1))
            p.CornerLMs(lm-1) = 1;
            p.CornerLMs(lm) = 1;
        end

        % store which dofs are interface dofs
        p.interfaceDOF{iNs_self}(xdof_self,1) = 1;
        p.interfaceDOF{iNs_self}(ydof_self,1) = 1;
        p.interfaceDOF{iNs_other}(xdof_other,1) = 1;
        p.interfaceDOF{iNs_other}(ydof_other,1) = 1;

        % store the LM index to every dof
        interfaceDOFlm{iNs_self}(xdof_self,1) = lm-1;
        interfaceDOFlm{iNs_self}(ydof_self,1) = lm;
        interfaceDOFlm{iNs_other}(xdof_other,1) = lm-1;
        interfaceDOFlm{iNs_other}(ydof_other,1) = lm;


    end
    assert(p.Nlm==lm,'Number of Lagrange Multipliers incorrect!');

end

%% set up DBCs (substructure and full)
for s = 1:p.Ns
    p.zerodofs{s} = [];
end

%% Te
% detect DOF subject to DBCs

%North
if (p.clamping == 1) % here the TOP side of the beam is clamped
    % feti:
    % select all uppermost nodes of all uppermost substructures
    if p.nonconforming==1
        s_end=p.Ns;
        s_step=1;
    else
        s_end=((p.Nsy-1)*p.Nsx+1);
        s_step=p.Nsx;
    end
    i=1;
    clamped_subs=[];
    if p.nonconforming==1
        for s = 1:s_step:s_end
            if p.poss(2,s)+p.sizes(2,s)-p.Height<=tol
                clamped_subs(i)=s;
                i=i+1;
                zeronodes{s}(1:p.Nely(s)+1) = 0;
                p.zerodofs{s}(1:2*(p.Nely(s)+1)) = 0;
                for j = 1:p.Nelx(s)+1
                    zeronodes{s}(j) = (p.Nelx(s)+1)*(p.Nely(s))+j;
                    p.zerodofs{s}(2*j-1) = zeronodes{s}(j)*2-1;
                    p.zerodofs{s}(2*j) = zeronodes{s}(j)*2;
                end
            end
        end
    else
        p.Nscorner(4)=p.Nsx*(p.Nsy-1)+1;
        p.Nscorner(3)=p.Nsy*p.Nsx;
        N_zeronodes=0;
        for s = p.Nscorner(4):1:p.Nscorner(3) % upermost substructures
            zeronodes{s}(1:p.Nelx(s)+1) = 0; % initialisieren(wie viele gibts)
            p.zerodofs{s}(1:2*p.Nelx(s)+1) = 0;
            for i = 1:p.Nelx(s)+1
                zeronodes{s}(i) = (p.Nelx(s)+1)*(p.Nely(s))+i; % uppermost nodes
                %disp(['Zeronode ' num2str(i) ': ' num2str(zeronodes{s}(i))])
                p.zerodofs{s}(2*i-1) = zeronodes{s}(i)*2-1; % entsprechendene DoFs
                p.zerodofs{s}(2*i) = zeronodes{s}(i)*2;
            end
            N_zeronodes=N_zeronodes+p.Nelx(s);
        end    
        % FULL
        % select all uppermost nodes
        s=1;
        zeronodesFull(1:p.Nelx(s)*p.Nsx+1) = 0; % initialisieren
        p.zerodofsFull(1:2*(p.Nelx(s)*p.Nsx+1)) = 0;
        for i = 1:p.Nelx(s)*p.Nsx+1
            zeronodesFull(i) = (p.Nelx(s)*p.Nsx+1)*(p.Nely(s)*p.Nsy) +i;
            p.zerodofsFull(2*i-1) = zeronodesFull(i)*2-1;
            p.zerodofsFull(2*i) = zeronodesFull(i)*2;
        end
    end

% East  
elseif (p.clamping == 2) % here the RIGHT side of the beam is clamped
    % feti:
    % select all rightmost nodes of all rightmost substructures
    for s = p.Nsx:p.Nsx:p.Nsy*p.Nsx % rightmost substructures
        zeronodes{s}(1:p.Nely+1) = 0; % initialisieren(wie viele gibts)
        p.zerodofs{s}(1:2*p.Nely+1) = 0;
        for i = 1:p.Nely+1
            zeronodes{s}(i) = i*(p.Nelx+1); % rightmost nodes
            p.zerodofs{s}(2*i-1) = zeronodes{s}(i)*2-1; % entsprechendene DoFs
            p.zerodofs{s}(2*i) = zeronodes{s}(i)*2;
        end
    end
    % FULL
    % select all rightmost nodes
    zeronodesFull(1:p.Nely*p.Nsy+1) = 0; % initialisieren
    p.zerodofsFull(1:2*(p.Nely*p.Nsy+1)) = 0;
    for i = 1:p.Nely*p.Nsy+1
        zeronodesFull(i) = i*(p.Nelx*p.Nsx);
        p.zerodofsFull(2*i-1) = zeronodesFull(i)*2-1;
        p.zerodofsFull(2*i) = zeronodesFull(i)*2;
    end

% South
elseif (p.clamping == 3) % here the BOTTOM side of the beam is clamped
    % feti:
    % select all bottommost nodes of all bottommost substructures
    for s = 1:p.Nsx % bottommost substructures
        zeronodes{s}(1:p.Nelx+1) = 0; % initialisieren(wie viele gibts)
        p.zerodofs{s}(1:2*p.Nelx+1) = 0;
        for i = 1:p.Nelx+1
            zeronodes{s}(i) = i; % bottommost nodes
            p.zerodofs{s}(2*i-1) = zeronodes{s}(i)*2-1; % entsprechendene DoFs
            p.zerodofs{s}(2*i) = zeronodes{s}(i)*2;
        end
    end
    % FULL
    % select all bottommost nodes
    zeronodesFull(1:p.Nelx*p.Nsx+1) = 0; % initialisieren
    p.zerodofsFull(1:2*(p.Nelx*p.Nsx+1)) = 0;
    for i = 1:p.Nelx*p.Nsx+1
        zeronodesFull(i) = i;
        p.zerodofsFull(2*i-1) = zeronodesFull(i)*2-1;
        p.zerodofsFull(2*i) = zeronodesFull(i)*2;
    end 

% West
elseif (p.clamping == 4) % here the LEFT side of the beam is clamped
    % feti:
    % select all leftmost nodes of all leftmost substructures
    if p.nonconforming==1
        s_end=p.Ns;
        s_step=1;
    else
        s_end=((p.Nsy-1)*p.Nsx+1);
        s_step=p.Nsx;
    end
    i=1;
    clamped_subs=[];
    for s = 1:s_step:s_end
        if p.nonconforming==1
            if p.poss(1,s)<=tol
                clamped_subs(i)=s;
                i=i+1;
                zeronodes{s}(1:p.Nely(s)+1) = 0;
                p.zerodofs{s}(1:2*(p.Nely(s)+1)) = 0;
                for j = 1:p.Nely(s)+1
                    zeronodes{s}(j) = 1+(j-1)*(p.Nelx(s)+1);
                    p.zerodofs{s}(2*j-1) = zeronodes{s}(j)*2-1;
                    if p.novertDBC~=1
                        p.zerodofs{s}(2*j) = zeronodes{s}(j)*2;
                    end
                end
            end
        else
            disp(['s: ' num2str(s)])
            disp(['p.Nsx: ' num2str(p.Nsx)])
            %if mod(s,p.Nsx)==0
            clamped_subs(i)=s;
            i=i+1;
            zeronodes{s}(1:p.Nely(s)+1) = 0;
            p.zerodofs{s}(1:2*p.Nely(s)+1) = 0;
            for j = 1:p.Nely(s)+1
                zeronodes{s}(j) = 1+(j-1)*(p.Nelx(s)+1);
                p.zerodofs{s}(2*j-1) = zeronodes{s}(j)*2-1;
                if p.novertDBC~=1
                    p.zerodofs{s}(2*j) = zeronodes{s}(j)*2;
                end
            end
            %end
        end
        
        j=1;
        while j<=size(p.zerodofs{s},2)
            if p.zerodofs{s}(j)==0
                p.zerodofs{s}(j)=[];
            else
                j=j+1;
            end
        end
    end
    % FULL
    % select all leftmost nodes
    if p.nonconforming~=1
        s=p.Nsy;
        zeronodesFull(1:p.Nely(s)*p.Nsy+1) = 0;
        p.zerodofsFull(1:2*(p.Nely(s)*p.Nsy+1)) = 0;
        for i = 1:p.Nely(s)*p.Nsy+1
            zeronodesFull(i) = 1+(i-1)*(p.Nelx(s)*p.Nsx+1);
            p.zerodofsFull(2*i-1) = zeronodesFull(i)*2-1;
            p.zerodofsFull(2*i) = zeronodesFull(i)*2;
        end
    end
else 
    display(['clamping case undefined']);
    return;
end

% feti:
% find LMs which correspond to DBC-dof
% those LMs (the corresponding lines) have to be deleted
% from B respectively all p.Bs{}
p.Ndof_s_Free = p.Ndof_s;

if p.nonconforming~=1
    foundlms(1:p.Nlm) = 0;
    p.Llms = eye(p.Nlm);
end
p.zerolms = [];
i_zerolms = 0;

%disp('Clamped Substructures:')
%disp(clamped_subs)
for s = 1:p.Ns
    if isempty(p.zerodofs{s}) %(length(p.zerodofs{s})<s || isempty(p.zerodofs{s}))
        continue;
    end
    j=1;
    while j<=length(clamped_subs)
        if clamped_subs(j)==s
            for i = 1:length(p.zerodofs{s})
                if ( ...
                        interfaceDOFlm{s}(p.zerodofs{s}(i)) > 0 && ...
                        foundlms(interfaceDOFlm{s}(p.zerodofs{s}(i))) == 0 ...
                   )

                    i_zerolms = i_zerolms + 1;
                    foundlms(interfaceDOFlm{s}(p.zerodofs{s}(i))) = 1;
                    rowtodel = interfaceDOFlm{s}(p.zerodofs{s}(i))-(i_zerolms-1);
                    p.Llms(rowtodel,:) = [];
                    p.zerolms(end+1) = interfaceDOFlm{s}(p.zerodofs{s}(i));
                end
            end

            % correct the number of dof of the substructure with DBC
            p.Ndof_s(s) = p.Ndof_s(s)-length(p.zerodofs{s});
            %disp(['InterfaceDOF Substructure ' num2str(s) ':'])
            %disp(size(p.interfaceDOF{s}))
            p.interfaceDOF{s}(p.zerodofs{s},:) = [];
        end
        j=j+1;
    end
end
% compute number of dof in the whole domain
p.Ndof = sum(p.Ndof_s);

% correct the number of lms with DBC
p.Nlm = p.Nlm - i_zerolms;

p.Nbdof = 0;
p.Nidof = 0;
for s=1:p.Ns
    p.bDOF{s} = [];
    p.iDOF{s} = [];
    for i=1:p.Ndof_s(s)
        if p.interfaceDOF{s}(i) == 1
            p.bDOF{s}(end+1) = i;
            p.Nbdof = p.Nbdof + 1;
        else
            p.iDOF{s}(end+1) = i;
            p.Nidof = p.Nidof + 1;
        end
    end
    
    p.L{s} = eye(p.Ndof_s_Free(s));
    p.L{s}(p.zerodofs{s},:) = [];
end

if p.nonconforming~=1
    p.Lfull = eye(p.NdofFull);
    p.Lfull(p.zerodofsFull,:) = [];

    p.W = p.Llms*p.W*p.Llms';
    p.Winv = inv(p.W);
end
%betaScales{s} = diag(p.Llms*betaScales{s});

end