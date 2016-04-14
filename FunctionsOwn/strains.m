function [p] = strains(X0,d,s,p)

% Function to postprocess the calculated displacements and calculate
% strains
el_jump=0;
for e=1:p.Nelx(s)*p.Nely(s)
    A=zeros(8,size(X0,1)); % Map global dofs of the element's upper nodes to local dofs of the element in opposite clock direction
    
    if mod(e-1,p.Nelx(s))==0 && e>1
        el_jump=el_jump+2;
    end

    for i=1:4
        A(i,el_jump+2*(e-1)+i)=1;
        if i<=2 
            A(4+i,2*(p.Nelx(s)+1)+2*(e-1)+el_jump+i+2)=1;
        else
            A(4+i,2*(p.Nelx(s)+1)+2*(e-1)+el_jump+i-2)=1;
        end
    end
    
    p.x_el{s}{e}=A*X0;
    p.d_el{s}{e}=A*d;
    
    el_x=p.x_el{s}{e}(1:2:size(p.x_el{s}{e},1));
    el_y=p.x_el{s}{e}(2:2:size(p.x_el{s}{e},1));
    
    xi=[-1 1 1 -1;
        -1 -1 1 1];
    
    p.eps{s}{e}=zeros(3,4);
    if p.cal_stress==1
        p.stress{s}{e}=zeros(4,4);
    end
    for node=1:4
        %N1=shapefunc(1,0,xi(1,node))*shapefunc(1,0,xi(2,node));
        %N2=shapefunc(2,0,xi(1,node))*shapefunc(1,0,xi(2,node));
        %N3=shapefunc(2,0,xi(1,node))*shapefunc(2,0,xi(2,node));
        %N4=shapefunc(1,0,xi(1,node))*shapefunc(2,0,xi(2,node));
        
        J=zeros(2);
        
        Nndxi=zeros(2,4);
        
        for i=1:2
            N1dxi=shapefunc(1,mod(i,2),xi(1,node))*shapefunc(1,i-1,xi(2,node));
            N2dxi=shapefunc(2,mod(i,2),xi(1,node))*shapefunc(1,i-1,xi(2,node));
            N3dxi=shapefunc(2,mod(i,2),xi(1,node))*shapefunc(2,i-1,xi(2,node));
            N4dxi=shapefunc(1,mod(i,2),xi(1,node))*shapefunc(2,i-1,xi(2,node));
            
            Nndxi(i,:)=[N1dxi N2dxi N3dxi N4dxi];
            
            Ndxi=[N1dxi 0 N2dxi 0 N3dxi 0 N4dxi 0;
                0 N1dxi 0 N2dxi 0 N3dxi 0 N4dxi];
            
            J(i,:)=(Ndxi*p.x_el{s}{e})'; % Jacobian matrix of element e at node
        end
        J_inv = inv(J);
        
        Ndx=J_inv*Nndxi;
        
        % Differential operator L*shapefunctions N:
        LN=[Ndx(1,1) 0 Ndx(1,2) 0 Ndx(1,3) 0 Ndx(1,4) 0;
            0 Ndx(2,1) 0 Ndx(2,2) 0 Ndx(2,3) 0 Ndx(2,4);
            Ndx(2,1) Ndx(1,1) Ndx(2,2) Ndx(1,2) Ndx(2,3) Ndx(1,3) Ndx(2,4) Ndx(1,4)];
        
        
        p.eps{s}{e}(:,node)=LN*p.d_el{s}{e};
        
        if p.cal_stress==1
            r=floor(e/p.Nelx(s))+1;
            % Constitutive matrix for plane stress
            if (~isempty(find(r == p.SteelRowNrsOddNsx,1))) || (~isempty(find(r == p.SteelRowNrsEvenNsx,1)))
                if p.plain==1
                    C=p.ESt/(1+p.nuSt^2)*[1 p.nuSt 0;
                        p.nuSt 1 0;
                        0 0 (1-p.nuSt)/2];
                elseif p.plain==2
                    C=p.ESt/((1+p.nuSt)*(1-2*p.nuSt))*[1-p.nuSt p.nuSt 0;
                        p.nuSt 1-p.nuSt 0;
                        0 0 (1-2*p.nuSt)/2];
                end
            else
                if p.plain==1
                    C=p.ERub/(1+p.nuRub^2)*[1 p.nuRub 0;
                        p.nuRub 1 0;
                        0 0 (1-p.nuRub)/2];
                elseif p.plain==2
                    C=p.ERub/((1+p.nuRub)*(1-2*p.nuRub))*[1-p.nuRub p.nuRub 0;
                        p.nuRub 1-p.nuRub 0;
                        0 0 (1-2*p.nuRub)/2];
                end
            end
            p.stress{s}{e}(1:3,node)=C*p.eps{s}{e}(:,node);
            p.stress{s}{e}(4,node)=sqrt(p.stress{s}{e}(1,node)^2+p.stress{s}{e}(2,node)^2-p.stress{s}{e}(1,node)*p.stress{s}{e}(2,node)+3*p.stress{s}{e}(3,node)^2);
        end
    end
    if p.plot_int~=0
        int=p.plot_int;
        if p.int(3,int)~=3  
            for n=1:size(p.posn_int{int},2)
                if p.posn_int{int}(4,n)==s
                    var_help=find(el_x==p.posn_int{int}(1,n));
                    if isempty(var_help)==0
                        for i=1:size(var_help,1)
                            if el_y(var_help(i))==p.posn_int{int}(2,n)
                                if p.cal_stress==1
                                    p.eps_int{s}(1,n)=p.stress{s}{e}(p.strain_dir,var_help(i));
                                    if p.int(3,int)==1
                                        p.eps_int{s}(2,n)=p.posn_int{int}(2,n);
                                    else
                                        p.eps_int{s}(2,n)=p.posn_int{int}(1,n);
                                    end
                                else
                                    p.eps_int{s}(1,n)=p.eps{s}{e}(p.strain_dir,var_help(i));
                                    if p.int(3,int)==1
                                        p.eps_int{s}(2,n)=p.posn_int{int}(2,n);
                                    else
                                        p.eps_int{s}(2,n)=p.posn_int{int}(1,n);
                                    end
                                end
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    %disp(['strains(e=' num2str(e) '):'])
    %disp(eps{e})
end

end