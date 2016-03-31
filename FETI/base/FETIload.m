function p = FETIload(p)

switch p.mode
    case 'dynamic'
        % loadcase (set force vector f)
        for s = 1:p.Ns
            p.fs{s}(1:p.Ndof_s_Free(s),1:length(p.t)) = 0;
            if p.nonconforming==1
                p.fsRevert{s}(1:p.Nn_s(s)) = 0;
            else
                p.fsRevert{s}(1:p.Nn_s) = 0;
            end            
        end
        % substructure-local node number of right-upper y-dof: p.Ndof_s_Free (last one)
        for i= p.Loadcase
            if i == 1 %rechts oben y
                p.fs{p.Ns}(p.Ndof_s_Free(p.Ns),:) = p.Loading; 
            elseif i == 2 %rechts unten y
                p.fs{p.Nsx}(2*(p.Nelx+1),:) = -p.Loading; 
            elseif i == 3 %links unten y
                p.fs{1}(2,:) = p.Loading;        
            elseif i == 4 %links oben y
                p.fs{p.Ns-p.Nsx+1}(2*((p.Nelx+1)*p.Nely+1),:) = p.Loading; 
            elseif i == 5 %rechts oben x
                p.fs{p.Ns}(p.Ndof_s_Free(p.Ns)-1,:) = p.Loading; 
            elseif i == 6 %rechts unten x
                p.fs{p.Nsx}(2*(p.Nelx+1)-1,:) = p.Loading; 
            elseif i == 7 %links unten x  
                p.fs{1}(1,:) = -p.Loading;      
            elseif i == 8 %links oben x
                p.fs{p.Ns-p.Nsx+1}(2*((p.Nelx+1)*p.Nely+1)-1,:) = -p.Loading;
            elseif i == 9 % mitte oben(Interfacegrenze) y
                p.fs{p.Ns-round(p.Nsx/2)+1}(p.Ndof_s_Free(p.Ns),:) = p.Loading; 
            elseif i == 10 %mitte oben2 y 
                p.fs{p.Ns-round(p.Nsx/2)+2}(p.Ndof_s_Free(p.Ns)-(round(p.Nelx/2))-2,:) = p.Loading; 
            else
                display(['loadcase undefined']);
                return;
            end
        end
    case 'static'
        %% loadcase (set force vector f)
        for s = 1:p.Ns
            p.fs{s}(1:p.Ndof_s_Free(s),1) = 0;
            if p.nonconforming==1
                p.fsRevert{s}(1:p.Nn_s(s)) = 0;
            else
                p.fsRevert{s}(1:p.Nn_s) = 0;
            end
        end
        % substructure-local node number of right-upper y-dof: Ndof_s (last one)
        switch p.Loadcase 
            case 1  % rechts oben y
                p.fs{p.Ns}(p.Ndof_s_Free(p.Ns)) = -p.bendforce;
            case 2  % rechts oben x
                p.fs{p.Ns}(p.Ndof_s_Free(p.Ns)-1) = p.axforce;
                %p.fs{p.Ns}(2*(p.Nelx+1)-1) = p.axforce;
            case 3
                for row = 1:(p.Nely+1)
                    p.fs{p.Ns}(2*(p.Nelx+1)*row-1) = -p.axforce;
                end
            case 4
                p.fs{p.Ns}(p.Ndof_s(p.Ns)-1) = p.axforce;
                p.fs{p.Ns}(2*(p.Nelx+1)-1) = -p.axforce;
            case 5
                for row = [1 2 p.Nely (p.Nely+1)]
                    p.fs{p.Ns}(2*(p.Nelx+1)*row-1) = -p.axforce;
                end
                for row = [ceil(p.Nely/2) (ceil(p.Nely/2)+1)]
                    p.fs{p.Ns}(2*(p.Nelx+1)*row-1) = 2*p.axforce;
                end
            case 6 % axial force-field 
                i=1;
                %j=1;
                %x=[];
                %y=[];
                if p.nonconforming~=1
                    p.Height=p.Nsy*p.Nely(1)*p.elHeight(1);
                    disp(['Height: ' num2str(p.Height)])
                    p.right_side=p.Nsx*[1:p.Nsy];
                    disp(p.right_side)
                end
                a=-4*(p.axforcefield_max-p.axforcefield_offset)/p.Height^2;
                b=4*(p.axforcefield_max-p.axforcefield_offset)/p.Height;
                k=0;
                for s=1:p.Ns
                    if s==p.right_side(i)
                        for row = 1:p.Nely(s)+1
                            if p.nonconforming ==1
                                y_pos=p.poss(2,s)+(row-1)*p.elHeight(s);
                            else
                                y_pos=k*p.Nely(s)*p.elHeight(1)+(row-1)*p.elHeight(1);
                            end
                            p.fs{s}(2*(p.Nelx(s)+1)*row-1) = a*y_pos^2+b*y_pos+p.axforcefield_offset;
                            %x(j)=y_pos;
                            %y(j)=a*y_pos^2+b*y_pos+p.axforcefield_offset;
                            %j=j+1;
                        end
                        i=i+1;
                        k=k+1;
                    end
                end
                %plot(x,y)
        end
end

end