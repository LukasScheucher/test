function FETIplot( p, FigureHandle, NodalPos, SubsToPlot, fsPost, description, xlim, ylim) %FigureHandle
%fetiPLOT Summary of this function goes here
%   Detailed explanation goes here

    hold on;
    title(description);
    
    %--% set(FigureHandle,'DeleteFcn','doc datacursormode');
    fastPlot = 0;
    if (fastPlot)
        
        if p.nonconforming~=1
            k=1;
        end

        for s = SubsToPlot
            clear X Y Z
            for r = 1:p.Nely(k) % r: row (substructure-local)
                for c = 1:p.Nelx(k) % c: column (substructure-local)
                    % for the bottom nodes, element row = nodes row
                    % for the left nodes, element column = nodes column

                    % left bottom node (local 1)
                    i = (r-1)*(p.Nelx(k)+1)+c; % i: node nr of the left bottom node
                    X(r,c) = NodalPos{s}(2*i-1);
                    Y(r,c) = NodalPos{s}(2*i);
                    % right bottom node: (local 2)
                    i = i+1;
                    X(r,c+1) = NodalPos{s}(2*i-1);
                    Y(r,c+1) = NodalPos{s}(2*i);
                    % right upper node: (local 3)
                    i = i+(p.Nelx(k)+1);
                    X(r+1,c+1) = NodalPos{s}(2*i-1);
                    Y(r+1,c+1) = NodalPos{s}(2*i);
                    % left upper node: (local 4)
                    i = i-1;
                    X(r+1,c) = NodalPos{s}(2*i-1);
                    Y(r+1,c) = NodalPos{s}(2*i);


                end
            end
            Z=zeros(size(X,1),size(X,2));
            hmesh = mesh(X,Y,Z);
            color = [0 101 189]./255;
            set(hmesh,'FaceAlpha',0.5);
            set(hmesh,'EdgeColor',[0 0 0]);
            set(hmesh,'EdgeAlpha',0.5);
            set(hmesh,'FaceColor',color);
            set(hmesh,'LineWidth',0.5); % 1.5 for png pictures
        end
    
    else % (plot with different color)
        for s = SubsToPlot
            if p.nonconforming~=1
                iNsy=int8(s/p.Nsy)+1;
                iNsx=mod(s,p.Nsx);
            end
            %s = (iNsy-1)*p.Nsx + iNsx;

            if isempty(find(SubsToPlot==s,1))
                continue;
            end

            %% Te
            if (p.NoPattern == 0)
                if p.nonconforming==1
                    if mod(s,2)
                        Odd = 1;
                    else
                        Odd = 0;
                    end
                else
                    if mod(iNsy,2) || ~p.ChangeForEvenNsy
                       NsxTest = iNsx;
                    else
                       NsxTest = iNsx+1;
                    end
                    if mod(NsxTest,2)
                       Odd = 1;
                    else
                       Odd = 0;
                    end
                end
            else                   
                if (~isempty(find(s == p.Odd,1)))
                    Odd = 1;              
                else
                    Odd = 0;           
                end

            end
            %%  Plot strains and stresses
            if p.cal_strains==1
                
                %axis equal
                %axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
                for e=1:p.Nelx(s)*p.Nely(s)
                    
                    el_x=p.x_el{s}{e}(1:2:size(p.x_el{s}{e},1))+p.d_el{s}{e}(1:2:size(p.x_el{s}{e},1));
                    el_y=p.x_el{s}{e}(2:2:size(p.x_el{s}{e},1))+p.d_el{s}{e}(2:2:size(p.x_el{s}{e},1));
                    
                    switch p.plot
                        case 'disp'
                            color=p.d_el{s}{e}(p.strain_dir:2:size(p.x_el{s}{e},1));
                        case 'strain'
                            color=p.eps{s}{e}(p.strain_dir,:);
                        case 'stress'
                            color=p.stress{s}{e}(p.strain_dir,:);
                    end
                    disp(color)
                    zz=zeros(size(el_x,1));
                    hmesh = patch(el_x,el_y,color,'CDataMapping','scaled');%,'EdgeColor','k','Marker','o','MarkerFaceColor','k');
                    %set(hmesh,'FaceColor','interp',color);
                end
                %set(hmesh,'EdgeColor',[0 0 0]);
                %set(hmesh,'EdgeAlpha',0.5);
                
                %set(hmesh,'LineWidth',0.5); % 1.5 for png pictures
                colormap(jet(128));                
                
                switch p.plot
                    case 'disp'
                        colorbar('location','eastoutside','YLim',[p.dis_min(p.strain_dir) p.dis_max(p.strain_dir)]);
                    case 'strain'
                        colorbar('location','eastoutside','YLim',[p.strain_min(p.strain_dir) p.strain_max(p.strain_dir)]);
                    case 'stress'
                        colorbar('location','eastoutside','YLim',[p.stress_min(p.strain_dir) p.stress_max(p.strain_dir)]);
                end
                
            else
                for r = 1:p.Nely(s) % r: row (substructure-local)
                    clear X Y Z
                    %--% datatipIndex = 0;
                    for c = 1:p.Nelx(s) % c: column (substructure-local)
                        %--% datatipIndex = datatipIndex + 4;

                        % for the bottom nodes, element row = nodes row
                        % for the left nodes, element column = nodes column

                        % left bottom node (local 1)
                        i = (r-1)*(p.Nelx(s)+1)+c; % i: node nr of the left bottom node
                        %--% datatip{datatipIndex + 1} = ['X: ' num2str(2*i-1) ', Y: ' num2str(2*i)];
                        X(1,c) = NodalPos{s}(2*i-1);
                        Y(1,c) = NodalPos{s}(2*i);
                        % right bottom node: (local 2)
                        i = i+1;
                        %--% datatip{datatipIndex + 2} = ['X: ' num2str(2*i-1) ', Y: ' num2str(2*i)];
                        X(1,c+1) = NodalPos{s}(2*i-1);
                        Y(1,c+1) = NodalPos{s}(2*i);
                        % right upper node: (local 3)
                        i = i+(p.Nelx(s)+1);
                        %--% datatip{datatipIndex + 3} = ['X: ' num2str(2*i-1) ', Y: ' num2str(2*i)];
                        X(2,c+1) = NodalPos{s}(2*i-1);
                        Y(2,c+1) = NodalPos{s}(2*i);
                        % left upper node: (local 4)
                        i = i-1;
                        %--% datatip{datatipIndex + 4} = ['X: ' num2str(2*i-1) ', Y: ' num2str(2*i)];
                        X(2,c) = NodalPos{s}(2*i-1);
                        Y(2,c) = NodalPos{s}(2*i);
                    end

                    % after every row, plot that row
                    Z=zeros(size(X,1),size(X,2));
                    %hmesh = mesh(X,Y,Z);


                    hmesh = patch(surf2patch(X,Y,Z)); % nearly 4x faster than mesh()
                    if mod(s,2)
                        if ( (Odd && ~isempty(find(r == p.SteelRowNrsOddNsx,1))) || (~Odd && ~isempty(find(r == p.SteelRowNrsEvenNsx,1))) )
                            color = [0 101 189]./255; % steel
                            %display(['s=' num2str(s) '; Odd=' num2str(Odd) '; r=' num2str(r) ' --> steel']);
                        else
                            %display(['s=' num2str(s) '; Odd=' num2str(Odd) '; r=' num2str(r) ' --> rubber']);
                            %color = [218 215 203]./255; % rubber
                            color = [186 183 171]./255; % rubber
                            %color = [162 173 0]./255; % rubber green
                        end
                        set(hmesh,'FaceAlpha',0.55);
                    else
                        if ( (Odd && ~isempty(find(r == p.SteelRowNrsOddNsx,1))) || (~Odd && ~isempty(find(r == p.SteelRowNrsEvenNsx,1))) )
                            color = [0 101 189]./255; % steel
                            %display(['s=' num2str(s) '; Odd=' num2str(Odd) '; r=' num2str(r) ' --> steel']);
                        else
                            %display(['s=' num2str(s) '; Odd=' num2str(Odd) '; r=' num2str(r) ' --> rubber']);
                            %color = [218 215 203]./255; % rubber
                            color = [186 183 171]./255; % rubber
                            %color = [162 173 0]./255; % rubber green
                        end
                        set(hmesh,'FaceAlpha',0.5);
                    end
                    set(hmesh,'EdgeColor',[0 0 0]);
                    set(hmesh,'EdgeAlpha',0.5);
                    set(hmesh,'FaceColor',color);
                    set(hmesh,'LineWidth',0.5); % 1.5 for png pictures

                    %--% dcm_obj = datacursormode(FigureHandle);
                    %--% set(dcm_obj,'UpdateFcn',{@myupdatefcn,datatip});
                end
            end
        end
    end
    
    
    
    set(gca,'DataAspectRatio',[1 1 1]);
    az = 0;
    el = 90;
    view(az,el);
    
    set(gca,'XLim',xlim);
    set(gca,'YLim',ylim);
    
    if (p.NoArrows == 0)
    ArrowNr = 0;
        for s = SubsToPlot
            disp(['Subs ' num2str(s)])
            disp(fsPost{s})
            % plot force arrows for this substructure
            % TODO :2: should be :Ndof_n:, rest is aswell hardcoded
            % for Ndof_n = 2 and x,y displ.
            
            %p.elHeight=p.elsize(s);
            if p.nonconforming==1
                j=s;
            else
                j=1;
            end
            ArrowLength = 1.5*p.elHeight(j);
            for i = 1:2:size(fsPost{s},1)
                disp(['i ' num2str(i)])
                if fsPost{s}(i) == 0 && fsPost{s}(i+1) == 0
                    continue;
                end
                % length = sqrt(fsPost{s}(i)^2 + fsPost{s}(i+1)^2);
                length = 1e4;
                direction = ( [fsPost{s}(i) fsPost{s}(i+1)]./length ).*ArrowLength;
                NodeNr = (i+1)/2;
%                 if (fsRevert{s}(NodeNr)==0 || fsRevert{s}(NodeNr)==1)
%                     % fsRevert: 0=Algorithm decides, 1=force Revert,
%                     %          -1=force no Revert
%                     % arrow start must be on the node
%                     start = [NodalPos{s}(2*NodeNr-1) NodalPos{s}(2*NodeNr)];
%                     stop = start + direction;
%                 else
                    stop = [NodalPos{s}(2*NodeNr-1) NodalPos{s}(2*NodeNr)];
                    start = stop - direction;
%                end
                ArrowNr = ArrowNr + 1;
                myarrows{ArrowNr} = arrow([start 0.01],[stop 0.01],6,'BaseAngle',60);
                %[arrowX, arrowY] = dsxy2figxy(gca, [start(1) stop(1)], [start(2) stop(2)]);
                %myarrows{ArrowNr} = annotation('arrow','X',arrowX,'Y',arrowY);
                %set(myarrows{ArrowNr}, 'Units', 'centimeters');
                
    %             diff = start - stop;
    %             display(num2str(sqrt(diff*diff')));
            end
        end
    end
                
    
    hold off;
%     for ArrowNr = 1:length(myarrows)
%         uistack(myarrows{ArrowNr},'top');
%         display(num2str(myarrows{ArrowNr}));
%     end
end

function txt = myupdatefcn(~,event_obj,t)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['I: ',num2str(I)],...
       ['T: ',t{I}]};
end