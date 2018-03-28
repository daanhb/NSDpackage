function [X, W] = choosePathV2(a,b,criticalPoints, CritPointEndPoint, FPindices, P, G, freq, Npts, numPathsSD, visuals, X_, W_, hFinite, R, ainf, binf)
       
    small=1E-14; %this amount needs to be small, but not too small.
    %All edges need to have some length, but if that length is negligable,
    % some silly connected paths can occur, especially when there are finite paths,
    % as an SD path can end at the beginning at another
    
    nodeCount=0;
    pathCount=0;
    if ainf
        nodeCount=1;
        node(nodeCount)=struct('z',a*R,'type','Q','CPindex',[],'infStart',true,'finiteEnd',[]);
        startNodes=1;
    else
        startNodes=[];
    end
    endNodes=[];
    
    %so many indices...
    for j=1:length(criticalPoints)
            nodeCount=nodeCount+1;
            node(nodeCount)=struct('z',criticalPoints(j),'type','C','CPindex',j,'infStart',false,'finiteEnd',[]);
            
            if ~ainf && node(nodeCount).z==criticalPoints(1)
                startNodes=[startNodes nodeCount];
            end
            if ~binf && node(nodeCount).z==criticalPoints(end)
                endNodes=[endNodes nodeCount];
            end
            
            %outerNodeCount=nodeCount; %keep track of node where the SD path starts
        for m=1:length(CritPointEndPoint{j})
            %nodeCount=nodeCount+1;
            pathCount=pathCount+1;
            if isempty(hFinite{j,m}) %if quad point isn't also another critical point (can happen for finite paths)
                nodeCount=nodeCount+1;
                node(nodeCount)=struct('z',CritPointEndPoint{j}(m),'type','Q','CPindex',j,'infStart',false,'finiteEnd',[]);
            else 
                nodeCount=nodeCount+1;
                node(nodeCount)=struct('z',CritPointEndPoint{j}(m), 'type', 'Q', 'CPindex', [j FPindices{j}(m)], 'infStart',false, 'finiteEnd', hFinite{j,m}(end));
            end
           
        end
    end
    if binf
        nodeCount=nodeCount+1;
        node(nodeCount)=struct('z',b*R,'type','Q','CPindex',[],'infStart',true,'finiteEnd',[]);
        endNodes=nodeCount;
    end
    
    A=inf(nodeCount);
    for j=1:nodeCount
       for ell=(j+1):nodeCount
           if strcmp(node(j).type,'C')
               if node(j).z==node(ell).z %nodes represent same critical point
                   A(ell,j)=small;
               elseif ismember(node(j).CPindex, node(ell).CPindex) %connected to same critical point
                   A(ell,j)=small;
               end
           elseif strcmp(node(j).type,'Q')
               if strcmp(node(ell).type,'Q')
                   A(ell,j)=pathCost( node(j).z, node(ell).z, G{1}, freq, Npts);
               elseif ismember(node(ell).CPindex, node(j).CPindex) 
                   A(ell,j)=small;
               end
           end
           A(j,ell)=A(ell,j);
       end
    end
    
    Gr = graph(A.','upper');
    cost=inf;
    for j=startNodes
        for ell=endNodes
            [path_, cost_] =shortestpath(Gr,j,ell);
            if cost_<cost
                cost=cost_;
                path=path_;
            end
        end
    end
    
    if cost>.1
        warning(sprintf('Cost of truncating SD path is %e, may not be optimal path',cost));
        warning('Using standard quadrature instead :-(')
        [X_,W_] = abortPathSearch(a,b,ainf,binf,freq,Npts);
        return;
    end
    
    %check for three consecutice Q types, doesn't make sense but can happen
    %in practice, delete the middle Q
    tripleCheck=true;
    while tripleCheck
        tripleCheck=false;
        for j=2:(length(path)-1)
            if strcmp(node(path(j-1)).type,'Q') && strcmp(node(path(j)).type,'Q') && strcmp(node(path(j+1)).type,'Q')
                path=[path(1:(j-1)) path((j+1):end)];
                tripleCheck=true;
                break;
            end
        end
    end
    
    %exlude the endpoints from the path order, as they were only artificial
    if ainf
       path=path(2:end);
    end
    if binf
       path=path(1:(end-1));
    end
    %also want to ignore the second path of any C-Q-C, as Q-C will be a
    %finite path
    finitePathNodesUsed=[];
    for j=2:(length(path)-1)
        if strcmp(node(path(j-1)).type,'C') && strcmp(node(path(j)).type,'Q') && strcmp(node(path(j+1)).type,'C')
            if length(node(path(j)).CPindex)<2
                error('The C-Q-C combo should mean Q is a finite path... something is wrong');
            end
            %determine if finite path starts or ends at this critical
            %point:
            if node(path(j)).CPindex(1)==node(path(j-1)).CPindex
                finitePathNodesUsed=[finitePathNodesUsed j];
            else
                finitePathNodesUsed=[finitePathNodesUsed j-1];
            end
        end
    end
    
    QQPathNodesUsed=[];
    for j=1:(length(path)-1)
        if  strcmp(node(path(j)).type,'Q') && strcmp(node(path(j+1)).type,'Q')
            QQPathNodesUsed=[QQPathNodesUsed j];
        end
    end
    
    %now construct pairs of nodes corresponding to each SD path actually
    %used:
    Pvals=[];
    for j=1:(length(path)-1)
        if ~ismember(j,[finitePathNodesUsed QQPathNodesUsed])
            Pvals=[Pvals; node(path(j)).z node(path(j+1)).z ];
        end
    end
    
    X=[];   W=[]; xv=[]; yv=[]; finalPath=[];
    for Pval = Pvals.'%[path(1:(end-1)).' path(2:end).']
        %Pval=[node(pathNodes(1)).z node(pathNodes(2)).z];
        if  ismember(Pval.',P,'rows')
            [~,ind]=ismember(Pval.',P,'rows');
            inOut=1;
        elseif ismember(fliplr(Pval.'),P,'rows')
            [~,ind]=ismember(fliplr(Pval.'),P,'rows');
            inOut=-1;
        else
            warning('Steepest descent path not found, using standard quadrature instead :-(');
            [X_,W_] = abortPathSearch(a,b,ainf,binf,freq,Npts);
            return;
        end
        W=[W; inOut*W_{ind}];
        X=[X; X_{ind};];   
        finalPath= [finalPath ind];
    end

    if visuals
        for fullIndex=1:numPathsSD
            if ismember(fullIndex,finalPath)
                plot(X_{fullIndex},'xm'); 
            else
                plot(X_{fullIndex},'oc'); 
            end
        end
        if ainf
            plot(a*R+0.0000000001i,'ro');
        end
        if binf
            plot(b*R+0.0000000001i,'ro');
        end
        %plot([a b],[0 0], 'b', 'LineWidth',2);
        if ~isempty(R)
            t=linspace(0,2*pi,1000);
            plot(R*exp(1i*t),'r--');
            plot(P(:,1)+0.0000000001i,'ko')
        end
    end
    
end

