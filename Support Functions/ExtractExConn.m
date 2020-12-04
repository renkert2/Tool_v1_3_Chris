function [ConnectE] = ExtractExConn(Comp,Graph)

% ConnectE(:,i)  = {[x y]    ;[E_x_m E_y_n]; Z};
% for x,y,... in {1,2,...,N} for N graphs
% E indicates Edge
% n,m indicate n-th and m-th edges of graphs x and y respectively
% Z is an element of {x,y}
% the first cell indicates which graphs are interconnected
% the second cell indicates which edges of the graphs are equivalent
% the third cell indicates which edge type is resultant in the equivalency

%%%%%%%% Here, "tail and head vertex" refer to the graph representation of
%%%%%%%% the block model

E = zeros(2); % Edge matrix for the block diagram
ind = 1; 
ConnectE = cell(3,size(E,1));

% loop through graphs to develop edge connection sets 
for i = 1:numel(Graph)
    for j = 1:numel(Graph(i).DownVertex)
        xT = i; % tail vertex of block digram
        xH = Graph(i).DownVertex(j); % head vertex of block diagram
        
        if xH ~= 0 % incase there is an unattached vertex
            idxH = find(xT == Graph(xH).Port); % the index of the head vertex port connected to the tail vertex
            idxT = find(xH == Graph(xT).Port); % the index of the tail vertex port connected to the head vertex
            
            portT = [Comp(xT).Edge(:).Port]; % list of ports associated with edges for tail vertex
            portH = [Comp(xH).Edge(:).Port]; % list of ports associated with edges for head vertex
            
            eConn = [find(idxT == portT) find(idxH == portH)]; %[tail vertex port associated edge , head vertex port associated edge]
            
            E(ind,1) = xT; % edge matrix tail
            E(ind,2) = xH; % edge matrix head
            ConnectE(:,ind) = {E(ind,:); eConn ;E(ind,1)};
            ind = ind + 1;
        end
    end
end

figure
G = digraph(E(:,1),E(:,2));
h = plot(G);
labelnode(h,1:numel(Graph),{Graph(:).Name})



end

