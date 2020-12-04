function [ConnectV] = ExtractVxConn(Comp,ConnectE)

V = zeros(4,2*2*size(ConnectE,2));

for i = 1:size(ConnectE,2)
    xh1 = Comp(ConnectE{1,i}(1)).E(ConnectE{2,i}(1),2);
    xt1 = Comp(ConnectE{1,i}(1)).E(ConnectE{2,i}(1),1);
    xh2 = Comp(ConnectE{1,i}(2)).E(ConnectE{2,i}(2),2);
    xt2 = Comp(ConnectE{1,i}(2)).E(ConnectE{2,i}(2),1);

%     ConnectV(1,2*i-1:2*i) = ConnectE(1,i);
%     ConnectV(2,2*i-1)     = {[xh1 xh2]};
%     ConnectV(2,2*i)       = {[xt1 xt2]};
        
    V(1,4*i-3:4*i) = repmat(ConnectE{1,i},1,2); % Graph index
    V(2,4*i-3:4*i-2) = [xh1 xh2]; % head index interconnection
    V(2,4*i-1:4*i-0) = [xt1 xt2]; % tail index interconnection
%     V(3,4*i-3:4*i-2) = 2*i-1; % connection number
%     V(3,4*i-1:4*i-0) = 2*i; % connection number
    
end

ind = 1;
for i = 1:2*size(ConnectE,2)
    if V(4,2*i) == 0
        V(4,2*i-1:2*i) = ind;
        idx1 = (V(1, 2*i+1:end) == V(1,2*i-1) & V(2, 2*i+1:end) == V(2,2*i-1));
        idx2 = (V(1, 2*i+1:end) == V(1,2*i-0) & V(2, 2*i+1:end) == V(2,2*i-0));
        idx = idx1 | idx2;
        NumCon = reshape(repmat(i+1:2*size(ConnectE,2),2,1),1,2*(2*size(ConnectE,2)-i));
        rep = NumCon.*idx;
        rep(rep == 0) = [];
        if any(idx)
            V(4,[2*rep 2*rep-1]) = ind;
        end
        ind = ind + 1;
    end
end

V_reduced = flipud(unique(flipud(V)','row')'); %% in improved code, this information should be sufficients for the system generation code

ConnectV = cell(2,ind-1);
for i = 1:ind-1
    idx = V_reduced(4,:) == i;
    ConnectV(:,i)  = {V_reduced(1,idx); V_reduced(2,idx)};    
end


    

end

