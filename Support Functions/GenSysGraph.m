 function [Sys] = GenSysGraph(Comp,ConnectV,ConnectE)
    % GenSysGraph is a function that translates component graphs into a
    % system level graph.
    
    %%% INPUTS
    % Comp - Cell Array of Comp_Graph classes
    % Connect - Information on graph interconnections

    %%% OUTPUTS
    % Sys - System Graph model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/29/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - The calculation for the System D and B matrices may be incorrect
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Initialize the System Graph
    Sys = Sys_Graph;
    
    %% Resolve/Merge Component Graph Informations
    [Comp,Sys.CType,Sys.PType,Sys.FvType] = ResolveComponents(Comp);
    
    %% Generate Component to Global Maps
   
    [V,E,E_mod]= Comp2SysMaps(Comp,ConnectV,ConnectE);   
    Sys.Nv = sum([Comp.Nv]);

    C_coeff = V*vertcat(Comp.C_coeff);
    Fv_coeff = V*vertcat(Comp.Fv_coeff);
    M = V*blkdiag(Comp.M)*E;
    x_init = V*vertcat(Comp.x_init);
    D = V*blkdiag(Comp.D)*V'; %
    D(:,~any(D,1)) = [];
    DynType = V*vertcat(Comp.DynType)./vecnorm(V,1,2);
    
    idx_d = C_coeff ~= 0;
    idx_a = C_coeff == 0;
    
    Sys.M       = [M(idx_d,:); M(idx_a,:)];
    Sys.C_coeff = [C_coeff(idx_d,:); C_coeff(idx_a,:)];
    Sys.Fv_coeff = [Fv_coeff(idx_d,:); Fv_coeff(idx_a,:)];
    Sys.x_init  = [x_init(idx_d,:); x_init(idx_a,:)];
    Sys.D       = [D(idx_d,:); D(idx_a,:)];
    Sys.DynType = [DynType(idx_d,:); DynType(idx_a,:)];

    Sys.P_coeff = E_mod'*vertcat(Comp.P_coeff);
    
    %%%%%%%%%%%%%%% START CALCULATE B_Comp to B_Sys %%%%%%%%%%%%
    B_Comp = [Comp.B];
    
    EB   = cell(numel(fieldnames(Comp(1).B)),1);
    uInt = cell(numel(fieldnames(Comp(1).B)),1);
    for i = 1:numel(fieldnames(Comp(1).B))
        BDiag = blkdiag(B_Comp.(['B',num2str(i)]));
        EB{i} = spones(E)'*BDiag;
        uIntIdx = EB{i}.*(1:size(BDiag,2)); % relate comp input index to sys edges
        uIntIdx(~(sum(EB{i},2) > 1),:) = []; % remove non coupled edges/inputs
        v = nonzeros(uIntIdx'); % indicies of zero elements to remove
        uInt{i} = reshape(v,2,size(uIntIdx,1))'; % remove zero elements form the matrix
    end
    
    uConnect = vertcat(uInt{:});
    uConnect = [uConnect zeros(size(uConnect,1),1)];
    ind = 1;
    for i = 1:size(uConnect,1)
        
        if uConnect(i,3) == 0
            uConnect(i,3) = ind;
            ind = ind + 1;   
            go = 1;
        
        
        srch = uConnect(i,[1 2]);
        while go <= length(srch)
            locs = [(uConnect(:,[1,2])==srch(go) & uConnect(:,3) == 0)];
            if sum(locs,'All') > 0
                uConnect(or(locs(:,1),locs(:,2)),3) = uConnect(i,3);
                locs(~or(locs(:,1),locs(:,2)),:) = 1;
                srch = [srch reshape(uConnect(~locs),1,[])];
            end
            go = go + 1;
        end
        end   
             
        
    end
    
    Uc = zeros(size(BDiag,1),ind-1);
    for i = 1:ind-1
        Uc(reshape(uConnect((uConnect(:,3) == i),1:2),1,[]),i) = 1;
    end
    
    Ucbar = zeros(size(BDiag,1),size(BDiag,1)-sum(sum(Uc)));
    idx_cbar = ((sum(Uc,2) ~= 1)).*(1:size(BDiag,1))';
    idx_cbar(idx_cbar == 0) = [];
    Ucbar(sub2ind(size(Ucbar),idx_cbar',1:numel(idx_cbar))) = 1;
    
    U = [Ucbar, Uc];
    
    B = cell(numel(fieldnames(Comp(1).B)),1);
    for i = 1:numel(fieldnames(Comp(1).B))
        B{i} =  EB{i}*U;
    end
    
    B_stacked = vertcat(B{:});
    uEmpty = ~any(B_stacked,1);
    
    for i = 1:numel(fieldnames(Comp(1).B))
        Sys.B.(['B',num2str(i)]) = full(spones(B{i}(:,~uEmpty)));
    end

    %%%%%%%%%%%%%%% END CALCULATE B_Comp to B_Sys %%%%%%%%%%%%
   
    
    %%%%%%%%%%%%%%% START CALCULATE Vertex_Map_Comp to Vertex_Map_Sys %%%%%%%%%%%%

% % % %     AllMaps = [];
% % % %     Map2V_ = [];
% % % %     for i = 1:numel(Comp)
% % % %         AllMaps = [AllMaps; vertcat(Comp(i).Vertex(:).Map)];
% % % %         Map2V = zeros(Comp(i).Nv+Comp(i).Nev,Comp(i).Nv+Comp(i).Nev);
% % % %         for j = 1:Comp(i).Nv+Comp(i).Nev
% % % %             if ~isempty(Comp(i).Vertex(j).Map.BkptIndex)
% % % %                 Map2V(j,Comp(i).Vertex(j).Map.BkptIndex) = 1:length(Comp(i).Vertex(j).Map.BkptIndex);
% % % %             end
% % % %         end
% % % %         Map2V_ = blkdiag(Map2V_,Map2V);
% % % %     end
% % % %     Locs = ([AllMaps.ProOrDiv] ~= 0);
% % % %     MapLocs = zeros(length(AllMaps),1);
% % % %     MapLocs(Locs') = (1:sum(Locs))';
% % % %     
% % % % 
% % % %     Map2V = eye(Sys.Nv+Sys.Nev); % since map_i impacts vertex_i
% % % %     V2Map = V*Map2V_*V';
% % % %     MapLocs_new = V*MapLocs;
% % % %     V2Map = [V2Map(idx_d,idx_d), V2Map(idx_d,idx_a); V2Map(idx_a,idx_d), V2Map(idx_a,idx_a)];
% % % %     Map2V = [Map2V(idx_d,idx_d), Map2V(idx_d,idx_a); Map2V(idx_a,idx_d), Map2V(idx_a,idx_a)];
% % % %     MapLocs_new = [MapLocs_new(idx_d); MapLocs_new(idx_a)];
% % % %     
% % % %     % remove extra maps
% % % %     idx_del = ~any(V2Map,2);
% % % %     Map2V(idx_del,:) = [];
% % % %     V2Map(idx_del,:) = [];
% % % %     MapLocs_new(idx_del) = [];
% % % %     
% % % %     Maps = AllMaps([AllMaps.ProOrDiv] ~= 0);
% % % %     Maps = Maps(MapLocs_new);
% % % %     Sys.vMaps = Maps;
% % % %     
% % % %     for i = 1:numel(Maps)
% % % %         Pts2Map = zeros(max(V2Map(i,:)),size(V2Map(i,:),2));
% % % %         [a,idx] = sort(V2Map(i,:));
% % % %         idx(a==0) = [];
% % % %         Pts2Map(sub2ind(size(Pts2Map),1:numel(idx),idx)) = 1;
% % % %         Sys.vMaps(i).Pts2Map = Pts2Map;
% % % %     end
% % % %     Sys.Map2v = Map2V';
    
    %%%%%%%%%%%%%%% END CALCULATE Vertex_Map_Comp to Vertex_Map_Sys %%%%%%%%%%%%

    
    %%%%%%%%%%%%%%% START CALCULATE VERTEX LOOKUP FUNCTION %%%%%%%%%%%%%%%
    Fv_Comp = [Comp.Fv];
%     F = cell(numel(fieldnames(Comp(1).F),1));
    for i = 1:numel(fieldnames(Comp(1).Fv))
        FvDiag = blkdiag(Fv_Comp.(['Fv',num2str(i)]));
        FvSys  = V*FvDiag*V';
%         F{i} = [FSys(idx_d,idx_d), FSys(idx_d,idx_a); FSys(idx_a,idx_d), FSys(idx_a,idx_a)];
        Sys.Fv.(['Fv',num2str(i)]) = [FvSys(idx_d,idx_d), FvSys(idx_d,idx_a); FvSys(idx_a,idx_d), FvSys(idx_a,idx_a)];
    end
%     F_stacked = vertcat(F{:}); 
%     fEmpty = ~any(F_stacked,1);
%     for i = 1:numel(fieldnames(Comp(1).F))
%         Sys.F.(['F',num2str(i)]) = full(spones(F{i}(:,~fEmpty)));
%     end
    
    %%%%%%%%%%%%%%% START CALCULATE VERTEX LOOKUP FUNCTION %%%%%%%%%%%%%%%

    
    Sys = MakeModifiedGraph(Sys);
    
    for i = 1:numel(Sys.CType)
        Sys.Capacitance(i) = CapacitanceType(Sys.CType(i).Type);
    end
    
    for i = 1:numel(Sys.PType)
        Sys.PowerFlow(i) = PowerFlowType(Sys.PType_mod(i).Type,numel(fieldnames(Sys.B)));
    end
    
    for i = 1:numel(Sys.FvType)
        Sys.VertexLookupFunction(i) = LookupFunctionType(Sys.FvType(i).Type,numel(fieldnames(Sys.Fv)));
    end
        

     %%% PLOTING
    plot(Sys)
    
    
    
end


function [Sys] = MakeModifiedGraph(Sys)
    % MakeModifiedGraph creates a modified graph as defined in Christopher
    % Aksland's Master's Thesis. This function improves on the original
    % definition by reducing the repmat() size of the modified graph to 3 
    % as opposed to the arbitrary 7 mentioned in the thesis.
    
    %%% INPUTS
    % Sys - System graph

    %%% OUTPUTS
    % Sys - System graph with modified graph matrices updated
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 8/31/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - reduced the size of the modified graph for systems where DynType is
    %   the same for all vertices
    % - try to remove the for loops
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate Modified Incidence matrix
    % type 1 elements' power flows are expressed in terms of power
    % type 2 elements' power flows are expressed in terms of generic flows
    Type1 = Sys.DynType == 1; % index of type 1 vertices
    Type2 = ~Type1; % index of type 2 vertices
    M_full = Sys.M; % Full incidence matrix
    M_head = -Sys.Heads'; % Heads portion of the incidence matrix
    M_tail = Sys.Tails'; % Tails portion of the incidence matrix
    
    M_full(Type2,:) = 0; % remove all type 2 elements from the full matrix
    M_head(Type1,:) = 0; % remove all type 1 elements
    M_tail(Type1,:) = 0; % remove all type 1 elements
    
    Sys.M_mod = [M_full M_head M_tail]; % Modified System Incidence Matrix
    
    
    % Generate Modified Power Flow Coefficient Matrix
    % the power flow vector will be replicated 3 times such that
    % P_mod = [P; P/xh; P/xt];
    E_coeff_full = Sys.P_coeff;
    E_coeff_head = Sys.P_coeff;
    E_coeff_tail = Sys.P_coeff;
    
    % loop through each graph edge to determine vertex types adjacent to
    % each edge. 
    for i = 1:Sys.Ne
        if sum(Sys.DynType(Sys.E(i,:))) == 4 % if there are two type 2 vertices on the edge, full power info is deleted 
            E_coeff_full(i,:) = 0;      
        elseif sum(Sys.DynType(Sys.E(i,:))) == 3 % one type 1 and one type 1 vertex on the edge
            idx = Sys.DynType(Sys.E(i,:)) == 2;
            if idx(1) 
                E_coeff_head(i,:) = 0; % delete head power coeff if the tail is type 2.
            else
                E_coeff_tail(i,:) = 0; % delete tail power coeff if the head is type 2.
            end
        else % if there are two type 1 vertices on the edge, delete the head and tail flow coeffs
            E_coeff_head(i,:) = 0;
            E_coeff_tail(i,:) = 0;
        end
    end
    
    % remove unused edge calculations to reduce the size of the simulated
    % system
    E_coeff_all = blkdiag(E_coeff_full,E_coeff_head,E_coeff_tail); % All power flow coefficients
    i_empty = all(E_coeff_all == 0,1); % unused power flow types (ie. empty columns of E_coeff_all
    
    % convert the Power Flow types to symbolic variables
    Ptype = [];
    for i = 1:numel(Sys.PType)
        Ptype = [Ptype; str2sym(Sys.PType(i).Type)];
    end
    xt = sym('xt'); xh = sym('xh'); % create tail and head symbolic variables for power flow reduction
    Ptype_all = [Ptype; Ptype/xh; Ptype/xt]; % create the full power flow type vector
    Ptype_all(i_empty) = []; % removed unused power flow types
    E_coeff_all(:,i_empty) = []; % removed unused power flow coefficients
    
    [Ptype_new, IC] = UniqueString(Ptype_all); % find the unique power flow representations
    
    % Create a map that condesnses E_coeff_all based on which power flow
    % types are unique
    all2new = zeros(numel(Ptype_all),numel(Ptype_new))'; 
    all2new(sub2ind(size(all2new),IC',1:numel(IC))) = 1; 
    
    E_coeff_new = E_coeff_all*all2new'; % condensed E_coeff matrix
    
    Sys.P_coeff_mod = E_coeff_new; % Modified power flow coefficient matrix
    Sys.PType_mod = Ptype_new; % Modified power flow type

end


function [Comp,CType,PType,FvType] = ResolveComponents(Comp)
    % ResolveComponents is a function that resolves differences in the
    % definitions of component graphs. Theses differences include
    % Capacitance Type, Power Flow Type, number of inputs incident per
    % edge, and number of lookup maps. This function will make sure that
    % component graphs are all defined in the same manner.
    
    %%% INPUTS
    % Comp - Cell Array of Comp_Graph classes

    %%% OUTPUTS
    % Comp - Cell Array of updated Comp_Graph classes
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/29/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    VAll = vertcat(Comp(:).Vertex);
    EAll = vertcat(Comp(:).Edge);

        
    % Resolve Capacitance Type
    CTypeAll = [Comp(:).CType];
    [CType, IC] = UniqueString(CTypeAll);
    Comp = MakeCoeffMatrix(Comp, IC,'C',{VAll(:).CapType},{VAll(:).Capacitance});
    
    % Resolve Power Flow Type
    PTypeAll = [Comp(:).PType];
    [PType, IC] = UniqueString(PTypeAll);
    Comp = MakeCoeffMatrix(Comp, IC,'P',{EAll(:).PowerType},{EAll(:).Coefficient});

    % Resolve Vertex Lookup Function Type
    FvTypeAll = [Comp(:).FvType];
    [FvType, IC] = UniqueString(FvTypeAll);
    Comp = MakeCoeffMatrix(Comp, IC,'Fv',{VAll(:).FcnType},{VAll(:).FcnCoefficient});

    % Resolve Edge Lookup Function Type
% % % %     FeTypeAll = [Comp(:).FeType];
% % % %     [FeType, IC] = UniqueString(FeTypeAll);
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %     % assume fcn coeff is 1 for all edges! May want to change that
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %     Comp = MakeCoeffMatrix(Comp, IC,'Fe',{EAll(:).FcnType},num2cell(ones(1,length({EAll(:).Coefficient})))); 
        
    % find the max number of inputs incident on an edge
    Nu = max(cellfun('size',{EAll(:).Input},2));
    for i = 1:length(Comp)
        for j = 1:Nu
               Comp(i).B.(['B',num2str(j)]) = zeros(Comp(i).Ne,Comp(i).Ne);
               for k = 1:length(Comp(i).Edge)
                   try
                       Comp(i).B.(['B',num2str(j)])(k,Comp(i).Edge(k).Input(j)) = 1;
                   end
               end
        end       
    end
    
    % find the max number of lookups on a vertex
    Nf = max(cellfun('size',{VAll(:).FcnBkpt},2));
    for i = 1:length(Comp)
        for j = 1:Nf
               Comp(i).Fv.(['Fv',num2str(j)]) = zeros(Comp(i).Nv+Comp(i).Nev,Comp(i).Nv+Comp(i).Nev);
               for k = 1:length(Comp(i).Vertex)
                   try
                       Comp(i).Fv.(['Fv',num2str(j)])(k,Comp(i).Vertex(k).FcnBkpt(j)) = 1;
                   end
               end
        end       
    end
    
% % % %     % find the max number of lookups on a Edge
% % % %     Nf = max(cellfun('size',{EAll(:).FcnBkpt},2));
% % % %     for i = 1:length(Comp)
% % % %         for j = 1:Nf
% % % %                Comp(i).Fe.(['Fe',num2str(j)]) = zeros(Comp(i).Nv+Comp(i).Nev,Comp(i).Nv+Comp(i).Nev);
% % % %                for k = 1:length(Comp(i).Vertex)
% % % %                    try
% % % %                        Comp(i).Fv.(['Fv',num2str(j)])(k,Comp(i).Vertex(k).FcnBkpt(j)) = 1;
% % % %                    end
% % % %                end
% % % %         end       
% % % %     end
    
    % update this to include lookup tables at some point
 
   
end

function [strUnique,IC] = UniqueString(strAll)
    % UniqueString is a function that finds unique strings represent
    % functions. This operates differently from strcmp. For example strcmp
    % would say that'u*x' ~= 'x*u' while UniqueString would say that 
    % 'u*x' == 'x*u'. This function operates be converting the strings to
    % symbolic variables and then establishing whether the results
    % equations are equal.
    
    %%% INPUTS
    % strAll - structure with strings or symbolic variables stored in the Type variable

    %%% OUTPUTS
    % strUnique - structure with unique strings stored in the Type variable
    % IC - Index vector that maps strUnique to strAll 
    %      (ie. strAll = strUnique(IC).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/29/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % convert the strings to symbolic equations if needed
    try
        allSym = [];
        for i = 1:length(strAll)
            allSym = [allSym; str2sym(strAll(i).Type)];
        end
    catch
        allSym = strAll;
    end
    allSymUnique = allSym; % we will use this variable to find the unique expressions
    
    % loop through the vector of symbolic equations and remove duplicate information
    i = 1;
    while i < length(allSymUnique)
        idxRep = [false(i,1); isAlways([allSymUnique(i)] == [allSymUnique(i+1:end)],'Unknown','false')];
        allSymUnique(idxRep) = []; % remove non-unique data
        i = i + 1;
    end
    
    % get the indices that map unique strings to all strings
    IC = zeros(length(allSym),1);
    for i = 1:length(allSymUnique)
        idx = isAlways([allSymUnique(i)] == [allSym],'Unknown','false');
        IC(idx) = i;
    end
    
    % output
    strArray = arrayfun(@char, allSymUnique, 'uniform', 0);
    strUnique = cell2struct(strArray,'Type',2);


end

function [Comp] = MakeCoeffMatrix (Comp,ic,indicator,type,val)
    % MakeCoeffMatrix is a function develops coefficient matrices for the
    % component graphs. Information stored in as Type and Value is
    % converted to matrices that are used the system generation code
    
    %%% INPUTS
    % Comp - Cell Array of Comp_Graph classes
    % ic - vector map of unique types to all types. This is the output of 
    %      the UniqueString function
    % indicator - String prefix for the type and coefficient matrix 'C' and
    %             'P' are currently supported
    % type - cell array of types
    % val - cell arracy of values

    %%% OUTPUTS
    % Comp - Cell array of Comp_Graph classes with updated coefficient
    % matrices
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/29/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize indices for height and width of the component coefficient
    iHT = 1;
    iWD = 1;
    
    % loop through and update coefficient matrix of component graphs
    for i = 1:length(Comp)
        ht    = size(Comp(i).([indicator,'_coeff']),1); % matrix height
        wd    = length(Comp(i).([indicator,'Type'])); % number of component types
    
        IC   =   ic(iWD:iWD+wd-1); % map from component to system type
        Type = type(iHT:iHT+ht-1); % component type
        Val  =  val(iHT:iHT+ht-1); % component value
        Comp(i).([indicator,'_coeff']) = zeros(ht,max(ic)); %initialize matrix size
        % fill the coefficient matrix using Type and Value information
        for j = 1:ht
            for k = 1:length(Val{j})
                Comp(i).([indicator,'_coeff'])(j,IC(str2num(Type{j}{k}))) = Val{j}(k);
            end
        end
        
        % update indicies
        iHT = iHT + ht;
        iWD = iWD + wd;
    end
    
end

function [V_op,E_op,E_op_mod] = Comp2SysMaps(Comp,ConnectV,ConnectE)

    e_tot = 0;
    v_tot = 0;
    Ne_mat = zeros(1,length(Comp));
    Nv_mat = zeros(1,length(Comp));
    chi_ubar = [];
    chi_lbar = [];
    for i = 1:length(Comp)
        Nv_mat(i) = v_tot;
        Ne_mat(i) = e_tot;
        chi_ubar = [chi_ubar; (1:1:Comp(i).Nv)'+v_tot];
        chi_lbar = [chi_lbar; (Comp(i).Nv+1:1:Comp(i).Nv+Comp(i).Nev)'+v_tot];
        v_tot = v_tot + Comp(i).Nv + Comp(i).Nev;
        e_tot = e_tot + Comp(i).Ne;
    end
    
    
chi = (1:1:v_tot)';
Xi = (1:1:e_tot)';

delta = cell(size(ConnectV,2),1);
delta_ubar = {};
delta_lbar = {};
chi_hat = [];
for i = 1:size(ConnectV,2)
    delta{i} = ConnectV{2,i} + Nv_mat(ConnectV{1,i});
    if isempty(my_setdiff(delta{i},chi_lbar))
        delta_lbar = [delta_lbar; delta{i}];
    else
        delta_ubar = [delta_ubar; delta{i}];
    end 
    chi_hat = [chi_hat; delta{i}'];
end

sigma = cell(size(ConnectE,2),1);
sigma_mod = cell(size(ConnectE,2),1);
Xi_hat = [];
for i = 1:size(ConnectE,2)
    sigma{i} = ConnectE{2,i} + Ne_mat(ConnectE{1,i});
    Xi_hat = [Xi_hat; sigma{i}'];
    sigma_mod{i} = sigma{i}(ConnectE{1,i} == ConnectE{3,i});
end

% for loops can be removed in future interations for speed up
ind = my_setdiff(chi_ubar,chi_hat);
V_sbar_cbar = zeros(v_tot,length(ind));
for i = 1:length(ind)
    V_sbar_cbar(ind(i),i) = 1;
end

V_sbar_c = zeros(v_tot,length(delta_ubar));
for i = 1:length(delta_ubar)
    V_sbar_c(delta_ubar{i},i) = ones(length(delta_ubar{i}),1);
end

V_s_c = zeros(v_tot,length(delta_lbar));
for i = 1:length(delta_lbar)
    V_s_c(delta_lbar{i},i) = ones(length(delta_lbar{i}),1);
end

ind = my_setdiff(chi_lbar,chi_hat);
V_s_cbar = zeros(v_tot,length(ind));
for i = 1:length(ind)
    V_s_cbar(ind(i),i) = 1;
end

ind = my_setdiff(Xi,Xi_hat);
E_cbar = zeros(e_tot,length(ind));
for i = 1:length(ind)
    E_cbar(ind(i),i) = 1;
end

E_c = zeros(e_tot,length(sigma));
for i = 1:length(sigma)
    E_c(sigma{i},i) = .5*[1; 1];
end

E_c_mod = zeros(e_tot,length(sigma_mod));
for i = 1:length(sigma_mod)
    E_c_mod(sigma_mod{i},i) = 1;
end
V_op = [V_sbar_cbar V_sbar_c V_s_c V_s_cbar]';
E_op = [E_cbar E_c];
E_op_mod = [E_cbar E_c_mod];


end