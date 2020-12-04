classdef Comp_Graph
    % Comp_Graph is a class for defining a component graph model. The
    % Comp_Graph class should be treated as a parent of a component block
    % model class and a child of the Sys_Graph class. This class adheres to
    % graph model definition standard that is usable by the system 
    % development algorithm. Note that models defined by Comp_Graph are not
    % simulatable. They must first be converted to a Sys_Graph object.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/7/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
    properties  % minimum properties needed to define graph
        
        % Number of Graph State Vertices
        Nv uint16 {mustBeNumeric, mustBeInteger}
        % Number of Graph Edges
        Ne uint16 {mustBeNumeric, mustBeInteger}
        % Number of Graph External Vertices
        Nev uint16 {mustBeNumeric, mustBeInteger}
        % Number of Graph External Edges
        Nee uint16 {mustBeNumeric, mustBeInteger}
        % Number of Graph Inputs
        Nu uint16 {mustBeNumeric, mustBeInteger}
        
        % Graph Edge Matrix:
        % FORMAT:
        % Column 1: Tail Vertex, Column 2: Head Vertex, Row: Edge
        E (:,2) uint16 {mustBeNumeric, mustBeInteger}
        
        % Graph vertex properties
        % Define the properties associated with each graph vertex
        % FORMAT:
        % Vertex(#,1) = Graph_Vertex(Type,Value(Capacitance), Initial Condition);
        % Vertex type defines whether the dynamics are calculate as
        % power flows or generic flows.
        % (1)-Power flow, (2)-Generic flow
        Vertex (:,1) Graph_Vertex        
        
        % Graph Edge Properties
        % Define the properties associated with each graph edge
        % FORMAT:
        % Edge(#,1) = Graph_Edge(Type,Value);
        % Type associates an edge with a Powerflow Type
        % Value is a constant coefficient associated with a powerflow
        % EXAMPLE:
        % Let PType(1).Type = 'u1' and PType(2).Type = 'xt*xh'
        % if Edge(2,1) = Graph_Edge({'1','2'},[3 -5]) then the power
        % flow equation 2 is given as P = 3*u1 -5*xt.
        Edge (:,1) Graph_Edge
        
        % Graph Capacitance Type Definitions
        % Add Comments
        CType
        
        % Graph Power Flow Type Definitions
        % Defines the type of constant conefficient power flow types in the
        % graph definition. The power flow equation should be provided in
        % terms of: Tail State 'xt', Head State 'xh',
        % Inputs 'u1',...,'uN' and compatible Symbolic toolbox functions.
        % FORMAT:
        % PType(Type#).Type = 'String representing powerflow type'
        PType
        
        % Graph Function Type Definitions vertex
        % add comments
        FvType
        
        % Graph Function Type Definitions edge
        % add comments
        FeType
        
    end
    properties
        
        % Graph Incidence Matrix
        M (:,:) {mustBeInteger, mustBeNumeric}
        
        % Input to Edge Mapping Matrix
        B (:,:)
        
        % External Edge to Vertex Matrix
        D (:,:) {mustBeInteger, mustBeNumeric}
        
        % Lookup Function to Vertex Matrix
        Fv (:,:) 
        
        % Lookup Function to Edge Matrix
        Fe (:,:) 
        
        % Capacitance Vector
        C_coeff (:,1) double {mustBeNumeric}
        
        % Initial Condition Vector
        x_init (:,1) double {mustBeNumeric}
        
        % Power Coefficient Matrix
        P_coeff (:,:) double {mustBeNumeric}
        
        % Lookup Function Coefficient Matrix vertex
        Fv_coeff (:,:) double {mustBeNumeric}
        
        % Lookup Function Coefficient Matrix edge
        Fe_coeff (:,:) double {mustBeNumeric}
        
        % Power Coefficient Matrix
        DynType (:,1) double {mustBeInteger}
% % % %         % Vertex to Vertex map lookup 
% % % %         v2Vmap (:,:) double {mustBeNumeric}
% % % %         
% % % %         % Vertex Maps
% % % %         VMap (:,1) Graph_Map
% % % %         
% % % %         % Vertex to Edge map lookups
% % % %         v2Emap (:,:)
% % % %         
% % % %         % Input to Edge map lookups
% % % %         u2Emap (:,:)
% % % %         
% % % %         % Edge Maps
% % % %         EMap (:,1) Graph_Map
    end
    
    methods
       
        % whenever the Edge matrix is updated, update the incidence matrix
        function obj = set.E(obj,v)
            obj.E = v;
            obj = obj.updateM;
        end
        % this function builds the incidence matrix
        function obj = updateM(obj)
            obj.M = zeros(obj.Nv+obj.Nev,obj.Ne);
            for i = 1:obj.Ne
                obj.M(obj.E(i,1),i) = 1;
                obj.M(obj.E(i,2),i) = -1;
            end
        end
        
        % whenever the vertex values change, update the capacitance vector 
        function obj = set.Vertex(obj,v)
            obj.Vertex = v;
            obj = obj.updateVertex;
        end
        % this function builds the Capacitance vector, IC vector, and
        % external edge to state map
        function obj = updateVertex(obj)
            obj.C_coeff    = zeros(obj.Nv+obj.Nev,length(obj.CType));
            for i = 1:obj.Nv+obj.Nev
                for j = 1:length(obj.Vertex(i).Capacitance)
                    obj.C_coeff(i,str2num(obj.Vertex(i).CapType{j})) = obj.Vertex(i).Capacitance(j);
                end
            end
            
            obj.x_init = [obj.Vertex(:).Initial]';
            
            % build the external edge to state map
            obj.D = zeros(obj.Nv+obj.Nev,obj.Nv+obj.Nev);
            for i = 1:length(obj.Vertex)
                obj.D(i,obj.Vertex(i).ExtEdge) = 1;
            end
            
            obj.DynType = [obj.Vertex(:).DynType]';
            obj.Fv_coeff    = zeros(obj.Nv+obj.Nev,length(obj.FvType));

            
            % I don't think we need to generate this until we start to
            % develop the system graph since all component graphs might have
            % different types of lookup maps
% % % %             obj.v2Vmap = zeros(obj.Nv+obj.Nev,obj.Nv+obj.Nev);
% % % %             for i = 1:length(obj.Vertex)
% % % %                 obj.v2Vmap(i,obj.Vertex(i).BkptIndex) = 1;
% % % %             end
% % % %             
% % % %             obj.VMap = [obj.Vertex(:).Map]';
        end

        % when the edge types or values change, update the P_coeff matrix
        % also update the input to edge structure
        function obj = set.Edge(obj,v)
            obj.Edge = v;
            obj = obj.updateP_coeff;
        end
        % the function updates P_coeff matrix
        function obj = updateP_coeff(obj)
            obj.P_coeff      = zeros(obj.Ne,length(obj.PType));
            for i = 1:obj.Ne
                for j = 1:length(obj.Edge(i).Coefficient)
                    obj.P_coeff(i,str2num(obj.Edge(i).PowerType{j})) = obj.Edge(i).Coefficient(j);
                end
            end
           
           % build the input to edge structure
           numU = max(cellfun('size',{obj.Edge(:).Input},2));
           for j = 1:numU
               obj.B.(['B',num2str(j)]) = zeros(obj.Ne,obj.Ne);
               for i = 1:length(obj.Edge)
                   try
                       obj.B.(['B',num2str(j)])(i,obj.Edge(i).Input(j)) = 1;
                   end
               end
           end
           
           
           obj.Fe_coeff    = zeros(obj.Ne,length(obj.FeType));

           % I don't think we need to generate this until we start to
           % develop the system graph since all component graphs might have
           % different types of lookup maps
% % % %            % Mapping from vertex and inputs to table dimension breakpoints
% % % %            for i = 1:length(obj.Edge)
% % % %                ND = length(obj.Edge(i).BkptIndex);
% % % %                for j = 1:ND
% % % %                    if obj.Edge(i).BkptType(j) == 1
% % % %                        obj.v2Emap.(['v2EMap',num2str(ND)])(i,obj.Edge(i).BkptIndex(j)) = 1;
% % % %                    elseif obj.Edge(i).BkptType(j) == 2
% % % %                        obj.u2Emap.(['u2EMap',num2str(ND)])(i,obj.Edge(i).BkptIndex(j)) = 1;
% % % %                    end
% % % %                end
% % % %            end
% % % %            
% % % %            % Edge maps
% % % %            obj.EMap = [obj.Edge(:).Map]';

        end
        
        function plot(obj)
            figure
            G = digraph(obj.E(:,1),obj.E(:,2));
            h = plot(G);
            labeledge(h,obj.E(:,1)',obj.E(:,2)',1:size(obj.E,1));       
        end
        
    end
    
    
        

% Things to add:
% Edge lookup maps, vertex lookup maps

end
