classdef MergeManifold < matlab.mixin.Copyable
    % NonpressurizedLiquidTank is a class the defines an
    % Nonpressurized Liquid tank model and its associated graph model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (NonCopyable)
        % minimum properties needed to define a the component object
        % The naming convention is chosen to match the block
        % model naming convention
        
        % Block Name
        Name char
        % Working Fluid
        fluid_props char
        % Number of Inflows
        N_in(1,1) uint16 {mustBeInteger}
        % Number of Outflows
        N_out(1,1) uint16 {mustBeInteger}
        % Component Graph
        Comp_Graph (1,1) Comp_Graph
        
    end
    
    
    
    methods 
        %         function obj = Graph
        %             % Constructor method
        %         end
        
        function obj = generateGraph(obj)
            
            % Number of Graph State Vertices
            obj.Comp_Graph.Nv = 1;
            % Number of Graph Edges
            obj.Comp_Graph.Ne = obj.N_in+obj.N_out;
            % Number of Graph External Vertices
            obj.Comp_Graph.Nev = obj.N_in+obj.N_out;
            % Number of Graph External Edges
            obj.Comp_Graph.Nee = 0;
            % Number of Graph Inputs
            obj.Comp_Graph.Nu = obj.N_in+obj.N_out;
            
            % Graph Edge Matrix:
            % FORMAT:
            % Column 1: Tail Vertex, Column 2: Head Vertex, Row: Edge
            obj.Comp_Graph.E = [[(2:obj.N_in+1)',ones(obj.N_in,1)]; ...
                                [ones(obj.N_out,1),(obj.N_in+2:obj.N_in+obj.N_out+1)']];
            
            % Graph Capacitance Type Definitions
            % Defines the type of constant coefficient capacitance types in the
            % graph definition. The capacitance equation should be provided in
            % terms of: State 'x' and compatible Symbolic toolbox functions.
            % FORMAT:
            % CType(Type#).Type = 'String representing capacitance type'
            CType(1).Type = '1';
            obj.Comp_Graph.CType = CType;
            
            % Graph Power Flow Type Definitions
            % Defines the type of constant conefficient power flow types in the
            % graph definition. The power flow equation should be provided in
            % terms of: Tail State 'xt', Head State 'xh',
            % Inputs 'u1',...,'uN' and compatible Symbolic toolbox functions.
            % FORMAT:
            % PType(Type#).Type = 'String representing powerflow type'
            PType(1).Type = 'xt*u1';
            obj.Comp_Graph.PType = PType;
            
            % Graph Capacitance Type Definitions
            % Defines the type of constant coefficient vertex function types in the
            % graph definition. The capacitance equation should be provided in
            % terms of parameters: 'p1',...,'pN' and compatible Symbolic toolbox functions.
            % FORMAT:
            % FvType(Type#).Type = 'String representing vertex function type'
            FvType(1).Type = '1';
            obj.Comp_Graph.FvType = FvType;
            
            % Graph Lookup Type Definitions
            % Unused. Don't change.
            FeType(1).Type = '1';
            obj.Comp_Graph.FeType = FeType;
            
             % Graph vertex properties
            % Define the properties associated with each graph vertex as
            % arguement-value pairs.
            % Vertex(#,1) = Graph_Vertex('arguement',value)
            % Valid arguments and format are
            % Capacitance - constant coefficient for CType function (vector)
            % FcnCoefficient - constant coefficient for FvType function (vector)
            % Initial - Vertex initial condition
            % DynType - value (1) for  power flow based dynamics or value (2) for generic flow based dynamics
            % CapType - Cell array of characters for capacitance calculation ex:CType(1) = 1 and CType(2) = x and CapType = [{'1' '2'}] => C = c_1 + c_2*x
            % FcnType - Cell array of characters for capacitance function calculation. Same format as CType
            % ExtEdge - index of external edge incident on the vertex (vector)
            % FcnBkpt - index of state or input breakpoints used to calulate the
            % FcnType (vector) ex: FcnType = 'p1*p2' with p1 = x5 and p2 = u3. Then FcnBkpt = [5 3] to represent state 1 and input 2 
            % BkptType - defines whether the breakpoint indices are states (1) or inputs (2) (vector)
            Vertex(1,1) = Graph_Vertex('DynType',1,'Capacitance',0,'CapType',[{'1'}],'FcnType',[{'1'}]);
            for i = 1:obj.N_in+obj.N_out
                Vertex(i+1,1) = Graph_Vertex('FcnCoefficient',0);
            end
            
            % Update the component graph vertex
            obj.Comp_Graph.Vertex = Vertex;
            
            % Graph Edge Properties
            % Define the properties associated with each graph edge as
            % arguement-value pairs
            % Edge(#,1) = Graph_Edge('arguement',value)
            % Valid arguments and format are
            % PowerType - Cell array of characters for power flow calculation ex:PType(1) = xt and PType(2) = xh and PType = [{'1' '2'}] => P = c_1*xt + c_2*xh
            % Coefficient - constant coefficient for PType function (vector)
            % Input - Index for the inputs affecting the edge (vector)
            % Port - Associates an edge with a in/out "port" of a graph (scalar)
            for i = 1:obj.N_in+obj.N_out
                Edge(i,1) = Graph_Edge('PowerType',[{'1'}],'Coefficient',1,'Input',i,'Port',i);
            end
            % Update the component graph edge
            obj.Comp_Graph.Edge = Edge;
            
            % Things to add:
            % Edge lookup maps, vertex lookup maps, nonlinear capacitance
        end
        
    end
end
