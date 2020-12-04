classdef Graph_Vertex < handle
    % Graph_Vertex class defines the properties relating to a vertex of a
    % dynamic graph. The properties are Type, Value (Capacitance), and
    % initial condition.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/7/2020 - Class creation
    % 11/22/2020 - Updated class constructor method to handle name-value
    %              pair arguements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % minimum properties needed to define the Graph Vertex object
        
        % Define the associated vertex constant coefficient capacitance
        Capacitance (1,1) double {mustBeNumeric} = 0
        
        % Define the associated vertex constant coefficient function
        FcnCoefficient (1,1) double {mustBeNumeric} = 1
        
        % Defines the associated vertex state initial conditions
        Initial (1,1) double {mustBeNumeric} = 0
        
        % Define whether the associatated vertex dyanmics are calculated
        % in terms of power flows or state flows
        DynType    (1,1) uint8  {mustBeNumeric, mustBeInteger} = 1
        
        % Determines how to calculate the capacitance for given vertex
        CapType   = [{'1'}]
        
        % Determines how to calculate the function for given vertex
        FcnType  = [{'1'}];
        
        % Defines the external edge associated with the vertex
        ExtEdge(:,1) uint16 
        
% % % %         % Breakpoint index
% % % %         BkptIndex(1,:) uint16 
% % % %         
% % % %         % Vertex Map
% % % %         Map (1,1) Graph_Map
        
% % % %         % Vertex Lookup Function
% % % %         Fcn (1,1) Lookup_Function
        
        % Lookup Function Breakpoint
        FcnBkpt = []
        
        % Lookup Function Breakpoint Type
        BkptType = []
        
    end
    

    
    methods
        % constructor method
        function obj = Graph_Vertex(varargin)
            
            obj = NameValueDef(obj,varargin{:}); % populate the object with user defined parameters

        end
    end    
    
end