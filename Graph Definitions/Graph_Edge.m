classdef Graph_Edge < handle
    % Graph_Edge class defines the properties relating to a edge of a
    % dynamic graph. The properties are Type and Value
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/7/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % minimum properties needed to define the Graph Vertex object
        
        % Defines which type of power flow is associated with each edge 
        PowerType = [{'1'}]
        
        % Defines the constant coefficient value associated with each edge
        Coefficient (1,:) {mustBeNumeric} = [0]

% % % %         % Define the associated vertex constant coefficient function
% % % %         FcnCoefficient (1,1) double {mustBeNumeric} = 1
        
        % Associates an input(s) index with the graph edge
        Input (1,:) = []
        
        % Associates an edge with a in/out "port" of a graph
        Port (1,:) {mustBeNumeric, mustBeInteger} = 0
        
% % % %         % Determines how to calculate the function for given vertex
% % % %         FcnType  = [{'1'}];
% % % %         
% % % %         % Breakpoint index and type
% % % %         % Together, these properties define what states and inputs are
% % % %         % breakpoints for the lookup map. BkptIndex in the index of the
% % % %         % assoicated breakpoint and BkptType is the type of the associated
% % % %         % breakpoint ([1]-state/sink, [2]-input).
% % % %         % ex: BkptIndex = [2 3] and BkptType = [1 2] indicates that
% % % %         % breakpoint 1 is "graph state 2" and breakpoint 2 is "graph input
% % % %         % 3"
% % % %         BkptIndex(1,:) uint16 
% % % %         BkptType(1,:) uint8 
% % % %         
% % % %         % Edge Map
% % % %         Map (1,1) Graph_Map
        
%         % Edge Lookup Function
%         Fcn (1,1) Lookup_Function
% % % %         
% % % %         % Lookup Function Breakpoint
% % % %         FcnBkpt = []
% % % %         
% % % %         % Lookup Function Breakpoint Type
% % % %         BkptType = []
% % % %         
    end
    

    
    methods
        % an edge object must be initialized with type and value defined
        function obj = Graph_Edge(varargin)
            
            obj = NameValueDef(obj,varargin{:}); % populate the object with user defined parameters

        end
    end    
    
end