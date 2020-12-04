classdef Linear_Graph < handle
    % Linear_Graph is a class for defining a linear graph model. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 9/14/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - Change PType from character structure to sring array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties  % minimum properties needed to define graph
        
        % Number of Graph State Vertices
        A double {mustBeNumeric}
        % Number of Graph Edges
        B double {mustBeNumeric}
        % Number of Graph External Vertices
        W double {mustBeNumeric}
        
    end
        
    
    methods
       
          
    end
    
   

end
