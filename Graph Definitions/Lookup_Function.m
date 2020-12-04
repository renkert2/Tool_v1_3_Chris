classdef Lookup_Function
    % Graph_Map class defines how lookup maps are incorporated into the
    % graph definition
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/24/2020 - Class creation
    % 9/25/2020 - added LookUp and LookUpDeriv functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % add lookup function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % Breakpoint index and type
        % Together, these properties define what states and inputs are
        % breakpoints for the lookup map. BkptIndex in the index of the
        % assoicated breakpoint and BkptType is the type of the associated
        % breakpoint ([1]-state/sink, [2]-input).
        % ex: BkptIndex = [2 3] and BkptType = [1 2] indicates that
        % breakpoint 1 is "graph state 2" and breakpoint 2 is "graph input
        % 3"
        BkptIndex(1,:) uint16 
        BkptType(1,:) uint8 
        
        % Lookup Type
        % This property specifies which lookup map function should be used.
        % Ex: LookupType = 2 inficates that the 2nd LookupFunctionType 
        % object should be used
        LookupType(1,1) uint8 {mustBeInteger}
        
        % Breakpoints to Map
        % This property defines a map such the the breakpoints bkpts for a 
        % lookup table can be given by
        % bkpts = Pts2Map*[x;u];
        Pts2Map(:,:) uint8 
    end
    
        methods
%         function obj = GraphMap
%             % Constructor method
%         end

        
        function val = LookUp(obj,pts)
            % LookUp calculates the value of a Graph_Map lookup table given
            % a set of lookup points
            
            %%% INPUTS
            % obj - a Graph_Map object
            % pts - a cell vector of lookup points (ordered ie. pts = {x1 x2 ... xN})
            
            %%% OUTPUTS
            % val - lookup value
            
    
            % saturate the lookup pts to maintain feasability
            for i = 1:length(pts)
                pts{i} = min(max(pts{i},obj.Bkpts{i}(1)),obj.Bkpts{i}(end));
            end
            
            % lookup the table value
            val = interpn(obj.Bkpts{:},obj.Map,pts{:});
            
            % invert the table value if the property enters by dvision
            if obj.ProOrDiv == 2
                val = 1/val;
            end            
        end
        
        function d_vals = LookUpDeriv(obj,pts,val)
            % LookUpDeriv calculates the derivate of a Graph_Map lookup 
            % table given a set of lookup points
            
            %%% INPUTS
            % obj - a Graph_Map object
            % pts - a cell vector of lookup points (ordered ie. pts = {x1 x2 ... xN})
            % val - value of the lookup table at pts (passed to avoid
            %       redundant lookups)
            
            %%% OUTPUTS
            % d_vals - map gradient w.r.t. each breakpoint (ordered ie. pts = {dg_dx1 ... dg_dxN})
            
            % initialize data 
            delta = inf; % percent to perturb the lookup pts;
            d_vals = zeros(length(pts),1); % initialize output vector
            
            % saturate the lookup pts to maintain feasability
            % also, get perturbation factor (smallest difference between
            % breakpoints on a lookup table axis)
            for i = 1:length(pts)
                pts{i} = min(max(pts{i},obj.Bkpts{i}(1)),obj.Bkpts{i}(end));
                delta = min(min(obj.Bkpts{i}(2:end)-obj.Bkpts{i}(1:end-1)),delta);
            end
            delta = delta/10; % divide perturbation factor by an order of magnitude
            
            % perturb the lookup points...
            pts_left = num2cell([pts{:}]-delta,1); % ... to the left
            pts_rght = num2cell([pts{:}]+delta,1); % ... to the right
            
            % saturate the perturbed lookup pts to maintain feasability
            for i = 1:length(pts)
                pts_left{i} = min(max(pts_left{i},obj.Bkpts{i}(1)),obj.Bkpts{i}(end));
                pts_rght{i} = min(max(pts_rght{i},obj.Bkpts{i}(1)),obj.Bkpts{i}(end));
            end
            
            % get lookup w.r.t each breakpoint
            for i = 1:length(pts)
                
                lookup_pts_left = [pts(1:i-1) pts_left(i) pts(i+1:end)]; % make left lookup points cell vector
                lookup_pts_rght = [pts(1:i-1) pts_rght(i) pts(i+1:end)]; % make right lookup points cell vector
                
                d_left = interpn(obj.Bkpts{:},obj.Map,lookup_pts_left{:}); % left side table value   
                d_rght = interpn(obj.Bkpts{:},obj.Map,lookup_pts_rght{:}); % right side table value
                d_vals(i)  = (d_rght - d_left) / (lookup_pts_rght{i}-lookup_pts_left{i}); % slope of the table
                
                if obj.ProOrDiv == 2
                    d_vals(i) = -d_vals(i)/(val^2); % derivative if lookup enters via division
                end
                
            end
            
          
        end
        
    end
end