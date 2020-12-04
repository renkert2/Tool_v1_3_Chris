classdef CapacitanceType < handle
    % PowerFlowType class defines a known constant coefficient power flow 
    % type using the MATLAB Symbolic toolbox. This class provides a 
    % function to caculate the specific power flow and its gradient
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    % 9/1/2020 - Updated the stored matlab functions to calculate the
    % inverse of the jacobian for easier use in the linearization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % minimum properties needed to define a the Power Flow object

        % Power Flow calculation (omitting the constant coefficient)
        % the power flow equation should be provided in terms of:
        % Tail State 'xt', Head State 'xh', Inputs 'u1',...,'uN' and
        % compatible Symbolic toolbox functions
        C char

    end
    
    properties %(Hidden)
        % Power Flow calculation function
        calcC (1,1)
        % Power flow jacobian with respect to the [tail, head, u1, ...,uN]
        calcJ (1,1)
    end
    
    methods
        function obj = CapacitanceType(value)
            obj.C = value;
        end
        
        % update the value of power flow calculations whenever the power
        % flow property P is changed.
        function set.C(obj,v)
            obj.C = v;
            obj.updateCalcC
        end
    end
    
    methods (Access = private)
        function updateCalcC(obj)
            
            try               
                % Initialize symbolic variables
                x = sym('x'); %state                                 
                c = sym('c'); g = sym('g'); % generic constant coefficient and nonlinear powerflow function terms
                dg_dx = sym('dg_dx'); %generic powerflow function gradients
                
                % evalulate the power flow equation in symbolic variables
                C_eqn = c*str2sym(obj.C);
                
                % calculate the power flow gradient using product rule
%                 J = -((jacobian(C_eqn,x)*g) + (C_eqn*[dg_dx]))/(C_eqn*g)^2;
                J = (jacobian(C_eqn,x)*g) + (C_eqn*dg_dx);
                             
                % Generate Matalab functions
                obj.calcC     = matlabFunction(g*C_eqn,'Vars',[{[c] [x] [g]}]);
%                 obj.calcJ     = matlabFunction(J,'Vars',[{[c] [x] [g] [dg_dx]}]);
                obj.calcJ     = matlabFunction(J,'Vars',[{[c] [x] [g] [dg_dx]}]);
                                
            catch
                warning('Make sure to define the power flow as a function of x and Symbolic Toolbox compatible functions.')
                obj.calcC = 0;
                obj.calcJ = 0;
            end
        end
        
        
    end
    
    
    
end