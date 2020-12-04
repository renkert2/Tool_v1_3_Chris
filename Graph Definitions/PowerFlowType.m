classdef PowerFlowType < handle
    % PowerFlowType class defines a known constant coefficient power flow 
    % type using the MATLAB Symbolic toolbox. This class provides a 
    % function to caculate the specific power flow and its gradient
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 7/6/2020 - Class creation
    % 9/2/2020 - add constant coefficient to the stored matlab functions
    % and the nonlinear power flow variable to the powerflow calculations 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    % give jacobian 3 outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % minimum properties needed to define a the Power Flow object

        % Power Flow calculation (omitting the constant coefficient)
        % the power flow equation should be provided in terms of:
        % Tail State 'xt', Head State 'xh', Inputs 'u1',...,'uN' and
        % compatible Symbolic toolbox functions
        P char

    end
    
    properties %(Hidden)
        % Power Flow calculation function
        calcP (1,1)
        % Power flow jacobian with respect to the [tail, head, u1, ...,uN]
        calcJ (1,1)
        % number of separate inputs to consider        
        Nu (1,1) 
    end
    
    methods
        function obj = PowerFlowType(value,nu)
            obj.Nu = nu;
            obj.P = value;
        end
        
        % update the value of power flow calculations whenever the power
        % flow property P is changed.
        function set.P(obj,v)
            obj.P = v;
            obj.updateCalcP
        end
    end
    
    methods (Access = private)
        function updateCalcP(obj)
            
            try               
                % Initialize symbolic variables
                xt = sym('xt'); xh = sym('xh'); u = sym('u',[obj.Nu,1]); %state and inputs                                 
                c = sym('c'); g = sym('g'); % generic constant coefficient and nonlinear powerflow function terms
                dg_dxt = sym('dg_dxt'); dg_dxh = sym('dg_dxh'); dg_du = sym('dg_du',[obj.Nu,1]); %generic powerflow function gradients
                
                % evalulate the power flow equation in symbolic variables
                P_eqn = c*str2sym(obj.P);
                
                % calculate the power flow gradient using product rule
                J = (jacobian(P_eqn,[xt;xh;u])*g) + (P_eqn*[dg_dxt dg_dxh dg_du.']);
                             
                % Generate Matalab functions
                obj.calcP     = matlabFunction(g*P_eqn,'Vars',[{[c] [xt] [xh]} num2cell(u)' {[g]}]);
                obj.calcJ     = matlabFunction(J,'Vars',[{[c] [xt] [xh]} num2cell(u)' {[g] [dg_dxt] [dg_dxh]} num2cell(dg_du)' ]);
                
            catch
                warning('Make sure to define the power flow as a function of xt, xh, u1,...,uN and Symbolic Toolbox compatible functions.')
                obj.calcP = 0;
                obj.calcJ = 0;
            end
        end
        
        
    end
    
    
    
end