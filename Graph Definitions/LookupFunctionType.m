classdef LookupFunctionType < handle
    % LookupFunction Type class defines a known constant coefficient Lookup 
    % type using the MATLAB Symbolic toolbox. This class provides a 
    % function to caculate the specific Function and its gradient. This
    % Lookup is intended to be used to capture power flow and capacitence
    % nonlinearities that exist outside the graph modeling framework.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 11/16/2020 - Class creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % potential improvements:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        % minimum properties needed to define a the LookupType object

        % LookupFunctionType calculation (omitting the constant 
        % coefficient) the look equation should be provided in terms of:
        % parameters 'p1',...,'pN' and
        % compatible Symbolic toolbox functions
        F char

    end
    
    properties %(Hidden)
        % Function calculation function
        calcF (1,1)
        % Function jacobian with respect to the [x1, ...,xN, u1, ...,uN]
        calcJ (1,1)
        % number of separate states to consider
        Np (1,1)
    end
    
    methods
        function obj = LookupFunctionType(value,np)
            obj.Np = np;
            obj.F = value;
        end
        
        % update the value of Lookup calculations whenever the power
        % flow property P is changed.
        function set.F(obj,v)
            obj.F = v;
            obj.updateCalcF
        end
    end
    
    methods (Access = private)
        function updateCalcF(obj)
            
            try               
                % Initialize symbolic variables
                p = sym('p',[obj.Np,1]); %parameters
                
                % evalulate the Function equation in symbolic variables
                F_eqn = str2sym(obj.F);
                
                % calculate the Function gradient
                J = jacobian(F_eqn,[p]);
                             
                % Generate Matalab functions
                obj.calcF     = matlabFunction(F_eqn,'Vars',[num2cell(p)]');
                obj.calcJ     = matlabFunction(J,'Vars',[num2cell(p)]');
                                
            catch
                warning('Make sure to define the power flow as a function of x and Symbolic Toolbox compatible functions.')
                obj.calcF = 0;
                obj.calcJ = 0;
            end
        end
        
        
    end
    
    
    
end