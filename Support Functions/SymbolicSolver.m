function [Sys] = SymbolicSolver(Sys)
    % SymbolicSolver creates a symbolic representation of a nonlinear graph
    % model. This function relies on the property that the algebraic state
    % equations are linear w.r.t. the algebraic states.
    % provided a DAE graph system in the form
    %               xd_dot = f(x_d,x_a,x_e,u,P_e)
    %                    0 = g(x_d,x_a,x_e,u,P_e)
    %                    y = h(x_d,x_a,x_e,u,P_e)
    % this function provides symbolic representations of the same system in
    % the form
    %               xd_dot = fc(x_d,x_e,u,P_e)
    %                    y = y(x_d,x_e,u,P_e) = [x_d; gc(x_d,x_e,u,P_e)].
    % Note, the output vector is just the state vector [x_d; x_a] expressed
    % in terms of dyanmic states, external states, inputs, and external
    % edges. The symbolic system representation is stored in the Symb
    % parameter of object Sys.
    
    %%% INPUTS
    % Sys - System graph model object

    %%% OUTPUTS
    % Sys - System graph model object
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 11/22/2020 - Function creation - Code adapted from Donald Docimo
    %             (University of Illinois at Urbana-Champaign)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define symbolic variables
idx_x_d = (sum(abs(Sys.C_coeff(1:Sys.Nv)),2) ~= 0);
idx_x_a = (sum(abs(Sys.C_coeff(1:Sys.Nv)),2) == 0);
idx_x_e = Sys.Nv+1:Sys.Nv+Sys.Nev;

x_d     = sym('x_d%d'    ,[sum(idx_x_d)        1]); % dynamic states
x_a     = sym('x_a%d'    ,[sum(idx_x_a)        1]); % algebraic states
x_e     = sym('x_t%d'    ,[length(idx_x_e)     1]); % external states
P_e     = sym('P_e%d'    ,[Sys.Nee             1]); % external edge flows
u       = sym('u%d'      ,[Sys.Nu              1]); % inputs

%% Get state vector filled up with symbolic varialbes
x_full(idx_x_d,1) = x_d;
x_full(idx_x_a,1) = x_a;
x_full(idx_x_e,1) = x_e;

%% Calculate power flows and capacitances
P = CalcP(Sys,x_full,u,ones(Sys.Ne,1)); % calculates power flows
C = CalcC(Sys,x_full); % calcualtes capacitance

%% Solve for the algebraic equations
eqnA(1:sum(idx_x_a),1) = -Sys.M_mod(idx_x_a,:)*P + Sys.D(idx_x_a,:)*P_e == 0; % system of algebraic equations
[A,B] = equationsToMatrix(eqnA,x_a); % convert eqnA to the form Ax=B
x_a_solution = linsolve(A,B); % find solution to the algebraic system

%% Solve for the dynamic equations
eqnD(1:sum(idx_x_d),1) = diag(C(idx_x_d))^-1*(-Sys.M_mod(idx_x_d,:)*P + Sys.D(idx_x_d,:)*P_e); % system of dynamic equations (
x_d_solution = subs(eqnD,x_a,x_a_solution); % plug in the algebraic system solution into the dynamic system equations

%% Create output 
y = [x_full(idx_x_d);x_a_solution]; %y = [x_d; x_a(x_d,x_e,u,P_e)]

%% Store symbolic calculations
Sys.Symb.fc = x_d_solution;
Sys.Symb.y = y;

%% Convert to functions
% % use to convert symbolic expressions into functions
x_d1     = sym('x_d%d'    ,[sum(idx_x_d)        1]); % dynamic states
x_e1     = sym('x_t%d'    ,[length(idx_x_e)     1]); % external states
P_e1     = sym('P_e%d'    ,[Sys.Nee             1]); % external edge flows
u1       = sym('u%d'      ,[Sys.Nu              1]); % inputs

% % use to save the algebraic solution in the graph structure
Sys.solveAlg = matlabFunction(x_a_solution,'Vars',[{[x_d1] [x_e1], [u1], [P_e1]}]); % algrebraic state solution
Sys.solveDyn = matlabFunction(x_d_solution,'Vars',[{[x_d1] [x_e1], [u1], [P_e1]}]); % dynamic state derivative solution
% Sys.solveOut = matlabFunction(y,           'Vars',[{[x_d1] [x_e1], [u1], [P_e1]}]); % system output solution
end


function [P] = CalcP(Sys,x0,u0,Pmap)
    % CalcP calculates the power flows of a graph model.
    
    %%% INPUTS
    % Sys  - System graph model object
    % x0   - state vector
    % u0   - input vector
    % Pmap - vector of lookup map values

    %%% OUTPUTS
    % P - Power flows
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 11/22/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P   = zeros(size(Sys.P_coeff_mod));
xt = Sys.Tails*x0; %tail states
xh = Sys.Heads*x0; %head states
Nu = numel(fieldnames(Sys.B)); % max number of inputs incident per edge
u = cell(1,Nu); % initialize edge input data.
for i = 1:Nu
    u{:,i} = repmat(Sys.B.(['B',num2str(i)])*u0,3,1); % inputs indident per edge
end

% calculate the powerflow along each edge. Note the 3x vector size from
% repmat required to simulate a multi-domain graph
for i = 1:size(Sys.P_coeff_mod,2)
    P(:,i) = Sys.PowerFlow(i).calcP(Sys.P_coeff_mod(:,i),repmat(xt,3,1),repmat(xh,3,1),u{:},repmat(Pmap,3,1));
end

% sum the powerflow coefficients
P = sum(P,2);

end

function [C] = CalcC(Sys,x0)
    % CalcC calculates the capacitance values of a graph model.
    
    %%% INPUTS
    % Sys  - System graph model object
    % x0   - state vector

    %%% OUTPUTS
    % C - Capacitance vector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 11/22/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential improvements 
    % - 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graph lookups
c   = sym(zeros(size(Sys.C_coeff)));
% caculate the capacitance of each vertex for each coefficient
for i = 1:size(Sys.C_coeff,2)
    c(:,i) = Sys.Capacitance(i).calcC(Sys.C_coeff(:,i),x0,1); % the 1 and 0 in these lines will need to be changed
end
c = sum(c,2); % sum across capacitance coefficients

% generic lookups
% initialize data
val = sym(zeros(size(Sys.Fv_coeff)));
Nf = numel(fieldnames(Sys.Fv));
f = cell(1,Nf);
for i = 1:Nf
    f{:,i} = Sys.Fv.(['Fv',num2str(i)])*x0;
end
% caculate the capacitance of each vertex for each lookup coefficient
for i = 1:size(val,2)
    val(:,i) = Sys.Fv_coeff(:,i).*Sys.VertexLookupFunction(i).calcF(f{:});
end
val = sum(val,2); % sum across capacitance coefficients

C = (val(1:Sys.Nv).*c(1:Sys.Nv)); % solve for the capacitance of each vertex

end