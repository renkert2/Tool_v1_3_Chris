%%
function Sys = LinearizeGraph(Sys,x0,u0)
    % LinearizeGraph linearizes the graph Sys about states x0 and inputs u0
    
    %%% INPUTS
    % Sys - System Graph to linearize
    % x0 - state linearization point (including sink states)
    % u0 - input linearization point

    %%% OUTPUTS
    % Sys - System Graph model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 9/1/2020 - Function creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% lookup maps about operating point
    Pmap.P = ones(Sys.Ne,1);
    Cmap.C = CalcCFcn(Sys,x0,ones(Sys.Nv+Sys.Nev,1));
%     Cmap.C = CLookup(Sys,x0,ones(Sys.Nv+Sys.Nev,1)); % for the lookup map
    
    %% lookup map derivates
    Pmap.P_xt = zeros(Sys.Ne,1);
    Pmap.P_xh = zeros(Sys.Ne,1);
    Pmap.P_u1 = zeros(Sys.Ne,1);
    Pmap.P_u2 = zeros(Sys.Ne,1);
    Cmap.C_x  = CalcDCFcn(Sys,x0,zeros(Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev));
%     Cmap.C_x  = DCLookup(Sys,x0,zeros(Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev));
    
    %% Calculate the Linear Edge coefficients
    P = CalcP(Sys,x0,u0,Pmap); % Calculate Power Flows
    [dP_dx, dP_du] = CalcDP(Sys,x0,u0,Pmap); % Calculate Power Flow Derivatives
    
    %% Calculate the Linear Vertex coefficients
    [C,dC_dx] = CalcC(Sys,x0,Cmap); % Calculate Capacitance and Capacitance Derivatives
    
    %% Linearize the Graph Model
    idxD = (Sys.C_coeff~= 0);
    idxA = (Sys.C_coeff(1:Sys.Nv)== 0);
    Sys.ContLinGraph.A = -[-C^-1*dC_dx*kron(C,eye(Sys.Nv+Sys.Nev))^-1*kron(Sys.M_mod(idxD,:)*P,eye(Sys.Nv+Sys.Nev)) + C^-1*Sys.M_mod(idxD,:)*dP_dx; Sys.M_mod(idxA,:)*dP_dx]; 
    Sys.ContLinGraph.B = -[C^-1*Sys.M_mod(idxD,:)*dP_du; Sys.M_mod(idxA,:)*dP_du];
    Sys.ContLinGraph.W =  [C^-1*Sys.D(idxD,:); Sys.D(idxA,:)];
    

end

%%
function [P] = CalcP(Sys,x0,u0,Pmap)

P   = zeros(size(Sys.P_coeff_mod));
xt = Sys.Tails*x0;
xh = Sys.Heads*x0;
Nu = numel(fieldnames(Sys.B));
u = cell(1,Nu);
for i = 1:Nu
    u{:,i} = repmat(Sys.B.(['B',num2str(i)])*u0,3,1);
end

for i = 1:size(P,2)
    P(:,i) = Sys.PowerFlow(i).calcP(Sys.P_coeff_mod(:,i),repmat(xt,3,1),repmat(xh,3,1),u{:},repmat(Pmap.P,3,1));
end

P = sum(P,2);

end

%%
function [dP_dx,dP_du] = CalcDP(Sys,x0,u0,Pmap)

d_.dP_dxt = zeros(size(Sys.P_coeff_mod));
d_.dP_dxh = zeros(size(Sys.P_coeff_mod));

xt = Sys.Tails*x0;
xh = Sys.Heads*x0;
Nu = numel(fieldnames(Sys.B));
u  = cell(1,Nu);
dP_map = cell(1,Nu);
for i = 1:Nu
    u{:,i} = repmat(Sys.B.(['B',num2str(i)])*u0,3,1);
    dP_map{:,i} = repmat(Pmap.(['P_u' num2str(i)]),3,1);
    d_.(['dP_du' num2str(i)]) = zeros(size(Sys.P_coeff_mod));
end

for i = 1:size(d_.dP_dxt,2)
    J      = Sys.PowerFlow(i).calcJ(Sys.P_coeff_mod(:,i),repmat(xt,3,1),repmat(xh,3,1),u{:},repmat(Pmap.P,3,1),repmat(Pmap.P_xt,3,1),repmat(Pmap.P_xh,3,1),dP_map{:});
    
    d_.dP_dxt(:,i) = J(:,1);
    d_.dP_dxh(:,i) = J(:,2);
    for j = 1:Nu
        d_.(['dP_du' num2str(j)])(:,i) = J(:,2+j);
    end
end

d_.dP_dxt = sum(d_.dP_dxt,2);
d_.dP_dxh = sum(d_.dP_dxh,2);
for i = 1:Nu
    d_.(['dP_du' num2str(i)]) = sum(d_.(['dP_du' num2str(i)]),2);
end

dP_dx = diag(d_.dP_dxt)*repmat(Sys.Tails,3,1) + diag(d_.dP_dxh)*repmat(Sys.Heads,3,1);
dP_du = zeros(Sys.Ne*3,Sys.Nu);
for i = 1:Nu
    dP_du = dP_du + diag(d_.(['dP_du' num2str(i)]))*repmat(Sys.B.(['B',num2str(i)]),3,1);
end

end

%%
function [C,dC_dx] = CalcC(Sys,x0,Cmap)

idxD = (Sys.C_coeff~= 0);

g = Cmap.C;
dg_dx = Cmap.C_x;

f   = zeros(size(Sys.C_coeff));
J = zeros(size(Sys.C_coeff)); 
df_dx = zeros(Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev,Sys.Nv+Sys.Nev);

for i = 1:size(Sys.C_coeff,2)
    f(:,i) = Sys.Capacitance(i).calcC(Sys.C_coeff(:,i),x0,1); % the 1 and 0 in these lines will need to be changed
    J(:,i) = Sys.Capacitance(i).calcJ(Sys.C_coeff(:,i),x0,1,0); % the 1 and 0 in these lines will need to be changed
end

f = sum(f,2);
J = sum(J,2);

for i = 1:Sys.Nv+Sys.Nev
    df_dx(i,i,i) = J(i);
end

C = diag(g(idxD).*f(idxD));
dC_dx = unfold(dg_dx(idxD,idxD,:),1)*kron(diag(f(idxD)),eye(Sys.Nv+Sys.Nev)) +  diag(g(idxD))*unfold(df_dx(idxD,idxD,:),1);

end

%%
function X = unfold(A,n)
    d = size(A);
    X = reshape(shiftdim(permute(A,[1 3 2]),n-1), d(n), []);
end
 
%%
function [Cmap] = CLookup(Sys,x0,Cmap)
    
    val = zeros(numel(Sys.vMaps),1);
    for i = 1:numel(Sys.vMaps)
        pts = num2cell(double(Sys.vMaps(i).Pts2Map)*x0);
        val(i) = Sys.vMaps(i).LookUp(pts);
    end
    
    idx = Sys.Map2v'*double(1:(Sys.Nv+Sys.Nev))';
    Cmap(idx) = val;

end

function [DCmap] = DCLookup(Sys,x0,DCmap)
    
    val = cell(numel(Sys.vMaps),1);
    for i = 1:numel(Sys.vMaps)
        pts = num2cell(double(Sys.vMaps(i).Pts2Map)*x0)';
        pt = Sys.vMaps(i).LookUp(pts);
        val{i} = Sys.vMaps(i).LookUpDeriv(pts,pt);
        idxI = double(Sys.vMaps(i).Pts2Map)*double(1:(Sys.Nv+Sys.Nev))';
        idxO = Sys.Map2v(:,i)'*double(1:(Sys.Nv+Sys.Nev))';
        DCmap(idxO,idxO,idxI) = val{i}; 
    end
    

end

function [val] = CalcCFcn(Sys,x0,Cmap)

val = zeros(size(Sys.F_coeff));
Nf = numel(fieldnames(Sys.F));
f = cell(1,Nf);
for i = 1:Nf
    f{:,i} = Sys.F.(['F',num2str(i)])*x0;
end

for i = 1:size(val,2)
    val(:,i) = Sys.F_coeff(:,i).*Sys.LookupFunction(i).calcF(f{:});
end

val = sum(val,2);


end

function [val] = CalcDCFcn(Sys,x0,Cmap)

% val = cell(
% 
% zeros(size(Sys.F_coeff));
% Nf = numel(fieldnames(Sys.F));
% f = cell(1,Nf);
% for i = 1:Nf
%     f{:,i} = Sys.F.(['F',num2str(i)])*x0;
% end
% 
% for i = 1:size(val,2)
%     val(:,i) = Sys.F_coeff(:,i).*Sys.LookupFunction(i).calcJ(f{:});
% end
% 
% val = sum(val,2);

 %%%%%%%%%%%%

Nf = numel(fieldnames(Sys.F));
f  = cell(1,Nf);
idxI = cell(1,Nf);
idxO = cell(1,Nf);
for i = 1:Nf
    f{:,i} = Sys.F.(['F',num2str(i)])*x0;
    idx = Sys.F.(['F',num2str(i)])*double(1:(Sys.Nv+Sys.Nev))';
    idxO{:,i} = find( idx ~= 0);
    idxI{:,i} = idx(idx ~= 0); 
    d_.(['dF_dx' num2str(i)]) = zeros(size(Sys.F_coeff));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % fix this code such that the correct dC_dx function is built
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(d_.dF_dx1,2)
    J      = Sys.F_coeff(:,i).*Sys.LookupFunction(i).calcJ(f{:});
    DCmap(idxO,idxO,idxI) = J{idxO};
    for j = 1:Nf
        d_.(['dF_dx' num2str(j)])(:,i) = J(:,j);
    end
end

for i = 1:Nu
    d_.(['dF_dx' num2str(i)]) = sum(d_.(['dF_dx' num2str(i)]),2);
end

dF_dx = zeros(Sys.Ne*3,Sys.Nu);
for i = 1:Nu
    dP_du = dP_du + diag(d_.(['dP_du' num2str(i)]))*repmat(Sys.B.(['B',num2str(i)]),3,1);
end


end

