%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
load("mesh2D_[-1,1,-1,1]\mesh2D_h0.mat");
% p: node coordinate
% e: boundary 
%   -node 1; -node 2 
%   -cumulate distance 1; -cumulate distance 2 
%   -which boundary
%   -?
%   -which domain
% t: triangle element
%   -node 1; -node 2; node 3;
%   -which domain
node = p; nNode = size(node,2);
elem = t(1:3,:); nElem = size(elem,2);
[conn,elemCn,edge] = genConn(elem,e([1,2,5],:)); 
nEdge = edge.nConn;
%% Real solution
syms x y t;
var = [x;y];

u = Fcn(var,exp(x+y+t));
f = diff(u,t,1) - diff(u,x,2)*2 - diff(u,y,2)*2;
%% transformation
syms l m x1 x2 x3 y1 y2 y3;
refVar = [l;m];

tf = Tf;
tf.refVar = [l;m];
tf.orgVar = [x;y];
tf.parm = [x1,x2,x3;y1,y2,y3];
B = [x2-x1, x3-x1; y2-y1, y3-y1];
b = [x1;y1];
tf.toRef = B * tf.refVar + b;
tf.toOrg = B \ (tf.orgVar - b);
%% Mesh
msh = Msh;
msh.node = node;
msh.elem = elem;
msh.edge = edge;
msh.conn = conn;
msh.tf = tf;
msh.gInt = GInt("D2P3");
msh.gIntI = GInt("D1P4");
msh.check;
%% p1-element
P1 = FE;
P1.nNode = nNode;
P1.node = node;
P1.elem = elem;
P1.edge = edge;
P1.base = Base(refVar,[1,l,m],[1,-1,-1;0,1,0;0,0,1]');
P1.tf = tf;
P1.baseTf = [1,1,1];
%% Assmeble matrix
Uh = P1; iu = 1;
d0u = [0;0]; dxu = [1;0]; dyu = [0;1];

Muv = DLF(Fcn(var,1),iu,iu,d0u,d0u);
Auv(1) = DLF(Fcn(var,2),iu,iu,dxu,dxu);
Auv(2) = DLF(Fcn(var,2),iu,iu,dyu,dyu);
Fv = SLF(f,iu,d0u);

[Stiff,~] = assemble(msh,Uh,Uh,Auv,SLF.empty,"Integral",GInt("D2P1"));
[Mass,~] = assemble(msh,Uh,Uh,Muv,SLF.empty,"Integral",GInt("D2P2"));
[~,Load0] = assemble(msh,FE.empty,Uh,DLF.empty,Fv.subs(t,0));
%% Time-dependent FEM solver
dltT = msh.h^2; %0.05;
T = 0:dltT:1;
nT = length(T);

% Initial condition
u0Fun = u.subs(t,0).getFun;
uh0 = zeros(Uh.nNode,1);
for iNode = 1:nNode
    uh0(iNode) = u0Fun(Uh.node(:,iNode));
end

% Algebra equation
theta = 1; % 0 - forward Eular scheme
             % 1 - backward Eular scheme
             % 1/2 - Crank-Nicolson scheme
A = Mass/dltT + theta*Stiff;
C = Mass/dltT - (1-theta)*Stiff;

for i = 2:nT
    [~,Load] = assemble(msh,FE.empty,Uh,DLF.empty,Fv.subs(t,T(i)));

    B = theta*Load + (1-theta)*Load0 + C*uh0;

    % Dirichlet boundary.
    uBd = u.subs(t,T(i)).getFun;
    
    for iEdge = 1:Uh.nEdge
        for iNode = 1:Uh.edge.nNode
            I = Uh.edge.node(iNode,iEdge);
            coord = Uh.node(:,I);
            switch Uh.edge.type(iEdge)
                case 1 % y = 1
                    A(I,:) = 0;
                    A(I,I) = 1;
                    B(I) = uBd(coord);
                case 2 % x = 1
                    A(I,:) = 0;
                    A(I,I) = 1;
                    B(I) = uBd(coord);
                case 3 % y = -1 
                    A(I,:) = 0;
                    A(I,I) = 1;
                    B(I) = uBd(coord);
                case 4 % x = -1
                    A(I,:) = 0;
                    A(I,I) = 1;
                    B(I) = uBd(coord);
            end
        end
    end

    uh = A\B;

    e_h0 = eNorm(msh,Uh,uh,u.subs(t,T(i)),d0u,"norm",2);
    e_inf = eNorm(msh,Uh,uh,u.subs(t,T(i)),d0u,"norm",inf);
    e_h1 = eNorm(msh,Uh,uh,u.subs(t,T(i)),{dxu,dyu},"norm",2);
    
    uh0 = uh;
    Load0 = Load;
end

fprintf('h: %f\nt: %f\ne_h0: %f\ne_inf: %f\ne_h1: %f\n',msh.h,dltT,e_h0,e_inf,e_h1);
