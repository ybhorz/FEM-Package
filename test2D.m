%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
clear; close all;
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
nEdge = edge.nConn; nConn = conn.nConn;
%% Real solution
syms x y;
var = [x;y];

% u = Fcn(var,x*y*(1-x/2)*(1-y)*exp(x+y)); % example 1
u = Fcn(var,exp(x+y)); % example 2,3
uFun = u.getFun;
f = -u.dif([2;0])-u.dif([0;2]); 

% f = Fcn([x;y], 1*(x>=-0.1 & x<=0.1 & y>=-0.1 & y<=0.1) ); % example 4
% f = Fcn([x;y], 1*(x>=0) ); % example 4
%% Simple tri-subdivision
% clear;
% node = [0 0 0 0.5 0.5 0.5 1 1 1; 0 0.5 1 0 0.5 1 0 0.5 1]; nNode = size(node,2);
% elem = [1 2 2 3 4 5 5 6; 4 4 5 5 7 7 8 8;2 5 3 6 5 8 6 9]; nElem = size(elem,2);
% edge = [1 2 3 6 9 8 7 4; 2 3 6 9 8 7 4 1;1 1 2 2 3 3 4 4]; nEdge = size(edge,2);
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

syms s;
EgParm = [x1,x2;y1,y2];

tfI = Tf;
tfI.refVar = s;
tfI.orgVar = [x;y];
tfI.parm = [x1,x2;y1,y2];
tfI.toRef = [(x2-x1)/2*s+(x2+x1)/2;
            (y2-y1)/2*s+(y2+y1)/2];
%% Mesh
msh = Msh;
msh.node = node;
msh.elem = elem;
msh.edge = edge;
msh.tf = tf;
msh.tfI = tfI;
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
%% P2-element
% Vetices + midpoints of connection.
P2 = FE;
P2.nNode = nNode + nConn;
P2.node = [node, (node(:,conn.node(1,:)) + node(:,conn.node(2,:)))/2]; 
P2.elem = [elem; abs(elemCn)+nNode];
P2.edge = edge;
P2.edge.node = [edge.node(1,:); nNode+(1:nEdge); edge.node(2,:)];

% P2.base(1:6) = Fcn;
% P2.base(1) = Fcn(refVar, 2*l^2 + 2*m^2 + 4*l*m - 3*l - 3*m + 1);
% P2.base(2) = Fcn(refVar, 2*l^2 - l);
% P2.base(3) = Fcn(refVar, 2*m^2 - m);
% P2.base(4) = Fcn(refVar, -4*l^2 - 4*l*m + 4*l);
% P2.base(5) = Fcn(refVar, 4*l*m);
% P2.base(6) = Fcn(refVar, -4*m^2 - 4*l*m + 4*m);
P2.base = Base(refVar,[1,l,m,l^2,l*m,m^2],[1,-3,-3,2,4,2; 0,-1,0,2,0,0; 0,0,-1,0,0,2; 0,4,0,-4,-4,0; 0,0,0,0,4,0; 0,0,4,0,-4,-4]');
P2.tf = tf;
P2.baseTf = ones(1,6);
%% Finite element space: Uh
Uh = P1;
trl = Uh; tst = Uh;
%% Varitional equation
Auv(1:2) = DLF;
Auv(1) = DLF(Fcn(var,1),1,1,[1;0],[1;0]);
Auv(2) = DLF(Fcn(var,1),1,1,[0;1],[0;1]);

Fv(1) = SLF(f,1,[0;0]);

%%%%%%%%%%%%%%%%%%% Nitsche method %%%%%%%%%%%%%%%%%%%%
% % Normal vector.
% n(1) = Fcn(var, 1/sqrt((x2-x1)^2+(y2-y1)^2)*(y2-y1), EgParm);
% n(2) = Fcn(var, - 1/sqrt((x2-x1)^2+(y2-y1)^2)*(x2-x1), EgParm);
% 
% Auv(3) = DLF(-n(1),1,1,[1;0],[0;0],"edge");
% Auv(4) = DLF(-n(2),1,1,[0;1],[0;0],"edge");
% 
% gamm = 100; 
% Auv(5) = DLF(Fcn(var,gamm/msh.h),1,1,[0;0],[0;0],"edge");
% Auv(6) = DLF(-n(1),1,1,[0;0],[1;0],"edge");
% Auv(7) = DLF(-n(2),1,1,[0;0],[0;1],"edge");
% 
% % Boundary
% uBd = u;
% Fv(2) = SLF(uBd*(gamm/msh.h),1,[0;0],"edge");
% Fv(3) = SLF(-n(1).*uBd,1,[1;0],"edge");
% Fv(4) = SLF(-n(2).*uBd,1,[0;1],"edge");

% % Robin boundary %example 3
Auv(3) = DLF(Fcn(var,1),1,1,[0;0],[0;0],"edge");

% Neuman boundary % example 2
% type 3: y = -1
% Fv(2) = SLF(Fcn(var,-exp(x-1)),1,[0;0],"edge",find(edge.type==3));

%% Assemble matrix
[Stiff,Load,nTrlNodes] = assemble(msh,trl,tst,Auv,Fv);
%% Dirichlet boundary 
% type 1: y = 1
% type 2: x = 1
% type 3: y = -1 
% type 4: x = -1

uBd1 = @(x) uFun([x(1);1]);
uBd2 = @(x) uFun([1;x(2)]);
uBd3 = @(x) uFun([x(1);-1]);
uBd4 = @(x) uFun([-1;x(2)]);

for iEdge = 1:Uh.nEdge
    for iNode = 1:Uh.edge.nNode
        I = Uh.edge.node(iNode,iEdge);
        coord = Uh.node(:,I);
        switch Uh.edge.type(iEdge)
            case 1 % y = 1
                Stiff(I,:) = 0;
                Stiff(I,I) = 1;
                Load(I) = uBd1(coord);
            case 2 % x = 1
                Stiff(I,:) = 0;
                Stiff(I,I) = 1;
                Load(I) = uBd2(coord);
%             case 3 % y = -1 % example 1, 4
%                 Stiff(I,:) = 0;
%                 Stiff(I,I) = 1;
%                 Load(I) = uBd3(coord);
            case 4 % x = -1
                Stiff(I,:) = 0;
                Stiff(I,I) = 1;
                Load(I) = uBd4(coord);
        end
    end
end
%% Solve equation
uh = Stiff\Load;
% tol = 1e-8;
% maxit = 100;
% Stiff = sparse(Stiff);
% [L,U] = ilu(Stiff,struct('type','ilutp','droptol',1e-6));
% uh = gmres(Stiff,Load,[],tol,maxit,L,U);

%% Error norm
e_h0 = eNorm(msh,Uh,uh,u,[0;0]);
e_inf = eNorm(msh,Uh,uh,u,[0;0],"norm",inf);
e_h1 = eNorm(msh,Uh,uh,u,{[1;0],[0;1]});
fprintf('h: %f e_h0: %f e_inf: %f e_h1: %f\n',msh.h,e_h0,e_inf,e_h1);
%% Plot
figure;
trisurf(elem',node(1,:)',node(2,:)',uh(1:nNode));

uVal = zeros(nNode,1);
for iEdge = 1:length(node)
    uVal(iEdge) = uFun(node(:,iEdge));
end
figure;
trisurf(elem',node(1,:)',node(2,:)',uVal);