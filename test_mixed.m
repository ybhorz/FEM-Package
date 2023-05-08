%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
clear; close all;
load("mesh2D_[-1,1,-1,1]\mesh2D_h1.mat");
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
[conn,elemCn,edge] = genConn(elem,e([1,2,5],:)); nEdge = edge.nConn;
%% Real solution
syms x y;
var = [x;y];

p = Fcn(var,(1-x)*(1+x)*(1-y)*(1+y));
u = -[p.dif([1;0]); p.dif([0;1])];
f = -p.dif([2;0]) - p.dif([0;2]);
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
msh.tf = tf;
msh.gInt = GInt('D2P3');
msh.gIntI = GInt('D1P4');
msh.check;
%% Finite element space: RT0
RT0 = FE;
RT0.nNode = conn.nConn;
RT0.node = (node(:,conn.node(1,:))+node(:,conn.node(2,:)))/2;
RT0.elem = elemCn;
RT0.edge = edge;
RT0.edge.node = 1:nEdge;
RT0.base = Base(refVar,[1,0,l;0,1,m],[0,-1,1;0,0,1;-1,0,1]');
RTtf = tf;
RTtf.refCoef = RTtf.getJdet.fun * inv(B);
RTtf.orgCoef = 1/RTtf.getJdet.fun * B;
RT0.tf = RTtf; RT0.baseTf = ones(1,RT0.nBase);
%% Finite element space: P0
P0 = FE;
P0.nNode = nElem;
P0.elem = 1:nElem;
org = msh.tf.getOrg;
P0.node = zeros([2,nElem]);
for iElem = 1:nElem
    parm = node(:,elem(:,iElem));
    P0.node(:,iElem) = org([1/3;1/3],parm);
end
P0.base = Base(refVar,1,1);
P0.tf = tf; P0.baseTf = 1;
%% Varitional equation
Uh = RT0; Ph = P0;
trls = [Uh,Ph]; tsts = trls;
iu = 1; ip = 2;
d0ux = [0,nan;0,nan];
d0uy = [nan,0;nan,0];
divu = [1,0;0,1];
d0p = [0;0];

Auv(1) = DLF(Fcn(var,1),iu,iu,d0ux,d0ux);
Auv(2) = DLF(Fcn(var,1),iu,iu,d0uy,d0uy);

Buq = DLF(Fcn(var,1),iu,ip,divu,d0p);
Bpv = DLF(Fcn(var,-1),ip,iu,d0p,divu);

Fq = SLF(f,ip,d0p);
%% Solve equation
[Stiff,Load,nTrlNodes] = assemble(msh,trls,tsts,[Auv,Bpv,Buq],Fq);

xh = Stiff\Load;

uh = xh(nTrlNodes(iu)+1:nTrlNodes(iu+1));
ph = xh(nTrlNodes(ip)+1:nTrlNodes(ip+1));
%% Error norm
ep_h0 = eNorm(msh,Ph,ph,p,d0p);
eu_h0 = eNorm(msh,Uh,uh,u,{d0ux,d0uy});
eu_div = eNorm(msh,Uh,uh,u,divu);
fprintf('h: %f e(p)_h0: %f e(u)_h0: %f e(u)_div: %f\n',msh.h,ep_h0,eu_h0,eu_div);
%% Plot
figure; view(3);
for iElem = 1:nElem
    patch(node(1,elem(:,iElem))',node(2,elem(:,iElem))',repmat(ph(iElem),[3,1]),ph(iElem));
end
title('Approximation: p');

figure; view(3);
pFun = p.tfRef(tf).getFun;
for iElem = 1:nElem
    parm = node(:,elem(:,iElem));
    pVal = pFun([1/3;1/3],parm);
    patch(node(1,elem(:,iElem))',node(2,elem(:,iElem))',repmat(pVal,[3,1]),pVal);
end
title('Real solution: p');