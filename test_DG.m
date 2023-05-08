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
[conn,elemCn,edge] = genConn(elem,e([1,2,5],:)); 
nEdge = edge.nConn; nConn = conn.nConn;
%% Real solution
syms x y;
var = [x;y];

% u = Fcn(var,(1-x)*(1+x)*(1-y)*(1+y));
u = Fcn(var,x*y*(1-x/2)*(1-y)*exp(x+y));
uFun = u.getFun;
f = - u.dif([2;0]) - u.dif([0;2]);

% f = Fcn([x;y], 1*(x>=-0.1 & x<=0.1 & y>=-0.1 & y<=0.1) );
% f = Fcn([x;y], 1*(x>=0) );
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
msh.conn = conn.part(nEdge+1:nConn); %! Note
msh.tf = tf;
msh.tfI = tfI;
msh.gInt = GInt("D2P3");
msh.gIntI = GInt("D1P4");
msh.check;
%% Finite element space: Uh
Uh = FE;
Uh.nNode = 3*nElem;
Uh.elem = reshape(1:Uh.nNode,[3,nElem]);
Uh.node = zeros(2,Uh.nNode);
for iElem = 1:Uh.nElem
    for iNode = 1:3
        Uh.node(:,Uh.elem(iNode,iElem)) = node(:,elem(iNode,iElem));
    end
end

Uh.edge = edge;
for iEdge = 1:nEdge
    iElem = edge.elem(iEdge);
    switch find(elemCn(:,iElem)==iEdge)
        case 1
            Uh.edge.node(:,iEdge) = Uh.elem([1,2],iElem);
        case 2
            Uh.edge.node(:,iEdge) = Uh.elem([2,3],iElem);
        case 3
            Uh.edge.node(:,iEdge) = Uh.elem([3,1],iElem);
    end
end

%!Note: Elements with one edge vertex are also included.
for iEdge = 1:nEdge
    for iElem = setdiff(1:nElem,Uh.edge.elem)
        for iNode = 1:3
            if ismember(elem(iNode,iElem),edge.node(:,iEdge))
                Uh.edge.node = [Uh.edge.node, [Uh.elem(iNode,iElem);0]];
                Uh.edge.elem = [Uh.edge.elem, iElem];
                Uh.edge.type = [Uh.edge.type, edge.type(iEdge)];
                break;
            end
        end
    end
end

% P1(1:3) = Fcn;
% P1(1) = Fcn(refVar,1-l-m);
% P1(2) = Fcn(refVar,l);
% P1(3) = Fcn(refVar,m);
% 
% Uh.base = P1;
Uh.base = Base(refVar,[1-l-m,l,m],eye(3));

Uh.tf = repmat(tf,[1,3]);
Uh.baseTf = 1:3;
Uh.check;

trl = Uh; tst = Uh;
%% Test: Uh
% figure;
% orgI = msh.tf.getOrg;
% alph = 0.5;
% for iElem = 1:Uh.nElem
%     vert = Uh.node(:,Uh.elem(:,iElem));
%     cent = orgI([1/3;1/3],vert); 
%     patch('Faces',[1,2,3],'Vertices',vert',...
%     'EdgeColor','black','FaceColor','none','LineWidth',0.1);
%     for iNode = 1:3
%         x = alph*vert(1,iNode) + (1-alph)*cent(1); 
%         y = alph*vert(2,iNode) + (1-alph)*cent(2);
%         text(x,y,num2str(Uh.elem(iNode,iElem)))
%     end
% end
% return;
%% Varitional equation
Auv(1:2) = DLF;
Auv(1) = DLF(Fcn(var,1),1,1,[1;0],[1;0]);
Auv(2) = DLF(Fcn(var,1),1,1,[0;1],[0;1]);

% Normal vector.
n(1) = Fcn;
n(1).var = [x;y];
n(1).parm = [x1,x2;y1,y2];
n(1).fun = 1/sqrt((x2-x1)^2+(y2-y1)^2) * (y2-y1);
n(2) = Fcn;
n(2).var = [x;y];
n(2).parm = [x1,x2;y1,y2];
n(2).fun = - 1/sqrt((x2-x1)^2+(y2-y1)^2) * (x2-x1);

Buv(1:8) = DLF;
iBuv = 0;
for trlSgn = [1,-1]
    for tstSgn = [1,-1]
        for dirc = 1:2
            iBuv = iBuv + 1;
            Buv(iBuv).coef = n(dirc) * (-1/2*tstSgn);
            Buv(iBuv).iTrl = trlSgn;
            Buv(iBuv).iTst = tstSgn;
            if dirc == 1
                Buv(iBuv).trlOrd = [1;0];
            else
                Buv(iBuv).trlOrd = [0;1];
            end
            Buv(iBuv).tstOrd = [0;0];
            Buv(iBuv).type = "conn";
        end
    end
end

gam = -1;
Bvu(1:8) = DLF;
iBvu = 0;
for trlSgn = [1,-1]
    for tstSgn = [1,-1]
        for dirc = 1:2
            iBvu = iBvu + 1;
            Bvu(iBvu).coef = n(dirc) * (-1/2*trlSgn*gam);
            Bvu(iBvu).iTrl = trlSgn;
            Bvu(iBvu).iTst = tstSgn;
            Bvu(iBvu).trlOrd = [0;0];
            if dirc == 1
                Bvu(iBvu).tstOrd = [1;0];
            else
                Bvu(iBvu).tstOrd = [0;1];
            end
            Bvu(iBvu).type = "conn";
        end
    end
end

sgm = 1;
hI = Fcn;
hI.var = [x;y];
hI.parm = [x1,x2;y1,y2];
hI.fun = sqrt((x1-x2)^2+(y1-y2)^2);

Cuv(1:4) = DLF;
iCuv = 0;
for trlSgn = [1,-1]
    for tstSgn = [1,-1]
        iCuv = iCuv + 1;
        Cuv(iCuv).coef = hI\sgm * trlSgn * tstSgn;
        Cuv(iCuv).iTrl = trlSgn;
        Cuv(iCuv).iTst = tstSgn;
        Cuv(iCuv).trlOrd = [0;0];
        Cuv(iCuv).tstOrd = [0;0];
        Cuv(iCuv).type = "conn";
    end
end

Fv = SLF(f,1,[0;0]);
%% Assemble matrix
[Stiff,Load,nTrlNodes] = assemble(msh,trl,tst,[Auv,Buv,Bvu,Cuv],Fv);
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
        if I == 0
            continue;
        end
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
            case 3 % y = -1
                Stiff(I,:) = 0;
                Stiff(I,I) = 1;
                Load(I) = uBd3(coord);
            case 4 % x = -1
                Stiff(I,:) = 0;
                Stiff(I,I) = 1;
                Load(I) = uBd4(coord);
        end
    end
end

uh = Stiff\Load;
%% Error Norm
e_h0  = eNorm(msh,Uh,uh,u,[0;0]);
e_h1 = eNorm(msh,Uh,uh,u,{[1;0],[0;1]});
fprintf('h: %f e_h0: %f e_h1: %f\n',msh.h,e_h0,e_h1);
%% Plot
figure; view(3);
for iElem = 1:nElem
    patch(node(1,elem(:,iElem))',node(2,elem(:,iElem))',uh(Uh.elem(:,iElem)),uh(Uh.elem(:,iElem)));
end
