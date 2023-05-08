function [node,elem,edge,conn,elemCn] = genTri(domn,nSub)
% genTri: generate triangulation.
% Author: Bohan Yang
% Email: ybh.orz@gmail.com

% Input arguments.
% domn: xmin, xmax, ymin, ymax
% nSub: number of subdivision in each direction.
%   nx = ny = nSub
%   nx = nSub(1), ny = nSub(2)

a1 = domn(1); a2 = domn(2); b1 = domn(3); b2 = domn(4);
if isscalar(nSub)
    nx = nSub;
    ny = nSub;
else
    nx = nSub(1);
    ny = nSub(2);
end

xI = linspace(a1,a2,nx+1);
yJ = linspace(b1,b2,ny+1);
nNode = (nx+1)*(ny+1);
node = zeros(2,nNode);
for iNode = 1:nNode
    [i,j] = ind2sub([nx+1,ny+1],iNode);
    node(:,iNode) = [xI(i);yJ(j)];
end

nRect = nx*ny; nElem = 2*nRect;
elem = zeros(3,nElem);
for iRect = 1:nRect
    [i,j] = ind2sub([nx,ny],iRect);
    iElem1 = 2*iRect - 1;
    elem(:,iElem1) = [  sub2ind([nx+1,ny+1],i,j);
                        sub2ind([nx+1,ny+1],i+1,j);
                        sub2ind([nx+1,ny+1],i,j+1)];
    iElem2 = 2*iRect;
    elem(:,iElem2) = [  sub2ind([nx+1,ny+1],i+1,j);
                        sub2ind([nx+1,ny+1],i+1,j+1);
                        sub2ind([nx+1,ny+1],i,j+1)];
end

elemCn = zeros(3,nElem);

nEdge = 2*nx + 2*ny;
edge = Conn;
edge.node = zeros(2,nEdge);
edge.elem = zeros(1,nEdge);
edge.type = zeros(1,nEdge);
for i = 1:nx
    iEdge = i;
    edge.node(:,iEdge) = [ sub2ind([nx+1,ny+1],i+1,ny+1);
                        sub2ind([nx+1,ny+1],i,ny+1)];
    edge.type(iEdge) = 1;
    iRect = sub2ind([nx,ny],i,ny); iElem = 2*iRect;
    edge.elem(iEdge) = iElem;
    elemCn(2,iElem) = iEdge;

    iEdge = nx+ny+i;
    edge.node(:,iEdge) = [ sub2ind([nx+1,ny+1],i,1);
                        sub2ind([nx+1,ny+1],i+1,1)];
    edge.type(iEdge) = 3;
    iRect = sub2ind([nx,ny],i,1); iElem = 2*iRect-1;
    edge.elem(iEdge) = iElem;
    elemCn(1,iElem) = iEdge;
end
for j = 1:ny
    iEdge = nx+j;
    edge.node(:,iEdge) = [ sub2ind([nx+1,ny+1],nx+1,j);
                        sub2ind([nx+1,ny+1],nx+1,j+1)];
    edge.type(iEdge) = 2;
    iRect = sub2ind([nx,ny],nx,j); iElem = 2*iRect;
    edge.elem(iEdge) = iElem;
    elemCn(1,iElem) = iEdge;

    iEdge = 2*nx+ny+j;
    edge.node(:,iEdge) = [ sub2ind([nx+1,ny+1],1,j+1);
                        sub2ind([nx+1,ny+1],1,j)];
    edge.type(iEdge) = 4;
    iRect = sub2ind([nx,ny],1,j); iElem = 2*iRect-1;
    edge.elem(iEdge) = iElem;
    elemCn(3,iElem) = iEdge;
end

nConn = (nx+1)*ny + nx*(ny+1) + nx*ny;
conn = Conn;
conn.node = zeros(2,nConn); conn.node(:,1:nEdge) = edge.node;
conn.elem = zeros(2,nConn); conn.elem(1,1:nEdge) = edge.elem;
conn.type = zeros(1,nConn); conn.type(1:nEdge) = edge.type;

iConn = nEdge;
for i = 1:nx-1
    for j = 1:ny
        iConn = iConn + 1;
        conn.node(:,iConn) = [sub2ind([nx+1,ny+1],i+1,j); sub2ind([nx+1,ny+1],i+1,j+1)];
        iRect1 = sub2ind([nx,ny],i,j); iRect2 = sub2ind([nx,ny],i+1,j);
        iElem1 = 2*iRect1; iElem2 = 2*iRect2-1;
        conn.elem(:,iConn) = [iElem1; -iElem2];
        elemCn(1,iElem1) = iConn;
        elemCn(3,iElem2) = -iConn;
    end
end
for i = 1:nx
    for j = 1:ny-1
        iConn = iConn + 1;
        conn.node(:,iConn) = [sub2ind([nx+1,ny+1],i+1,j+1); sub2ind([nx+1,ny+1],i,j+1)];
        iRect1 = sub2ind([nx,ny],i,j); iRect2 = sub2ind([nx,ny],i,j+1);
        iElem1 = 2*iRect1; iElem2 = 2*iRect2-1;
        conn.elem(:,iConn) = [iElem1; -iElem2];
        elemCn(2,iElem1) = iConn;
        elemCn(1,iElem2) = -iConn;
    end
end
for j = 1:ny
    for i = 1:nx
        iConn = iConn + 1;
        conn.node(:,iConn) = [sub2ind([nx+1,ny+1],i,j+1),sub2ind([nx+1,ny+1],i+1,j)];
        iRect1 = sub2ind([nx,ny],i,j); iRect2 = sub2ind([nx,ny],i,j);
        iElem1 =2*iRect1; iElem2 = 2*iRect2-1;
        conn.elem(:,iConn) = [iElem1; -iElem2];
        elemCn(3,iElem1) = iConn;
        elemCn(2,iElem2) = -iConn;
    end
end

end