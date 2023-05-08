function [conn,elemCn,edge] = genConn(elem,edge)
% genConn: generate connection and element-connection information.
%   !Note: 2 dimension.
% Author: Bohan Yang
% Email: ybh.orz@gmail.com

% Input arguments.
% elem: vertices of mesh elements (By column).
% edge: vertices of each part of edge (By column).
%       -node 1
%       -node 2
%       -type (~=0): for example 1,2,3,4 (edge)

% Output arguments.
% conn (Conn): connection.
%       Note: First nEdge connections are edges.
%       -node 1
%       -node 2
%       -connected element 1 (+)
%       -connected element 2 (-)
%       -type: ~=0 (edge) \ 0 (connection)
% elemCn: element-connection (By column)
%       -connection 1 (+\-)
%       -connection 2 (+\-)
%       -connection 3 (+\-)
%       Note: node 1 -> conn 1 -> node 2 -> conn2 -> node 3 -> conn 3 -> node 1
% edge (Conn): edge connection.
%       -node 1
%       -node 2
%       -connected element
%       -type

nElem = size(elem,2);
nEdge = size(edge,2);
% First nEdge connections are edges.
iConn = nEdge;
conn = Conn;
conn.node = zeros(2,nEdge);
conn.elem = zeros(2,nEdge);
conn.type = zeros(1,nEdge);
elemCn = zeros(size(3,nElem));

% Find element connected to edge.
for iElem = 1:nElem
    for jEdge = 1:nEdge
        inter = elem(:,iElem)' == edge(1:2,jEdge);
        if sum(inter,"all") == 2
            edgeI = sum(inter);
            if all(edgeI==[1,1,0])
                elemCn(1,iElem) = jEdge;
                connI = elem(1:2,iElem);
            elseif all(edgeI==[0,1,1])
                elemCn(2,iElem) = jEdge;
                connI = elem(2:3,iElem);
            elseif all(edgeI==[1,0,1])
                elemCn(3,iElem) = jEdge;
                connI = elem([3,1],iElem);
            end
            conn.node(:,jEdge) = connI;
            conn.elem(:,jEdge) = [iElem;0];
            conn.type(jEdge) =  edge(3,jEdge);
        end
    end
end
edge = conn;
edge.elem = edge.elem(1,:);

% Find connected elements.
for iElem = 1:nElem
    for jElem = iElem+1:nElem
        inter = elem(:,iElem)' == elem(:,jElem);
        if sum(inter,"all") == 2
            iConn = iConn + 1;
            edgeI = sum(inter);
            if all(edgeI==[1,1,0])
                elemCn(1,iElem) = iConn;
                connI = elem(1:2,iElem);
            elseif all(edgeI==[0,1,1])
                elemCn(2,iElem) = iConn;
                connI = elem(2:3,iElem);
            elseif all(edgeI==[1,0,1])
                elemCn(3,iElem) = iConn;
                connI = elem([3,1],iElem);
            end
            edgeJ = sum(inter,2);
            if all(edgeJ==[1,1,0]')
                elemCn(1,jElem) = -iConn;
            elseif all(edgeJ==[0,1,1]')
                elemCn(2,jElem) = -iConn;
            elseif all(edgeJ==[1,0,1]')
                elemCn(3,jElem) = -iConn;
            end
           
            conn.node = [conn.node, connI];
            conn.elem = [conn.elem, [iElem;-jElem]];
            conn.type = [conn.type, 0];
        end
    end
end