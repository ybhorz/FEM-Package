classdef Msh
    % Msh : Mesh information.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com
    
    properties
        node (:,:); % Coordinates of mesh nodes (By column).
        elem (:,:); % Vertices (node index) of mesh elements (By column). 
                    % Note: Nodes are arranged in a certain order.(counter-clockwise)
        edge Conn;  % Edge elements and Edge nodes.
                    % Note: nodes are arranged counter-clockwise in the element 
                    %! Note: 2 dimension
        conn Conn;  % Connection nodes, connected elements.
                    % conn.elem > 0: nodes are arranged counter-clockwise in this element.
                    % conn.elem < 0: nodes are arranged clockwise in this element.
                    % conn.elem = 0: no connected element.
                    % Note: edge may be included in conn.
        %refElem Rgn = Rgn; % Region of reference element. 
        tf Tf; % Transformation between mesh element and reference element.
        tfI Tf; % Transformation from edge to reference interval [-1,1]. %! 2 dimension
        gInt GInt; % Guass numerical integration in reference element.
        gIntI GInt; % Guass numerical integration in reference interval [-1,1]. %! 2 dimension
    end
    properties (Dependent)
        dim; % Dimension.
        nNode; % The number of all mesh nodes.
        nElem; % The number of mesh elements.
        nEdge; % The number of edge parts.
        nConn; % The number of connection.
        h; % The size of mesh.
    end
    methods
        function disp(msh)
            divLn('#',64,inputname(1));
            fprintf('node: [%d  %d]\n',msh.dim,msh.nNode);
            fprintf('elem: [%d  %d]\n',size(msh.elem,1),msh.nElem);
            if ~isempty (msh.edge)
                divLn('=',64,'edge');
                disp(msh.edge);
            end
            if ~isempty(msh.conn)
                divLn('=',64,'conn');
                disp(msh.conn);
            end
%             divLn('=',64,'refElem');
%             disp(msh.refElem);
            divLn('=',64,'tf');
            disp(msh.tf);
            if ~isempty(msh.tfI)
                divLn('=',64,'tfI');
                disp(msh.tfI);
            end
            divLn('=',64,'gInt');
            disp(msh.gInt);
            if ~isempty(msh.gIntI)
                divLn('=',64,'gIntI');
                disp(msh.gIntI);
            end
            divLn('#',64);
        end
        function dim = get.dim(msh)
            dim = size(msh.node,1);
        end
        function nNode = get.nNode(msh)
%             nNode(1) = size(msh.elem,1);
            nNode = size(msh.node,2);
        end
        function nElem = get.nElem(msh)
            nElem = size(msh.elem,2);
        end
        function nEdge = get.nEdge(msh)
            nEdge = msh.edge.nConn;
        end
        function nConn = get.nConn(msh)
            nConn = msh.conn.nConn;
        end
        function h = get.h(msh)
            h = 0;
            Jdet = msh.tf.getJdet.getFun;
            for iElem = 1:msh.nElem
                parm = msh.node(:,msh.elem(:,iElem));
                h = max(h, Jdet([1/3;1/3],parm)^(1/msh.dim));
            end
        end
        function check(msh)
            if ~isempty(msh.edge)
                msh.edge.check;
            end
            if ~isempty(msh.conn)
                msh.conn.check;
            end
            msh.tf.check;
%             if ~isempty(msh.tfI)
%                 msh.tfI.check;
%             end
            msh.gInt.check;
            if ~isempty(msh.gIntI)
                msh.gIntI.check;
            end
            if any( diff([msh.dim, msh.tf.nVar, msh.gInt.dim]) ~= 0 )
                error('Msh: the dimension is inconsistent - "node", "tf", "gInt".');
            end
            if msh.tf.nParm(2) ~= size(msh.elem,1)
                error('Msh: the numer of variables is inconsistent - "tf.parm", "elem".');
            end
        end
    end
end