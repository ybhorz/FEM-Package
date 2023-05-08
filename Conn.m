classdef Conn
    % Conn: Connection
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com

    properties
        node (:,:); % Vertices of connection (By coloumn);
                    % Note: Nodes are arranged in a certain order.
        elem (:,:); % Connected element (By coloumn);
                    % elem > 0: nodes are arranged counter-clockwise in this element.
                    % elem < 0: nodes are arranged clockwise in this element.
                    % elem = 0: no connected element.
%         locl (:,:); % 
        type; % Type of connection. For example: 1\2\3\4 for edge, 0 for interior connection.
    end
    properties(Dependent)
        nConn; % Total number of connection.
        nNode; % Number of nodes in one connection.
        nElem; % Number of connected element in one connection.
    end
    methods
        function disp(conn)
            fprintf('node: [%d  %d]\n',conn.nNode,conn.nConn);
            fprintf('elem: [%d  %d]\n',conn.nElem,size(conn.elem,2));
        end
        function conn = set.type(conn,type)
            if isvector(type)
                conn.type = type;
            else
                error('Conn: "type" must be vector.');
            end
        end
        function nConn = get.nConn(conn)
            nConn = size(conn.node,2);
        end
        function nNode = get.nNode(conn)
            nNode = size(conn.node,1);
        end
        function nElem = get.nElem(conn)
            nElem = size(conn.elem,1);
        end
        function conn = part(conn,subs)
        % Take a part of conn;
            conn.node = conn.node(:,subs);
            conn.elem = conn.elem(:,subs);
            conn.type = conn.type(:,subs);
        end
        function check(conn)
            if ~isempty(conn.elem) && size(conn.elem,2)~=conn.nConn
                error('Conn: the number of variables is inconsistent - "node", "elem".');
            end
            if ~isempty(conn.type) && size(conn.type,2)~=conn.nConn
                error('Conn: the number of variables is inconsistent - "node", "type".');
            end
        end
    end
end