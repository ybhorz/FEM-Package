classdef FE
    % FE : Finite element space.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com
    
    properties
        nNode; % The number of FE nodes.
        node (:,:); % Node variables (node coordinates for lagrange element, while have no sense for other element.)
        elem (:,:); % Indices of FE nodes in mesh element. (By coloumn)
                    % i.e. indices of the base functions in mesh element.
                    % Note: FE nodes are arranged as the same order of base functinos.
                    %! Note: negative index means the corresponding base function take nagetive value.
        edge Conn;  % Edge elements and Edge nodes.
                    % Note: nodes are arranged counter-clockwise in the element 
                    %! Note: 2 dimension
        base Base; % Basis functions in reference element. (By coloumn for vector-FE)
        tf Tf; % Transformation between basis funciton in mesh element and basis function in reference element.
        baseTf; % Correspondence between base and transformation.
    end
    properties (Dependent)
        dim; % Dimension.
        nElem; % The number of elements.
        nEdge; % The number of edge parts.
        nBase; % The number of base functions in one element.
        nTf; % The number of transformation.
    end
    methods
        function disp(fE)
            for iFE = 1:length(fE)
                divLn('#',64,[inputname(1),'(',num2str(iFE),')']);
                fprintf('node: %d\n',fE(iFE).nNode);
                fprintf('elem: [%d %d]\n', fE(iFE).nBase,fE(iFE).nElem);
                if ~isempty(fE(iFE).edge)
                    divLn('=',64,'edge');
                    disp(fE(iFE).edge);
                end
                disp(fE.base);
                for iTf = 1:fE(iFE).nTf
                    divLn('=',64,['tf(',num2str(iTf),')']);
                    disp(fE(iFE).tf(iTf));
                end
            end
            divLn('#',64);
        end
        function fE = set.tf(fE,tf)
            if ~isvector(tf)
                error('FE : "tf" must be a vector.');
            end
            fE.tf = tf;
        end
        function fE = set.baseTf(fE,baseTf)
            if isvector(baseTf)
                fE.baseTf = baseTf;
            else
                error('FE : "baseTf" must be a vector.');
            end
        end
        function dim = get.dim(fE)
            dim = fE.base.nVar;
            if any(dim ~= [fE.tf.nVar])
                error('FE: the dimension is inconsistent - "base", "tf".');
            end
        end
        function nElem = get.nElem(fE)
            nElem = size(fE.elem,2);
        end
        function nEdge = get.nEdge(fE)
            nEdge = fE.edge.nConn;
        end
        function nBase = get.nBase(fE)
            nBase = fE.base.nBase;
        end
        function nTf = get.nTf(fE)
            nTf = length(fE.tf);
        end
        function fE = plus(fE1,fE2)
            % Overload + operator
            fE = FE;
            fE.nNode = fE1.nNode + fE2.nNode;
            if ~isempty(fE1.node)
                fE.node = fE1.node;
                if ~isempty(fE2.node)
                    fE.node = [fE.node, fE2.node];
                else
                    fE.node = [fE.node, nan(fE2.dim,fE2.nNode)];
                end
            elseif ~isempty(fE2.node)
                fE.node = [nan(fE1.dim,fE1.nNode), fE2.node];
            end
            fE.elem = [fE1.elem; sign(fE2.elem).*(fE1.nNode+abs(fE2.elem))];
            fE.base = fE1.base + fE2.base;
            fE.tf = [fE1.tf, fE2.tf];
            fE.baseTf = [fE1.baseTf, fE1.nTf+fE2.baseTf];
            if ~isempty(fE1.edge)
                fE.edge = fE1.edge;
                if ~isempty(fE2.edge)
                    fE.edge.node = [fE.edge.node; fE2.nNode+fE2.edge.node];
                end
            elseif ~isempty(fE2.edge)
                fE.edge = fE2.edge;
            end
        end
        function fE = power(fE,r)
            % Overload .^ operator.
            if r ~= 2
                error('FE.power: "r" must be "2".');
            end
            nNode0 = fE.nNode;
            fE.nNode = 2*nNode0;
            if ~isempty(fE.node)
                fE.node= [fE.node,fE.node];
            end
            fE.elem = [fE.elem; sign(fE.elem).*(nNode0+abs(fE.elem))];
            fE.base = fE.base.^2;
            fE.tf = fE.tf;
            fE.baseTf = [fE.baseTf,fE.baseTf];
            if ~isempty(fE.edge)
                fE.edge.node = [fE.edge.node; nNode0+fE.edge.node];
            end
        end
        function check(fE)
            if ~isempty(fE.node)
%                 if fE.dim ~= size(fE.node,1)
%                     error('FE: the size of "node" is inconsistent with "dim".')
%                 end
                if fE.nNode ~= size(fE.node,2)
                    error('FE: the length of "node" is inconsistent with "nNode".');
                end
            end
            if ~isempty(fE.edge)
                fE.edge.check;
            end
            fE.base.check;
            if any(diff([size(fE.elem,1), fE.base.nBase, length(fE.baseTf)]))
                error('FE: the length is inconsistent - "elem", "base", "baseTf".');
            end
            for iTf = 1:fE.nTf
                fE.tf(iTf).check;
            end
        end
    end
end

