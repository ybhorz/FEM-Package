function [Stiff,Load,nTrlNodes] = assemble(msh,trls,tsts,Auvs,Fvs,options)
% assemble: Assemble the total stiffness matrix and load vector.
% Author: Bohan Yang
% Email: ybh.orz@gmail.com

% Input arguments.
arguments
    msh Msh; % Mesh information.
    trls FE; % Trail function spaces.
    tsts FE; % Test function spaces.
    Auvs DLF; % Double linear functions (used to assemble stiffness matrix).
    Fvs SLF; % Single linear functions (used to assemble load vector).
    options.MatrixType {mustBeMember(options.MatrixType,{'full','sparse'})} = 'full';
    options.Integral GInt = GInt.empty; % Gauss Inegral
    options.IntvlInt GInt = GInt.empty; % Interval Integral
end

% Checking and local argument.
if ~isvector(trls) && ~isempty(trls)
    error('assemble: "trls" must be a vector.');
else
    nTrl = length(trls); % The number of trial function spaces.
end
if ~isvector(tsts)
    error('assemble: "tsts" must be a vector.');
else
    nTst = length(tsts); % The number of test function spaces.
end
if isvector(Auvs) || isempty(Auvs)
    nAuv = length(Auvs); % The number of double linear functions.
else
    error('assemble: "Auvs" must be a vector.');
end
if isvector(Fvs) || isempty(Fvs)
    nFv = length(Fvs); % The number of single linear functions.
else
    error('assemble: "Fvs" must be a vector.');
end

msh.check;
dim = msh.dim; % Dimension of the problem.
for iTrl = 1:nTrl
   trls(iTrl).check;
   if dim ~= trls(iTrl).dim
      error('assemeble: the dimension is inconsistent - "msh", "trls".'); 
   end
end
for iTst = 1:nTst
   tsts(iTst).check;
   if dim ~= tsts(iTst).dim 
        error('assemble: the dimension is inconsistent - "msh", "tsts".');
   end
end
for iAuv = 1:nAuv
    Auvs(iAuv).check;
    if dim ~= Auvs(iAuv).dim
        error('assembel: the dimension is inconsistent - "msh", "Auvs".');
    end
end
for iFv = 1:nFv
    Fvs(iFv).check;
    if dim ~= Fvs(iFv).dim
        error('assemble: the dimension is inconsistent - "msh", "Fvs".');
    end
end

nElem = msh.nElem; % The number of element. 
for iTrl = 1:nTrl
   if nElem ~= trls(iTrl).nElem
      error('assemble: the "nElem" is inconsistent - "msh", "trls".'); %? May be elimited.
   end
end
for iTst = 1:nTst
   if nElem ~= tsts(iTst).nElem
      error('assemble: the "nElem" is inconsistent - "msh", "tsts".'); %? May be elimited.
   end
end

nTrlNodes = zeros(nTrl+1,1); % Number of node variables in trial function spaces. (Cumulative numbers)
for iTrl = 1:nTrl
    nTrlNodes(iTrl+1) = nTrlNodes(iTrl) + trls(iTrl).nNode;
end

nTstNodes = zeros(nTst+1,1); % Number of node variables in test function spaces. (Cumulative numbers)
for iTst = 1:nTst
    nTstNodes(iTst+1) = nTstNodes(iTst) + tsts(iTst).nNode;
end

if isempty(options.Integral)
    gInt = msh.gInt;
else
    gInt = options.Integral;
end

if isempty(options.IntvlInt)
    gIntI = msh.gIntI;
else
    gIntI = options.IntvlInt;
end

% Output arguments.
% Stiffness matrix.
switch options.MatrixType
    case 'full'
        Stiff = zeros(nTstNodes(end),nTrlNodes(end)); 
    case 'sparse'
        Stiff = sparse(nTstNodes(end),nTrlNodes(end));
end
Load = zeros(nTstNodes(end),1); % Load vector;

% Generate stiffness matrix.
% For all double linear functions.
for iAuv = 1:nAuv
    Auv = Auvs(iAuv);
    if Auv.iTrl > nTrl
        error('assemble: "Auv.iTrl" beyond "nTrl".');
    end
    if Auv.iTst > nTst
        error('assemble: "Fv.iTst" beyond "nTst".');
    end
    switch Auv.type
        case 'elem' % If integration domain is element.
            tf = msh.tf; % Transformation.
            Jdet = tf.getJdet.getFun; % Jacobi determinant (function).
            coef = Auv.coef.tfRef(tf).getFun; % Transform coefficient function to reference element.
            iTrl = Auv.iTrl; % Which trial function space.
            iTst = Auv.iTst; % Which test function space.
            trlOrd = Auv.trlOrd; % The order of derivative of trial function.
            tstOrd = Auv.tstOrd; % The order of derivative of test function.
            trl = trls(iTrl); % Trail function space.
            tst = tsts(iTst); % Test function space.
            % 1.Get the base function in mesh element by transformation.
            % 2.Get the derivative of base function.
            % 3.Transform back to reference element.
            trlBases = cell(1,trl.nTf); % All trail base functions (under different transform).
            for iTrlTf = 1:trl.nTf % For different transformation
                trlTf = trl.tf(iTrlTf);
                trlBases{iTrlTf} = trl.base.getFcn.tfOrg(trlTf).dif(trlOrd).tfRef(tf).getFun;
            end
            tstBases = cell(1,tst.nTf); % All test base functions (under different transform).
            for iTstTf = 1:tst.nTf % For different transformation
                tstTf = tst.tf(iTstTf); 
                tstBases{iTstTf} = tst.base.getFcn.tfOrg(tstTf).dif(tstOrd).tfRef(tf).getFun;
            end
            % Combination of trail base function and test base function.
            if isempty(Auv.comb)
                nComb = trl.nBase*tst.nBase;
                comb = zeros(2,nComb);
                iComb = 0;
                for iTrlNode = 1:trl.nBase
                    for iTstNode = 1:tst.nBase
                        iComb = iComb + 1;
                        comb(:,iComb) = [iTrlNode;iTstNode];
                    end
                end
            else
                comb = Auv.comb;
                nComb = size(comb,2);
            end
            % For all combinations of trail base function and test base function.
            for iComb = 1:nComb
                iTrlNode = comb(1,iComb);
                iTstNode = comb(2,iComb);
            % For all base functions in trial function space.
%             for iTrlNode = 1:trl.nBase
                trlBase = trlBases{trl.baseTf(iTrlNode)}; % Trail base function.
                trlCoef = trl.base.coef(:,iTrlNode); % Cofficient of trial base function.
                % For all base functions in test function space.
%                 for iTstNode = 1:tst.nBase
                    tstBase = tstBases{tst.baseTf(iTstNode)}; % Test base function.
                    tstCoef = tst.base.coef(:,iTstNode); % Cofficient of test base funceion.
                    % For all integrated elements.
                    if strcmp(Auv.domn,'all')
                        domn = 1:nElem;
                    else
                        domn = Auv.domn;
                    end
                    for iElem = domn 
                        parm = msh.node(:,msh.elem(:,iElem)); % The parameters about mesh element. (vertices)
                        % Substitue parameters to function (defined in referrece element).
                        switch dim
                            case 1
                                AuvFun = @(x) coef(x,parm).*trlBase(x,parm,trlCoef).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x;
                            case 2
                                AuvFun = @(x) coef(x,parm).*trlBase(x,parm,trlCoef).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x(1)*x(2);
                            case 3
                                AuvFun = @(x) coef(x,parm).*trlBase(x,parm,trlCoef).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x(1)*x(2)*x(3);
                            otherwise
                                error('assemble: dimension exceed 3.');
                        end
                        % Use Gauss integration.
                        AuvVal = gInt.getVal(AuvFun);
                        % Assemble to stiffness matrix.
                        I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                        J = nTrlNodes(iTrl) + abs(trl.elem(iTrlNode,iElem));
                        %! Note: For base function located in connection, the sign may change between connected elements.
                        sgn = sign(tst.elem(iTstNode,iElem)) *  sign(trl.elem(iTrlNode,iElem));
                        Stiff(I,J) = Stiff(I,J) + sgn * AuvVal ;
                    end
%                 end
            end
        case 'edge'
            orgI = msh.tfI.getOrg; % Edge transformation. %! 2 dimension
            JnormI = msh.tfI.getJnorm.getFun; % Jacobi norm of edge transformation.
            coef = Auv.coef.getFun;
            iTrl = Auv.iTrl;
            iTst = Auv.iTst;
            trlOrd = Auv.trlOrd;
            tstOrd = Auv.tstOrd;
            trl = trls(iTrl);
            tst= tsts(iTst);
            trlBases = cell(1,trl.nTf); % All trail base functions (under different transform).
            for iTrlTf = 1:trl.nTf % For different transformation
                trlTf = trl.tf(iTrlTf);
                trlBases{iTrlTf} = trl.base.getFcn.tfOrg(trlTf).dif(trlOrd).getFun;
            end
            tstBases = cell(1,tst.nTf); % All test base functions (under different transform).
            for iTstTf = 1:tst.nTf % For different transformation
                tstTf = tst.tf(iTstTf); 
                tstBases{iTstTf} = tst.base.getFcn.tfOrg(tstTf).dif(tstOrd).getFun;
            end
            for iTrlNode = 1:trl.nBase
                trlBase = trlBases{trl.baseTf(iTrlNode)}; % Trail base function.
                trlCoef = trl.base.coef(:,iTrlNode);
                for iTstNode = 1:tst.nBase
                    tstBase = tstBases{tst.baseTf(iTstNode)}; % Test base function.
                    tstCoef = tst.base.coef(:,iTstNode); 
                    % For all integrated edges.
                    if strcmp(Auv.domn,'all')
                        domn = 1:msh.nEdge;
                    else
                        domn = Auv.domn;
                    end
                    for iEdge = domn 
                        iElem = msh.edge.elem(iEdge);
                        EgParm = msh.node(:,msh.edge.node(:,iEdge)); % Parameter about edge.
                        ElParm = msh.node(:,msh.elem(:,iElem)); % Parameter about element.
                        switch dim %! 2 dimension
                            case 2
                                % Transform to reference interval [-1,1].
%                                 x1 = parm1(1,1); y1 = parm1(2,1); x2 = parm1(1,2); y2 = parm1(2,2);
%                                 AuvFun = @(t) coef([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm1)...
%                                     .*trlBase([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm2)...
%                                     .*tstBase([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm2)...
%                                     .*sqrt((x2-x1)^2+(y2-y1)^2)/2 + 0*t(1);
                                AuvFun = @(t) coef(orgI(t,EgParm),EgParm) .* trlBase(orgI(t,EgParm),ElParm,trlCoef) .* tstBase(orgI(t,EgParm),ElParm,tstCoef) .* JnormI(t,EgParm) + 0*t;
                        end
                        % Use Gauss integration in [-1,1].
                        AuvVal = gIntI.getVal(AuvFun);
                        I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                        J = nTrlNodes(iTrl) + abs(trl.elem(iTrlNode,iElem));
                        %! Note: For base function located in connection, the sign may change between connected elements.
                        sgn = sign(tst.elem(iTstNode,iElem)) *  sign(trl.elem(iTrlNode,iElem));
                        Stiff(I,J) = Stiff(I,J) + sgn* AuvVal;
                    end
                end
            end
        case 'conn'
            orgI = msh.tfI.getOrg; % Edge transformation. %! 2 dimension
            JnormI = msh.tfI.getJnorm.getFun; % Jacobi norm of edge transformation.
            coef = Auv.coef.getFun;
            iTrl = Auv.iTrl;
            iTst = Auv.iTst;
            trlOrd = Auv.trlOrd;
            tstOrd = Auv.tstOrd;
            trl = trls(abs(iTrl));
            tst= tsts(abs(iTst));
            trlBases = cell(1,trl.nTf); % All trail base functions (under different transform).
            for iTrlTf = 1:trl.nTf % For different transformation
                trlTf = trl.tf(iTrlTf);
                trlBases{iTrlTf} = trl.base.getFcn.tfOrg(trlTf).dif(trlOrd).getFun;
            end
            tstBases = cell(1,tst.nTf); % All test base functions (under different transform).
            for iTstTf = 1:tst.nTf % For different transformation
                tstTf = tst.tf(iTstTf); 
                tstBases{iTstTf} = tst.base.getFcn.tfOrg(tstTf).dif(tstOrd).getFun;
            end
            for iTrlNode = 1:trl.nBase
                trlBase = trlBases{trl.baseTf(iTrlNode)}; % Trail base function.
                trlCoef = trl.base.coef(:,iTrlNode);
                for iTstNode = 1:tst.nBase
                    tstBase = tstBases{tst.baseTf(iTstNode)}; % Test base function.
                    tstCoef = tst.base.coef(:,iTstNode); 
                    % For all integrated connections.
                    if strcmp(Auv.domn,'all')
                        domn = 1:msh.nConn;
                    else
                        domn = Auv.domn;
                    end
                    for iConn = domn 
                        cnElem = msh.conn.elem(:,iConn);
                        iTrlElem = abs(cnElem(sign(cnElem)==sign(iTrl))); % Trial function element.
                        iTstElem = abs(cnElem(sign(cnElem)==sign(iTst))); % Test function element.
                        if isempty(iTrlElem) || isempty(iTstElem) % If there is no required connected element, skip this connection.
                            continue;
                        end
                        EgParm = msh.node(:,msh.conn.node(:,iConn)); % Parameter about connection.
                        trlParm = msh.node(:,msh.elem(:,iTrlElem)); % Parameter about trial function element.
                        tstParm = msh.node(:,msh.elem(:,iTstElem)); % Parameter about test function element.
                        switch dim %! 2 dimension
                            case 2
                                % Transform to reference interval [-1,1].
%                                 x1 = parm1(1,1); y1 = parm1(2,1); x2 = parm1(1,2); y2 = parm1(2,2);
%                                 AuvFun = @(t) coef([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm1)...
%                                     .*trlBase([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm2)...
%                                     .*tstBase([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm2)...
%                                     .*sqrt((x2-x1)^2+(y2-y1)^2)/2 + 0*t(1);
                                AuvFun = @(t) coef(orgI(t,EgParm),EgParm) .* trlBase(orgI(t,EgParm),trlParm,trlCoef) .* tstBase(orgI(t,EgParm),tstParm,tstCoef) .* JnormI(t,EgParm) + 0*t;
                        end
                        % Use Gauss integration in [-1,1].
                        AuvVal = gIntI.getVal(AuvFun);
                        I = nTstNodes(abs(iTst)) + abs(tst.elem(iTstNode,iTstElem));
                        J = nTrlNodes(abs(iTrl)) + abs(trl.elem(iTrlNode,iTrlElem));
                        %! Note: For base function located in connection, the sign may change between connected elements.
                        sgn = sign(tst.elem(iTstNode,iTstElem)) *  sign(trl.elem(iTrlNode,iTrlElem));
                        Stiff(I,J) = Stiff(I,J) + sgn* AuvVal;
                    end
                end
            end
        case 'refElem' % If integration domain is reference element.
            tf = msh.tf; % Transformation.
            coef = Auv.coef.tfRef(tf).getFun; % Transform coefficient function to reference element.
            iTrl = Auv.iTrl; % Which trial function space.
            iTst = Auv.iTst; % Which test function space.
            trlOrd = Auv.trlOrd; % The order of derivative of trial function.
            tstOrd = Auv.tstOrd; % The order of derivative of test function.
            trl = trls(iTrl); % Trail function space.
            tst = tsts(iTst); % Test function space.
            trlBase = trl.base.getFcn.dif(trlOrd).getFun; % Trail base function.
            tstBase = tst.base.getFcn.dif(tstOrd).getFun; % Test base function.
            % Combination of trail base function and test base function.
            if isempty(Auv.comb)
                nComb = trl.nBase*tst.nBase;
                comb = zeros(2,nComb);
                iComb = 0;
                for iTrlNode = 1:trl.nBase
                    for iTstNode = 1:tst.nBase
                        iComb = iComb + 1;
                        comb(:,iComb) = [iTrlNode;iTstNode];
                    end
                end
            else
                comb = Auv.comb;
                nComb = size(comb,2);
            end
            % For all combinations of trail base function and test base function.
            for iComb = 1:nComb
                iTrlNode = comb(1,iComb);
                iTstNode = comb(2,iComb);
            % For all base functions in trial function space.
%             for iTrlNode = 1:trl.nBase
                trlCoef = trl.base.coef(:,iTrlNode); % Cofficient of trial base function.
                % For all base functions in test function space.
%                 for iTstNode = 1:tst.nBase
                    tstCoef = tst.base.coef(:,iTstNode); % Cofficient of test base funceion.
                    % For all integrated elements.
                    if strcmp(Auv.domn,'all')
                        domn = 1:nElem;
                    else
                        domn = Auv.domn;
                    end
                    for iElem = domn 
                        parm = msh.node(:,msh.elem(:,iElem)); % The parameters about mesh element. (vertices)
                        % Substitue parameters to function (defined in referrece element).
                        switch dim
                            case 1
                                AuvFun = @(x) coef(x,parm).*trlBase(x,[],trlCoef).*tstBase(x,[],tstCoef) + 0*x;
                            case 2
                                AuvFun = @(x) coef(x,parm).*trlBase(x,[],trlCoef).*tstBase(x,[],tstCoef) + 0*x(1)*x(2);
                            case 3
                                AuvFun = @(x) coef(x,parm).*trlBase(x,[],trlCoef).*tstBase(x,[],tstCoef) + 0*x(1)*x(2)*x(3);
                            otherwise
                                error('assemble: dimension exceed 3.');
                        end
                        % Use Gauss integration.
                        AuvVal = gInt.getVal(AuvFun);
                        % Assemble to stiffness matrix.
                        I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                        J = nTrlNodes(iTrl) + abs(trl.elem(iTrlNode,iElem));
                        %! Note: For base function located in connection, the sign may change between connected elements.
                        sgn = sign(tst.elem(iTstNode,iElem)) *  sign(trl.elem(iTrlNode,iElem));
                        Stiff(I,J) = Stiff(I,J) + sgn * AuvVal ;
                    end
%                 end
            end
    end
end
% For all single linear functions.
for iFv = 1:nFv
    Fv = Fvs(iFv);
    switch Fv.type
        case 'elem' % If integration domain is element.
            tf = msh.tf; % Transformation.
            Jdet = tf.getJdet.getFun; % Jacobi determinant (function).
            load = Fv.load.tfRef(tf).getFun; % Transform load function to referrence element.
            iTst = Fv.iTst; % Which test function space.
            tstOrd = Fv.tstOrd; % The order of derivative of test function.
            tst = tsts(iTst); % Test space.
            tstBases = cell(1,tst.nTf); % All test base functions (under different transform).
            for iTstTf = 1:tst.nTf % For different transformation
                tstTf = tst.tf(iTstTf); 
                tstBases{iTstTf} = tst.base.getFcn.tfOrg(tstTf).dif(tstOrd).tfRef(tf).getFun;
            end
            % For all base functions in test function space.
            for iTstNode = 1:tst.nBase
                tstBase = tstBases{tst.baseTf(iTstNode)}; % Test base function.
                tstCoef = tst.base.coef(:,iTstNode); % Cofficient of test base funceion.
                % For all integrated element.
                if strcmp(Fv.domn,'all')
                    domn = 1:nElem;
                else
                    domn = Fv.domn;
                end
                for iElem = domn
                    parm = msh.node(:,msh.elem(:,iElem)); % The parameters about mesh element. (vertices)
                    % Substitue parameters to function (defined in referrece element).
                    switch dim
                        case 1
                            FvFun = @(x) load(x,parm).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x(1);
                        case 2
                            FvFun = @(x) load(x,parm).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x(1)*x(2);
                        case 3
                            FvFun = @(x) load(x,parm).*tstBase(x,parm,tstCoef).*Jdet(x,parm) + 0*x(1)*x(2)*x(3);
                        otherwise
                    end
                    % Use Gauss integration.
                    FvVal = gInt.getVal(FvFun);
                    % Assemble to stiffness matrix.
                    I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                    sgn = sign(tst.elem(iTstNode,iElem));
                    Load(I) = Load(I) + sgn * FvVal;
                end
            end
        case 'edge' % If integration domain is edge.
            orgI = msh.tfI.getOrg; % Edge transformation. %! 2 dimension
            JnormI = msh.tfI.getJnorm.getFun; % Jacobi norm of edge transformation.
            load = Fv.load.getFun;
            iTst = Fv.iTst;
            tstOrd = Fv.tstOrd;
            tst = tsts(iTst);
            tstBases = cell(1,tst.nTf); % All test base functions (under different transform).
            for iTstTf = 1:tst.nTf % For different transformation
                tstTf = tst.tf(iTstTf); 
                tstBases{iTstTf} = tst.base.getFcn.tfOrg(tstTf).dif(tstOrd).getFun;
            end
            for iTstNode = 1:tst.nBase
                tstBase = tstBases{tst.baseTf(iTstNode)}; % Test base function.
                tstCoef = tst.base.coef(:,iTstNode); 
                % For all integrated edge element.
                if strcmp(Fv.domn,'all')
                    domn = 1:msh.nEdge;
                else
                    domn = Fv.domn;
                end
                for iEdge = domn
                    iElem = msh.edge.elem(iEdge);
                    EgParm = msh.node(:,msh.edge.node(:,iEdge)); % Parameter about edge.
                    ElParm = msh.node(:,msh.elem(:,iElem));  % Parameter about element.
                    switch dim %! 2 dimension
                        case 1
                        case 2
                            % Transform to reference interval [-1,1]
%                             x1 = parm1(1,1); y1 = parm1(2,1); x2 = parm1(1,2); y2 = parm1(2,2);
%                             FvFun = @(t) load([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm1).*...
%                                 tstBase([(x2-x1)/2*t+(x2+x1)/2;(y2-y1)/2*t+(y2+y1)/2],parm2)... 
%                                 .* sqrt((x2-x1)^2+(y2-y1)^2)/2 + 0*t(1);
                            FvFun = @(t) load(orgI(t,EgParm),EgParm) .* tstBase(orgI(t,EgParm),ElParm,tstCoef) .* JnormI(t,EgParm) + 0*t;
                        case 3
                    end
                    % Use Gauss integration in [-1,1]
                    FvVal = gIntI.getVal(FvFun);
                    I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                    sgn = sign(tst.elem(iTstNode,iElem));
                    Load(I) = Load(I) + sgn * FvVal;
                end
            end   
        case 'refElem' % If integration domain is element.
            tf = msh.tf; % Transformation.
            load = Fv.load.tfRef(tf).getFun; % Transform load function to referrence element.
            iTst = Fv.iTst; % Which test function space.
            tstOrd = Fv.tstOrd; % The order of derivative of test function.
            tst = tsts(iTst); % Test space.
            tstBase = tst.base.getFcn.dif(tstOrd).getFun; % Test base function.
            % For all base functions in test function space.
            for iTstNode = 1:tst.nBase
                tstCoef = tst.base.coef(:,iTstNode); % Cofficient of test base funceion.
                % For all integrated element.
                if strcmp(Fv.domn,'all')
                    domn = 1:nElem;
                else
                    domn = Fv.domn;
                end
                for iElem = domn
                    parm = msh.node(:,msh.elem(:,iElem)); % The parameters about mesh element. (vertices)
                    % Substitue parameters to function (defined in referrece element).
                    switch dim
                        case 1
                            FvFun = @(x) load(x,parm).*tstBase(x,[],tstCoef) + 0*x(1);
                        case 2
                            FvFun = @(x) load(x,parm).*tstBase(x,[],tstCoef) + 0*x(1)*x(2);
                        case 3
                            FvFun = @(x) load(x,parm).*tstBase(x,[],tstCoef) + 0*x(1)*x(2)*x(3);
                        otherwise
                    end
                    % Use Gauss integration.
                    FvVal = gInt.getVal(FvFun);
                    % Assemble to stiffness matrix.
                    I = nTstNodes(iTst) + abs(tst.elem(iTstNode,iElem));
                    sgn = sign(tst.elem(iTstNode,iElem));
                    Load(I) = Load(I) + sgn * FvVal;
                end
            end
    end
end

end
