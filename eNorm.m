function eNorm = eNorm(msh,fE,appr,real,order,options)
% eNorm: Get the norm of error between approximation and real solutions. 
%   Note: Approxiamtion is represeted as the combination of coefficients and base functions in finte element space.
% Author: Bohan Yang
% Email: ybh.orz@gmail.com

% Input arguments.
arguments
    msh Msh; % Mesh information.
    fE FE; % Finie element space.
    appr; % Approximation solution. (Coefficients of the base functions)
    real Fcn; % Real solutions.
    order; % Order of derivative.
           % For more than one set, use cell array.
    options.norm = 2; % p-norm. (Include inf-norm)
    options.domain = 'all'; % Domain
    options.type {mustBeMember(options.type,{'absolute','relative'})} = 'absolute'; % Absolute error or relative error.
end
% if ~isa(msh,'Msh')
%     error('eNorm: "msh" must be "Msh".');
% end
% if ~isa(fE,'FE')
%     error('eNorm: "fE" must be "FE".');
% end
% if ~isa(real,'Fcn')
%     error('eNorm: "real" must be "Fcn".');
% end

% Checking.
msh.check;
fE.check;
real.check;
if ~isvector(appr)
    error('eNorm: "coef" must be a vector.');
end

dim = msh.dim;
if any(diff([dim,fE.dim,real.nVar]))
    error('eNorm: the dimension is inconsistent - "msh", "fE", "real".');
end

nElem = msh.nElem; 
if nElem ~= fE.nElem
    error('eNorm: the number is inconsisten - "msh", "fE".');
end

if fE.nNode ~= length(appr)
    error('eNorm: the number is inconsistent - "fEs", "coef".');
end

% Default options.
p = options.norm;
if strcmp(options.domain,'all')
    domn = 1:nElem;
else
    domn = options.domain;
end


% Error norms witch different orders. 
if iscell(order)
    if isvector(order)
        nOrder = length(order);
    else
        error('eNorm: "order" must be vector.');
    end
else
    order = {order};
    nOrder = 1;
end
eNorms = zeros(1,nOrder);
uNorms = zeros(1,nOrder);

% Transformation
tf = msh.tf;
Jdet = tf.getJdet.getFun;
% Symblic cofficients of base funcions.
nBase = fE.nBase;
c = sym('cFun',[1,nBase]);
% For all orders.
for iOrder = 1:nOrder
    % Get error function (defined in reference element).
    uFcn = real.dif(order{iOrder}).tfRef(tf);
    err = uFcn; % e = u
    for iBase = 1:nBase
        fETf = fE.tf(fE.baseTf(iBase));
        err = err - fE.base.getFcn(iBase).tfOrg(fETf).dif(order{iOrder}).tfRef(tf) * c(iBase); % e = u - sum c_i phi_i
    end
    err = err.getFun;
    uFcn = uFcn.getFun;
    % For all elements.
    for iElem = domn
        if isinf(p) % inf-norm
            parm = msh.node(:,msh.elem(:,iElem));
            coef = sign(fE.elem(:,iElem)) .* appr(abs(fE.elem(:,iElem)));
            switch dim
                case (1)
                    errFun = @(x) abs(err(x,parm,coef)) + 0*x(1);
                case (2)
                    errFun = @(x) abs(err(x,parm,coef)) + 0*x(1)*x(2);
                case (3)
                    errFun = @(x) abs(err(x,parm,coef)) + 0*x(1)*x(2)*x(3);
            end
            uFun = @(x) abs(uFcn(x,parm));
            % Take the max error on Gauss integration points.
            for iPoint = 1:msh.gInt.nPoint
                eNorms(iOrder) = max(eNorms(iOrder),errFun(msh.gInt.point(:,iPoint)));
                uNorms(iOrder) = max(uNorms(iOrder),uFun(msh.gInt.point(:,iPoint)));
            end
        else % p-norm
            parm = msh.node(:,msh.elem(:,iElem));
            coef = sign(fE.elem(:,iElem)) .* appr(abs(fE.elem(:,iElem)));
            switch dim
                case (1)
                    errFun = @(x) abs(err(x,parm,coef)).^p .*Jdet(x,parm) + 0*x(1);
                case (2)
                    errFun = @(x) abs(err(x,parm,coef)).^p .*Jdet(x,parm) + 0*x(1)*x(2);
                case (3)
                    errFun = @(x) abs(err(x,parm,coef)).^p .*Jdet(x,parm) + 0*x(1)*x(2)*x(3);
            end
            uFun = @(x) abs(uFcn(x,parm)).^p .*Jdet(x,parm);
            eNorms(iOrder) = eNorms(iOrder) + msh.gInt.getVal(errFun);
            uNorms(iOrder) = uNorms(iOrder) + msh.gInt.getVal(uFun);
        end
    end
end
% Total error norm.
if isinf(p)
    eNorm = sum(abs(eNorms));
    uNorm = sum(abs(uNorms));
else
    eNorm = sum(eNorms)^(1/p);
    uNorm = sum(uNorms)^(1/p);
end

if strcmp(options.type,'relative')
    eNorm = eNorm/uNorm;
end

end