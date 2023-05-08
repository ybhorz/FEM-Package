classdef Fcn
    % Fcn: Function.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com
    
    % Fcn is a vector space over sym. Fcn = {f1(var,parm)*c1 + f2(var,parm)*c2 + ...}
    % Operator: 
    %   (Fcn,sym) -> Fcn : * \
    %   (Fcn,Fcn) -> Fcn : + - .* .\
    %   Fcn -> Fcn : - .^
    properties
        var sym; % Variables.(By coloumn)
        parm (:,:) sym; % Parameters.(By coloumn)
        fun sym; % Function.(By coloumn)
                 % Note: may be vector-function.
        coef sym; % Coefficient.(By coloumn)
    end
    properties (Dependent)
        nVar; % The number of variables.
        nParm; % The size of parameters.
        nRang; % The length of vetor-valued function.
    end
    methods
        function fcn = Fcn(var,fun,parm)
            if nargin >= 2
                fcn.var = var;
                fcn.fun = fun;
                if nargin >= 3
                    fcn.parm = parm;
                end
            end
        end
        function disp(fcn)
            nFcn = length(fcn);
            if nFcn == 1
                fprintf('var: \n'), disp(fcn.var);
                if ~isempty(fcn.parm)
                   fprintf('parm: \n'), disp(fcn.parm);
                end
                fprintf('fun: \n'), disp(fcn.fun);
                if ~isempty(fcn.coef)
                   fprintf('coef: \n'), disp(fcn.coef);
                end
            else
                for iFcn = 1:nFcn
                    divLn('=',64,[inputname(1),'(',num2str(iFcn),')']);
                    fprintf('var: \n'), disp(fcn(iFcn).var);
                    if ~isempty(fcn(iFcn).parm)
                       fprintf('parm: \n'), disp(fcn(iFcn).parm);
                    end
                    fprintf('fun: \n'), disp(fcn(iFcn).fun);
                    if ~isempty(fcn(iFcn).coef)
                       fprintf('coef: \n'), disp(fcn(iFcn).coef);
                    end
                end
                divLn('=',64);
            end
        end
        function fcn = set.var(fcn,var)
            if ~isvector(var)
                error('Fcn: "var" must be a vector.');
            end
            fcn.var = var;
        end
        function fcn = set.fun(fcn,fun)
            if isvector(fun)
                fcn.fun = fun;
            else
                error('Fcn: "fun" must be vector.');
            end
        end
        function nVar = get.nVar(fcn)
            nVar = length(fcn.var);
        end
        function nParm = get.nParm(fcn)
            nParm = size(fcn.parm);
        end
        function nRang = get.nRang(fcn)
            nRang = length(fcn.fun);
        end
        function fcn = tfOrg(fcn,tf)
            % Transform to original element.
            if ~isa(tf,'Tf')
                error('Fcn.tfOrg: "tf" must be "Tf".');
            end
            fcn = fcn.sub(tf.toOrg,tf.orgVar,tf.parm,tf.orgCoef);
        end
        function fcn = tfRef(fcn,tf)
            % Transform to reference element.
            if ~isa(tf,'Tf')
                error('Fcn.tfRef: "tf" must be "Tf".');
            end
            fcn = fcn.sub(tf.toRef,tf.refVar,tf.parm,tf.refCoef);
        end
        function fcn = dif(fcn,order)
            % Get the derivative of "order". (By coloumn)
            % Note: function may be vector-valued, and result is the sum of derivates of every components.
            % Note: if non-positive integer appears, the corresponding components will be omitted.
            if size(order,1) ~= fcn.nVar || size(order,2) ~= fcn.nRang
                error('Fcn.dif: "order" must be compatible with "fun".');
            end
            for iRang = 1:fcn.nRang
                if all(order(:,iRang) >= 0)
                    for iVar = 1:fcn.nVar
                        fcn.fun(iRang) = diff(fcn.fun(iRang),fcn.var(iVar),order(iVar,iRang));
                    end
                else
                    fcn.fun(iRang) = 0;
                    continue;
                end
            end
            fcn.fun = sum(fcn.fun);
        end
        function func = getFun(fcn)
            % Get the function.
            % func = @(var(1),var(2),...,parm)
%             vars = [num2cell(fcn.var),{fcn.parm}];
            if isempty(fcn.coef)
                vars = {fcn.var,fcn.parm};
            else
                vars = {fcn.var,fcn.parm,fcn.coef};
            end
            func = matlabFunction(fcn.fun,'Vars',vars);
        end
        function check(fcn)
            
        end
%%%%%%%%%%%%%%%%% Overload some functions %%%%%%%%%%%%%%%%%
        function fcn = subs(fcn,var,val)
            % Overload subs.
            fcn.fun = subs(fcn.fun,var,val);
        end
        function fcn = diff(fcn,var,order)
            % Overload diff.
            if nargin == 2
                fcn.fun = diff(fcn.fun,var);
            else
                fcn.fun = diff(fcn.fun,var,order);
            end
        end
        function fcn = dot(fcn1,fcn2)
            % Overload dot.
            fcn = fcn1;
            fcn.fun = sum(fcn1.fun.*fcn2.fun);
        end
%%%%%%%%%%%%%%%%% Overload some operators %%%%%%%%%%%%%%%%%
        function fcn = plus(fcn1,fcn2)
            % Overload + operator.
            %? 如何判断符号变量相同（变量、参数）
            fcn = Fcn;
            fcn.var = fcn1.var;
            fcn.parm = fcn1.parm;
            fcn.fun = fcn1.fun + fcn2.fun;
            fcn.coef = [fcn1.coef;fcn2.coef];
        end
        function fcn = minus(fcn1,fcn2)
            % Overload - operator.
            %? 如何判断符号变量相同（变量、参数）
            fcn = Fcn;
            fcn.var = fcn1.var;
            fcn.parm = fcn1.parm;
            fcn.fun = fcn1.fun - fcn2.fun;
            fcn.coef = [fcn1.coef;fcn2.coef];
        end
        function fcn = uminus(fcn)
            % Overload - operator.
            fcn.fun = -fcn.fun;
        end
        function fcn = mtimes(fcn,coef)
            % Overload * operator.
            % Note: fcn * coef
            if isnumeric(coef)
                fcn.fun = fcn.fun*coef;
            elseif isa(coef,"sym")
                fcn.fun = fcn.fun*coef;
                fcn.coef = [fcn.coef;coef];
            else
                error('Fcn: "coef" must be number or sym.');
            end
        end
        function fcn = mldivide(fcn,coef)
            % Overload \ operator.
            % Note: fcn \ coef
            if isnumeric(coef)
                fcn.fun = fcn.fun\coef;
%             elseif isa(coef,"sym")
%                 fcn.fun = fcn.fun\coef;
%                 fcn.coef = [fcn.coef;coef];
            else
                error('Fcn: "coef" must be number.');
            end
        end
        function fcn = times(fcn1,fcn2)
            % Overload .* operator.
            fcn = Fcn;
            fcn.var = fcn1.var;
            fcn.parm = fcn1.parm;
            fcn.fun = fcn1.fun .* fcn2.fun;
            fcn.coef = fcn1.coef;
        end
        function fcn = ldivide(fcn1,fcn2)
            % Overload .\ operator.
            fcn = Fcn;
            fcn.var = fcn1.var;
            fcn.parm = fcn1.parm;
            fcn.fun = fcn1.fun .\ fcn2.fun;
            fcn.coef = fcn1.coef;
        end
        function fcn = power(fcn,r)
            % Overload .^ operator.
            fcn.fun = (fcn.fun).^r;
        end
        function fcn = vertcat(fcn1,fcn2)
            % Overload [ ; ] operator
            fcn = fcn1;
            fcn.fun = [fcn1.fun;fcn2.fun];
        end
    end
    methods (Access = private)
        function fcn = sub(fcn,toNew,newVar,parm,coef)
            % Replace "var" in "fun" with "toNew", new variables is "newVar", parameters is "parm"
            fcn.fun = coef * subs(fcn.fun,fcn.var,toNew);
            fcn.var = newVar;
            fcn.parm = parm;
        end
    end
end

