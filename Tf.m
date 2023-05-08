classdef Tf
    % Tf: Transformation between original element and reference element.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com

    properties
        refVar sym; % Reference variables.(By coloumn)
        orgVar sym; % Original variables.(By coloumn)
        parm (:,:) sym; % Parameters. The coordinates of nodes of original mesh element (By column).
        toRef sym; % Transformation from orgrinal element to reference element.(By coloumn)
        toOrg sym; % Transformation form reference element to original element.(By coloumn)
        refCoef sym = 1; % Coefficient of transformation to reference element. (Scalar or Matrix)
        orgCoef sym = 1; % Coefficient of transformation to origianl element. (Scalar or Matrix)
    end
    properties (Dependent)
        nVar; % The number of variables.
        nParm; % The size of parameters.
    end
    methods
        function disp(tf)
            nTf = length(tf);
            if nTf == 1
                fprintf('refVar: \n'), disp(tf.refVar);
                fprintf('orgVar: \n'), disp(tf.orgVar);
                fprintf('parm: \n'), disp(tf.parm);
                fprintf('toRef: \n'), disp(tf.toRef);
                fprintf('toOrg: \n'), disp(tf.toOrg);
                fprintf('refCoef: \n'), disp(tf.refCoef);
                fprintf('orgCoef: \n'), disp(tf.orgCoef);
            else
                for iTf = 1:nTf
                    divLn('=',64,[inputname(1),'(',num2str(iTf),')']);
                    fprintf('refVar: \n'), disp(tf(iTf).refVar);
                    fprintf('orgVar: \n'), disp(tf(iTf).orgVar);
                    fprintf('parm: \n'), disp(tf(iTf).parm);
                    fprintf('toRef: \n'), disp(tf(iTf).toRef);
                    fprintf('toOrg: \n'), disp(tf(iTf).toOrg);
                    fprintf('refCoef: \n'), disp(tf(iTf).refCoef);
                    fprintf('orgCoef: \n'), disp(tf(iTf).orgCoef);
                end
                divLn('=',64);
            end
        end
        function tf = set.refVar(tf,refVar)
            if ~isvector(refVar)
                error('Tf: "refVar" must be a vector.');
            end
            tf.refVar = refVar;
        end
        function tf = set.orgVar(tf,orgVar)
            if ~isvector(orgVar)
                error('Tf: "orgVar" must be a vector.');
            end
            tf.orgVar = orgVar;
        end
        function tf = set.toRef(tf,toRef)
            if ~isvector(toRef)
                error('Tf: "toRef" must be a vector.');
            end
            tf.toRef = toRef;
        end
        function tf = set.toOrg(tf,toOrg)
            if ~isvector(toOrg)
                error('Tf: "toOrg" must be a vector.');
            end
            tf.toOrg = toOrg;
        end
        function nVar = get.nVar(tf)
            if any(diff( [length(tf.refVar), length(tf.orgVar), length(tf.toRef), length(tf.toOrg), size(tf.parm,1)] ) ~= 0)
                error('Tf: the number of variables is inconsistent - "refVar", "orgVar", "toRef", "toVar".');
            end
            nVar = length(tf.toRef);
        end
        function nParm = get.nParm(tf)
            nParm = size(tf.parm);
        end

        function Jdet = getJdet(tf)
            % Get the determinant of the Jacobi matrix of the transformation from original element to reference element.
            Jdet = Fcn;
            Jdet.var = tf.refVar;
            Jdet.parm = tf.parm;
            Jdet.fun = det(jacobian(tf.toRef,tf.refVar));
        end
        function Jnorm = getJnorm(tf)
            % Get the norm of the Jacobi matrix of the transformation from original element to reference element.
            Jnorm = Fcn;
            Jnorm.var = tf.refVar;
            Jnorm.parm = tf.parm;
            Jnorm.fun = norm(jacobian(tf.toRef,tf.refVar));
        end
%         function Jdet = getJ(tf)
%             % Get the function of determinant of the Jacobi matrix of the transformation from original element to reference element.
%             % Jdet = @(refVar(1),refVar(2),...,parm)
%             vars = {tf.refVar,tf.parm};
% %             vars = [num2cell(tf.refVar),{tf.parm}];
%             Jdet = matlabFunction(det(jacobian(tf.toRef,tf.refVar)),'Vars',vars);
%         end
        function ref = getRef(tf)
           % Get the function which transform original element to reference element.
           vars = {tf.orgVar,tf.parm};
%            vars = [num2cell(tf.orgVar),{tf.parm}];
           ref = matlabFunction(tf.toOrg,'Vars',vars);
        end
        function org = getOrg(tf)
           % Get the function which transform reference element to original element.
           vars = {tf.refVar,tf.parm};
%            vars = [num2cell(tf.refVar),{tf.parm}];
           org = matlabFunction(tf.toRef,'Vars',vars);
        end
        function check(tf)
            tf.nVar;
        end
    end
end

