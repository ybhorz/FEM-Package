classdef Base
    % Base: Base function
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com
    
    properties
        var sym; % Variables.(By coloumn)
        Fs sym; % Finite function spase.(By column)
        coef;   % Base functions are represented as the linear combination of "FS" and "coef".
                % Base(i) = Fs * coef(:,i)
    end
    properties (Dependent)
        nVar; % The number of variables.
        nRang; % The length of vector-valued function.
        nFun; % The number of functions.
        nBase; % The number of base functions.
    end
    methods
        function base = Base(var,FS,coef)
            if nargin >= 3
                base.var = var;
                base.Fs = FS;
                base.coef = coef;
            end
        end
        function disp(base)
            for iBase = 1:base.nBase
                divLn('=',64,['base(',num2str(iBase),')']);
                disp(base.Fs*base.coef(:,iBase));
            end
        end
        function base = set.var(base,var)
            if ~isvector(var)
                error('Base: "var" must be a vector.');
            end
            base.var = var;
        end
        function nVar = get.nVar(base)
            nVar = length(base.var);
        end
        function nRang = get.nRang(base)
            nRang = size(base.Fs,1);
        end
        function nFun = get.nFun(base)
            nFun = size(base.Fs,2);
        end
        function nBase = get.nBase(base)
            nBase = size(base.coef,2);
        end
        function baseFcn = getFcn(base,iBase)
            % Get the Fcn of base(i)
            baseFcn = Fcn;
            baseFcn.var = base.var;
            if nargin == 1
                cBase = sym('cBase',[base.nFun,1]);
                baseFcn.fun = base.Fs * cBase;
                baseFcn.coef = cBase;
            elseif nargin == 2
                baseFcn.fun = base.Fs * base.coef(:,iBase);
            end
        end
        function base = plus(base1,base2)
            % Overload + operator
            base = Base;
            base.var = base1.var;
            base.Fs = [base1.Fs,base2.Fs];
            base.coef = blkdiag(base1.coef,base2.coef);
        end
        function base = power(base,r)
            % Overload .^ operator.
            if r ~= 2
                error('Base.power: "r" must be "2".');
            end
            base.Fs = [base.Fs,zeros(size(base.Fs));zeros(size(base.Fs)),base.Fs];
            base.coef = [base.coef,zeros(size(base.coef));zeros(size(base.coef)),base.coef];
        end
        function check(base)
            if size(base.Fs,2) ~= size(base.coef,1)
                error('Base: the size of inconsistent - "Fs", "coef".');
            end
        end
    end
end

