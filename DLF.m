classdef DLF
    % DLF: Double linear function.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com

    properties
        coef Fcn = Fcn; % Coefficient function.
        iTrl = 1; % The indice of the trial function space. (which variable)
                  %! Note: when integraiton type is "conn", 
                  %         positive value means trial function is taken from positive connected element.
                  %         negative value means trial function is taken from negative connected element.
        iTst = 1; % The indice of the test function space. (which variable)
                  %! Note: when integraiton type is "conn", 
                  %         positive value means test function is taken from positive connected element.
                  %         negative value means test function is taken from negative connected element.
        trlOrd (:,:); % The order of the derivative of trial function. (By coloumn)
                      % Note: trial function may be vector-valued, and result is the sum of derivates of every components.
                      % Note: if non-positive integer appears, the corresponding components will be omitted.
        tstOrd (:,:); % The order of the derivative of test function. (By coloumn)
                      % Note: test function may be vector-valued, and result is the sum of derivates of every components.
                      % Note: if non-positive integer appears, the corresponding components will be omitted.
        type {mustBeMember(type,{'elem','edge','conn','refElem'})} = 'elem'; % Integration type: elem, edge, conn, refElem
        domn = 'all'; % Integration domain (Indices about elem, edge, conn). 
        comb; % The combination of trail base function and test base function. (By coloumn);
    end
    properties (Dependent)
        dim; % Dimension.
    end
    methods
        function dLF = DLF(coef,iTrl,iTst,trlOrd,tstOrd,type,domn,comb)
            if nargin == 3
                % "iTrl" represent "trlOrd"
                % "iTst" represent "tstOrd"
                dLF.coef = coef;
                dLF.iTrl = 1;
                dLF.iTst = 1;
                dLF.trlOrd = iTrl;
                dLF.tstOrd = iTst; 
            elseif nargin >= 5
                dLF.coef = coef;
                dLF.iTrl = iTrl;
                dLF.iTst = iTst;
                dLF.trlOrd = trlOrd;
                dLF.tstOrd = tstOrd;
                if nargin >= 6
                    dLF.type = type;
                    if nargin >= 7
                        dLF.domn = domn;
                        if nargin >= 8
                            dLF.comb = comb;
                        end
                    end
                end
            end
        end
        function disp(dLF)
            for iDLF = 1:length(dLF)
                divLn('#',64,[inputname(1),'(',num2str(iDLF),')']);
                divLn('-',64,'coef');
                disp(dLF(iDLF).coef);
                divLn('-',64);
                fprintf('iTrl: \n'), disp(dLF(iDLF).iTrl);
                fprintf('iTst: \n'), disp(dLF(iDLF).iTst);
                fprintf('trlOrd: \n'), disp(dLF(iDLF).trlOrd);
                fprintf('tstOrd: \n'), disp(dLF(iDLF).tstOrd);
                fprintf('type: %s\n',dLF(iDLF).type);
                if ~strcmp(dLF(iDLF).domn,'all')
                    fprintf('domn: [1,%d]\n',length(dLF(iDLF).domn));
                end
                if ~isempty(dLF(iDLF).comb)
                    fprintf('comb: \n'), disp(dLF(iDLF).comb);
                end
            end
            divLn('#',64);
        end
%         function dLF = set.trlOrd(dLF,trlOrd)
%             if ~isvector(trlOrd)
%                 error('DLF: "trlOrd" must be integer vector.');
%             end
%             dLF.trlOrd = trlOrd;
%         end
%         function dLF = set.tstOrd(dLF,tstOrd)
%             if ~isvector(tstOrd)
%                 error('DLF: "tstOrd" must be integer vector.');
%             end
%             dLF.tstOrd = tstOrd;
%         end
%         function dLF = set.domn(dLF,domn)
%             if isvector(domn)
%                 dLF.domn = domn;
%             else
%                 error('DLF: "domn" must be vector.');
%             end
%         end
        function dim = get.dim(dLF)
            if any( diff([dLF.coef.nVar, size(dLF.trlOrd,1), size(dLF.tstOrd,1)]) ~= 0 )
                error('DLF: the dimension is inconsistent - "coef", "trlOrd", "tstOrd".');
            end
            dim = dLF.coef.nVar;
        end

        function dLF = subs(dLF,var,val)
            % Take var = val in coef.
            for i = 1:length(dLF)
                dLF(i).coef = dLF(i).coef.subs(var,val);
            end
        end
        function check(dLF)
            dLF.coef.check;
            dLF.dim;
        end
    end
end

