classdef SLF
    % SLF: Single linear function.
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com

    properties
        load Fcn = Fcn; % Load function.
        iTst = 1; % The indice of the test function space. (which variable)
        tstOrd; % The order of the derivative of test function. (By coloumn)
        type {mustBeMember(type,{'elem','edge','refElem'})} = 'elem'; % Integration type: elem, edge,'refElem'
        domn = 'all'; % Integration domain (Indices about elem, edge). 
    end
    properties (Dependent)
        dim; % Dimension;
    end
    methods
        function sLF = SLF(load,iTst,tstOrd,type,domn)
            if nargin == 2
                % "iTst" represent "tstOrd"
                sLF.load = load;
                sLF.iTst = 1;
                sLF.tstOrd = iTst;
            elseif nargin >= 3
                sLF.load = load;
                sLF.iTst = iTst;
                sLF.tstOrd = tstOrd;
                if nargin >= 4
                    sLF.type = type;
                    if nargin >= 5
                        sLF.domn = domn;
                    end
                end
            end
        end
        function disp(sLF)
            for iSLF = 1:length(sLF)
                divLn('#',64,[inputname(1),'(',num2str(iSLF),')']);
                divLn('-',64,'load');
                disp(sLF(iSLF).load);
                divLn('-',64);
                fprintf('iTst: \n'), disp(sLF(iSLF).iTst);
                fprintf('tstOrd: \n'), disp(sLF(iSLF).tstOrd);
                fprintf('type: %s\n',sLF(iSLF).type);
                if ~strcmp(sLF(iSLF).domn,'all')
                    fprintf('domn: [1,%d]\n',length(sLF(iSLF).domn));
                end
            end
            divLn('#',64);
        end
%         function sLF = set.tstOrd(sLF,tstOrd)
%             if ~isvector(tstOrd)
%                 error('SLF: "tstOrd" must be integer vector.');
%             end
%             sLF.tstOrd = tstOrd;
%         end
%         function sLF = set.domn(sLF,domn)
%             if isvector(domn)
%                 sLF.domn = domn;
%             else
%                 error('SLF: "domn" must be vector.');
%             end
%         end
        function dim = get.dim(sLF)
            if sLF.load.nVar == size(sLF.tstOrd,1)
                dim = sLF.load.nVar;
            else
                error('SLF: dimension is inconsistent: "load", "tstOrd".');
            end
        end
        function sLF = subs(sLF,var,val)
            % Overload subs.
            for i = 1:length(sLF)
                sLF(i).load = sLF(i).load.subs(var,val);
            end
        end
        function check(sLF)
            sLF.load.check;
        end
    end
end

