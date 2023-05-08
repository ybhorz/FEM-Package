classdef GInt
    % GInt: Gauss numerical integration  
    % Author: Bohan Yang
    % Email: ybh.orz@gmail.com

    properties
        point (:,:); % (By coloumn)
        coef;
    end
    properties(Dependent)
        dim;
        nPoint;
    end
    methods
        function gInt = GInt(type)
            % D: dimension
            % P: polynomial degree
            % Integration domain: 1D [-1,1]
            %                     2D triangle with vetices (0,0) (1,0) (0,1)  
            arguments
                type {mustBeMember(type,{'void','D1P2','D1P3','D1P4','D1P5','D2P1','D2P2','D2P3','D2P4'})};
            end
            switch type
                case 'D1P2'
                    gInt.point = [0.577350269189626,0.577350269189626];
                    gInt.coef = [1,1];
                case 'D1P3'
                    gInt.point = [-0.774596669241483,0,0.774596669241483]; 
                    gInt.coef = [0.555555555555556,0.888888888888889,0.555555555555556];
                case 'D1P4'
                    gInt.point = [-0.861136,0.861136,-0.339981,0.339981];
                    gInt.coef = [0.347855,0.347855,0.652145,0.652145];
                case 'D1P5'
                    gInt.point = [0.906180,-0.906180,0.538469,-0.538469,0]; 
                    gInt.coef = [0.236927,0.236927,0.478629,0.478629,0.568889];
                case 'D2P1'
                    gInt.point = [1/3;1/3]; 
                    gInt.coef = 1/2;
                case 'D2P2'
                    gInt.point = [1/3,1/3;2/15,2/15;11/15,2/15;2/15,11/15]';
                    gInt.coef = 1/2*[-27/48,25/48,25/48,25/48];
                case 'D2P3'
                    gInt.point = [1/3,1/3;1/2,0;0,1/2;1/2,1/2;0,0;1,0;0,1]'; 
                    gInt.coef = 1/2*[27/60,8/60,8/60,8/60,3/60,3/60,3/60];
                case 'D2P4'
                    a1 = 0.05971587; a2 = 0.79742699; b1 = 0.47014206; b2 = 0.10128651;
                    gInt.point = [1/3,1/3; b1,b1; a1,b1; b1,a1; b2,b2; a2,b2; b2,a2]'; 
                    gInt.coef = 1/2*[0.225,0.13239415,0.13239415,0.13239415,0.12593918,0.12593918,0.12593918];
            end
        end
        function disp(gInt)
           fprintf('point: \n'), disp(gInt.point);
           fprintf('coef: \n'), disp(gInt.coef);
        end
        function tf = set.coef(tf,coef)
            if ~isvector(coef)
                error('GInt: "coef" must be a vector.');
            end
            tf.coef = coef;
        end
        function dim = get.dim(gInt)
            dim = size(gInt.point,1);
        end
        function nPoint = get.nPoint(gInt)
            nPoint = size(gInt.point,2);
        end
        function val = getVal(gInt,f)
            val = 0;
            for iPoint = 1:gInt.nPoint
                val = val + gInt.coef(iPoint)*f(gInt.point(:,iPoint));
            end
        end
        function check(gInt)
            if length(gInt.coef) ~= size(gInt.point,2)
                error('GInt: the number of variables is inconsistent - "coef", "point".');
            end
        end
    end
end

