classdef MMF14 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF14 test function
    
    properties
        POS; % Pareto optimal set
    end
    methods
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = 3; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            M = obj.M;
            g = 2 - (sin(2*pi.*X(:,end))).^2;
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
            % PF is a unit sphere (approx) when g=1 (sin^2 = 1)
            % sin(2*pi*z) = 1 => 2*pi*z = pi/2 => z=0.25 or 0.75
            % Generate points on a unit sphere
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            
            % Sample POS
            n = ceil(sqrt(N));
            [x,y] = meshgrid(linspace(0,1,n));
            x1 = x(:); x2 = y(:);
            obj.POS = [x1, x2, repmat(0.25, length(x1), 1); x1, x2, repmat(0.75, length(x1), 1)];
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        function score = CalMetric(obj,metName,Population)
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end