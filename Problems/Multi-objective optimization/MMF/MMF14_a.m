classdef MMF14_a < PROBLEM
    % <multi> <real> <multimodal>
    % MMF14_a test function
    
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
            num_of_peak = 2;
            x_g = X(:,end) - 0.5*sin(pi*X(:,end-1));
            g = 2 - (sin(num_of_peak*pi.*(x_g+1/(2*num_of_peak)))).^2;
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            
            % Sample POS: z = 0.5*sin(pi*y) + {0, 0.5}
            n = ceil(sqrt(N));
            [x,y] = meshgrid(linspace(0,1,n));
            x1 = x(:); y1 = y(:);
            % Subregion 1
            z1 = 0.5*sin(pi*y1);
            % Subregion 2
            z2 = 0.5 + 0.5*sin(pi*y1);
            
            obj.POS = [x1, y1, z1; x1, y1, z2];
            % Ensure bounds
            obj.POS = max(0, min(1, obj.POS));
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