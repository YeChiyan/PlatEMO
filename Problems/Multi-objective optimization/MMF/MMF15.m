classdef MMF15 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF15 test function
    
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
            g = 2 - exp(-2*log10(2).*((X(:,end)-0.1)/0.8).^2).*(sin(num_of_peak*pi.*X(:,end))).^2;
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            
            % Sample POS
            n = ceil(sqrt(N));
            [x1,x2] = meshgrid(linspace(0,1,n));
            obj.POS = [x1(:), x2(:), repmat(0.25, length(x1(:)), 1); x1(:), x2(:), repmat(0.75, length(x1(:)), 1)];
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