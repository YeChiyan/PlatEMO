classdef MMF16_l3 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF16_l3 test function
    
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
            [N, ~] = size(X);
            g = zeros(N, 1);
            num_of_g_peak = 2;
            num_of_l_peak = 2;
            
            z = X(:, end);
            low_idx = z >= 0 & z < 0.5;
            high_idx = z >= 0.5 & z <= 1;
            
            if any(low_idx)
                g(low_idx) = 2 - (sin(2*num_of_g_peak*pi.*z(low_idx))).^2;
            end
            if any(high_idx)
                g(high_idx) = 2 - exp(-2*log10(2).*((z(high_idx)-0.1)/0.8).^2).*(sin(2*num_of_l_peak*pi.*z(high_idx))).^2;
            end
            
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(X(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(X(:,M-1:-1:1)*pi/2)];
        end
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            
            % Sample POS: z = 0.125 and z = 0.375
            n = ceil(sqrt(N));
            [x,y] = meshgrid(linspace(0,1,n));
            x1 = x(:); y1 = y(:);
            obj.POS = [x1, y1, repmat(0.125, length(x1), 1); x1, y1, repmat(0.375, length(x1), 1)];
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
