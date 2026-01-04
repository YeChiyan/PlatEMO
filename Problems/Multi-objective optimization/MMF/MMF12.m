classdef MMF12 < PROBLEM
    % <multi> <real> <multimodal>
    % MMF12 test function
    
    properties
        POS; % Pareto optimal set
    end
    methods
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [0, 0];
            obj.upper    = [1, 1];
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(~,X)
            PopObj(:,1) = X(:,1);
            num_of_peak = 2;
            Alfa = 2;
            q = 4;
            g = 2 - ((sin(num_of_peak*pi.*X(:,2))).^6).*(exp(-2*log10(2).*((X(:,2)-0.1)/0.8).^2));
            h = 1 - (PopObj(:,1)./g).^Alfa - (PopObj(:,1)./g).*sin(2*pi*q*PopObj(:,1));
            PopObj(:,2) = g.*h;
        end
        function R = GetOptimum(obj,N)
            % Typical PF shape for this kind of h
            % We'll use a numerical or high-density sampling for g=1
            x = linspace(0, 1, N)';
            % PS at y=0.25 and y=0.75
            obj.POS = [x, repmat(0.25, N, 1); x, repmat(0.75, N, 1)];
            
            % PF estimation
            g = 1;
            R(:,1) = x;
            R(:,2) = g * (1 - (x./g).^2 - (x./g).*sin(2*pi*4*x));
            % Filter non-dominated
            [FrontNo, ~] = NDSort(R, 1);
            R = R(FrontNo == 1, :);
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