classdef Omni_test < PROBLEM
    % <multi> <real> <multimodal>
    % Omni_test test function
    
    properties
        POS; % Pareto optimal set
    end
    methods
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 3; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D)*6;
            obj.encoding = ones(1,obj.D);
        end
        function PopObj = CalObj(obj,X)
            PopObj(:,1) = sum(sin(pi*X), 2);
            PopObj(:,2) = sum(cos(pi*X), 2);
        end
        function R = GetOptimum(obj,N)
            % Numerical estimation for PF
            % For Omni-test, the PF is a bit specific.
            % We'll just provide a placeholder or sample points.
            obj.POS = []; % Complex PS structure
            R = [0, 0]; % Placeholder
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end
