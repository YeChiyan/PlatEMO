classdef SYM_PART_rotated < PROBLEM
    % <multi> <real> <multimodal>
    % SYM-PART-rotated test function
    
    %------------------------------- Reference --------------------------------
    % Y. Liu, G. G. Yen, and D. Gong, A multi-modal multi-objective
    % evolutionary algorithm using two-archive and recombination strategies,
    % IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
    %--------------------------------------------------------------------------
    
    properties
        POS; % Pareto optimal set
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = [-20, -20];
            obj.upper    = [20, 20];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(~,X)
            a = 1; b = 10; c = 8;
            w = pi/4; % rotation angle
            [N, ~] = size(X);
            PopObj = zeros(N, 2);
            for i = 1:N
                % Apply rotation
                xx1 = cos(w)*X(i,1)-sin(w)*X(i,2);
                xx2 = sin(w)*X(i,1)+cos(w)*X(i,2);
                
                % Cell-shifting logic
                temp_t1 = sign(xx1)*ceil((abs(xx1)-(a+c/2))/(2*a+c));
                temp_t2 = sign(xx2)*ceil((abs(xx2)-b/2)/b);
                t1 = sign(temp_t1)*min(abs(temp_t1),1);
                t2 = sign(temp_t2)*min(abs(temp_t2),1);
                x1 = xx1-t1*(c+2*a);
                x2 = xx2-t2*b;
                
                PopObj(i,1) = (x1+a)^2 + x2^2;
                PopObj(i,2) = (x1-a)^2 + x2^2;
            end
        end
        %% Generate Pareto optimal solutions
        function R = GetOptimum(obj,N)
            % Generate points in Pareto optimal set
            % Rotated 9 Pareto sets
            a = 1; b = 10; c = 8;
            npos = 100;
            x1_base = linspace(-a, a, npos)';
            x2_base = zeros(npos, 1);
            t1_vals = [-1, 0, 1] * (c+2*a);
            t2_vals = [-1, 0, 1] * b;
            pos_rotated = [];
            for t1 = t1_vals
                for t2 = t2_vals
                    tx = x1_base + t1;
                    ty = x2_base + t2;
                    % [X1; X2] = [cos w, sin w; -sin w, cos w] * [xx1; xx2]
                    X1 = cos(pi/4)*tx + sin(pi/4)*ty;
                    X2 = -sin(pi/4)*tx + cos(pi/4)*ty;
                    pos_rotated = [pos_rotated; X1, X2];
                end
            end
            obj.POS = pos_rotated;
            
            % Base PF (unaffected by rotation)
            x = linspace(-1, 1, N)';
            R(:,1) = (x+1).^2;
            R(:,2) = (x-1).^2;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            % Generate POS for IGDX calculation
            % Rotated 9 Pareto sets
            if isempty(obj.POS)
                a = 1; b = 10; c = 8;
                N = 100;
                x1_base = linspace(-a, a, N)';
                x2_base = zeros(N, 1);
                
                t1_vals = [-1, 0, 1] * (c+2*a);
                t2_vals = [-1, 0, 1] * b;
                
                pos_rotated = [];
                for t1 = t1_vals
                    for t2 = t2_vals
                        tx = x1_base + t1;
                        ty = x2_base + t2;
                        % Rotate back to find original coordinates
                        % x' = x cos w - y sin w
                        % y' = x sin w + y cos w
                        % w is pi/4 clockwise? No, original was counter-clockwise rotation of coordinates?
                        % Original code applied rotation to X to get xx1,xx2 which are used in local calculation.
                        % So the PS occurs when xx2=0 and -a<=xx1<=a (plus shifts).
                        % To get X from xx, we use the inverse rotation:
                        % [X1; X2] = [cos w, sin w; -sin w, cos w] * [xx1; xx2]
                        X1 = cos(pi/4)*tx + sin(pi/4)*ty;
                        X2 = -sin(pi/4)*tx + cos(pi/4)*ty;
                        pos_rotated = [pos_rotated; X1, X2];
                    end
                end
                obj.POS = pos_rotated;
            end
            
            switch metName
                case 'IGDX'
                    score = feval(metName,Population,obj.POS);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
    end
end