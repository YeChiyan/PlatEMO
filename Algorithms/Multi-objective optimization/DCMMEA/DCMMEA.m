classdef DCMMEA < ALGORITHM
% <multi/multimodal> <real/binary/permutation><multimodal> <constrained>
% Balancing exploration and exploitation in dynamic constrained multimodal multi-objective co-evolutionary algorithm
% Swarm and Evolutionary Computation, 89, 101652, 2024
% Guoqing Li (li241700@126.com)

    methods
        function main(Algorithm,Problem)

            %% Generate random population
            Population1 = Problem.Initialization(Problem.N/2);
            Population2 = Problem.Initialization(Problem.N/2);
            
%                 x1=std(Population1.cons);
%                 y1=std(Population1.objs);
%                 H=1;
%                 epsilon_x1=H* sum(x1)/Problem.M;
%                 epsilon_y1=H* sum(y1)/Problem.M;
                [epsilon_x1,epsilon_y1]=parameters(Population1,Population2,Problem);
            
            
            %% Calculate fitness of populations
            [D_Dec,D_Pop,~]    = CalFitness(Population1.objs,Population1.cons,Population1.decs,epsilon_x1,Problem);
            [D,P,~]    = CalFitness2(Population2.objs,Population2.decs,epsilon_y1,Problem);
            
            len1=length(Population1);
            len2=length(Population2);
            
           

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                           
                MatingPool1 = TournamentSelection(2,len1,D_Dec);
                MatingPool2 = TournamentSelection(2,len2,D);
                Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                [Population1,D_Dec,D_Pop] = EnvironmentalSelection([Population1,Offspring1,Offspring2],len1,epsilon_x1,Problem);
                [Population2,D,P] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],len2,epsilon_y1,Problem);
                
                
                len1=ceil(Problem.N/2*(1+sin((Problem.FE/Problem.maxFE)*0.5*pi)));
                len2=ceil(Problem.N/2*(1-sin((Problem.FE/Problem.maxFE)*0.5*pi)));
                if len2<5
                    len2=5;
                end
                [epsilon_x1,epsilon_y1]=parameters(Population1,Population2,Problem);
%                 num=Problem.FE/Problem.maxFE;
%                 if num<0.5
%                     H=1-0.5*exp((num-0.5)/0.1);
% %                     H=0.5*exp(((Problem.FE/Problem.maxFE)-1)/0.5);
%                 else
%                     H=0.5*exp(-(num-0.5)/0.1);
% %                     H=1-0.5*exp(((Problem.FE/Problem.maxFE)-1)/0.5);
%                 end
%                 x1=std(Population1.cons);
%                 y1=std(Population2.objs);
%                 epsilon_x1=0.5*H* sum(x1)/Problem.M;
%                 epsilon_y1=0.5*H* sum(y1)/Problem.M;
               
                                
%                 epsilon_k =H* sum(min(Population1.objs,[],1),2)/Problem.M;
            [D_Dec,D_Pop,~]    = CalFitness(Population1.objs,Population1.cons,Population1.decs,epsilon_x1,Problem);
            [D,P,~]    = CalFitness2(Population2.objs,Population2.decs,epsilon_y1,Problem);
            end
        end
    end
end