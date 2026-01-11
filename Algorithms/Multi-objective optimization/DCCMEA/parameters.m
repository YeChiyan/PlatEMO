function [detal1,detal2]=parameters(Population1,Population2,Problem)
N=length(Population1);
Con=Population1.cons;
Con = sum(max(0,Con),2);
Con=sort(Con);
ConX1=Con(1:N-1);
ConX2=Con(2:N);
ConX=ConX2-ConX1;
detal1=mean(ConX);

N2=length(Population2);
Obj=Population2.objs;
for i=1:Problem.M
    O=Obj(:,i);
    O=sort(O);
    O1=O(1:N2-1,1);
    O2=O(2:N2,1);
    O=O2-O1;
    m(i)=mean(O);    
end
detal2=mean(m);

num=Problem.FE/Problem.maxFE;
if num<0.5
      H=1-0.5*exp((num-0.5)/0.1);
else
      H=0.5*exp(-(num-0.5)/0.1);
end

detal1=H*detal1;
detal2=H*detal2;


end