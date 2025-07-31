function [duty,iterations,Pos] = PSO_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num   pbest  v;
if isempty(num)
    num=6; % population quantity
end
if isempty(p)
    p=zeros(1,num);
    Pold=1;
    pbest=zeros(1,num);
    counter=1;
    u=1;
    iteration=0;
    iter_max=15;
    v=zeros(1,num);
end
if isempty(dc)
    %% Initialize the population
    dc=zeros(1,num);
    dc = linspace(0, 1, 8);   
    dc = dc(2:end-1);
    dbest=0;
end

Pos=dc';
iterations=iteration;
if iterations<=iter_max
    if(counter>=1 && counter<=1500)
        duty=dc(u);
        counter=counter+1;
        return;
    end

    if(u>=1 && u<=num)
        if((vpv*ipv)>p(u))
            p(u) = vpv*ipv;
            pbest(u)=dc(u);
        end
    end

    u=u+1;
    if(u<num+1)
        duty=dc(u);
        counter=1;
        return;
    end
    u=1;
    counter=1;
    iteration=iteration+1;
    [~,i]=max(p);
    dbest=pbest(i);
    dc=PSO_UpdateDuty(dc',v',pbest,dbest,num);
    duty=dc(u);
    Pold=vpv*ipv;
    return;
else
    duty=dbest;
    P=vpv*ipv;
    if abs(P-Pold)>0.05*Pold  
        iteration=0;
        dc = linspace(0, 1, 8);   
        dc = dc(2:end-1);
    end
    Pold=P;
    return;
end
end

%% PSO-MPPT
function D=PSO_UpdateDuty(dc,v,pbest,dbest,num)
D=zeros(1,num);
for i=1:num
    w=0.1;
    c1=1.2;
    c2=1.2;
    v(i) = (w*v(i))+(c1*rand(1)*(pbest(i)-dc(i)))+(c2*rand(1)*(dbest-dc(i)));
    D(i)=dc(i)+v(i);
end

for i=1:num
    dup_uplim=(D(i)-1)<0;
    dup_lowlim=D(i)>0;
    D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
end

end



    
