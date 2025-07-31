function [duty,iterations,Pos] = GWO_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num pbest alpha alpha_pos Beta_pos Beta_score Delta_score Delta_pos v;
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
if (isempty(Beta_pos))
    alpha_pos=zeros(1,1);
    Beta_pos=zeros(1,1);
    Delta_pos=zeros(1,1);
    Beta_score=0;
    Delta_score=0;
    alpha=0;
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
        p(u)=vpv*ipv;
        if p(u)>alpha
            alpha=p(u); % Update alpha
            alpha_pos=dc(u);
        end
        if p(u)<alpha && p(u)>Beta_score
            Beta_score=p(u); % Update beta
            Beta_pos=dc(u);
        end

        if p(u)<alpha && p(u)<Beta_score && p(u)>Delta_score
            Delta_score=p(u); % Update delta
            Delta_pos=dc(u);
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
    dc=GWO_UpdateDuty(iteration,alpha_pos,dc',Beta_pos,Delta_pos,num);
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

%% GWO-MPPT
function D=GWO_UpdateDuty(l,alpha_pos,D_cur,Beta_pos,Delta_pos,num)
D=zeros(1,num);
a=2-2*((l-1)/14);
for kk=1:num
    r1=rand(); 
    r2=rand();
    A1=2*a*r1-a;
    C1=2*r2;
    D_alpha=abs(C1*alpha_pos-D_cur(kk));
    X1=alpha_pos-A1*D_alpha;
    r1=rand();
    r2=rand();
    A2=2*a*r1-a;
    C2=2*r2;
    D_beta=abs(C2*Beta_pos-D_cur(kk));
    X2=Beta_pos-A2*D_beta;
    r1=rand();
    r2=rand();

    A3=2*a*r1-a;
    C3=2*r2;

    D_delta=abs(C3*Delta_pos-D_cur(kk));
    X3=Delta_pos-A3*D_delta;

    X4=(X1+X2+X3)/3;
    D(kk)=X4;
end

for i=1:num
    dup_uplim=(D(i)-1)<0;
    dup_lowlim=D(i)>0;
    D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
end

end