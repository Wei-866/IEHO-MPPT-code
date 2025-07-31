function [duty,iterations,Pos] = RIME_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num dc_old p_old;
if isempty(num)
    num=6; % population quantity
end
if isempty(p)
    p=zeros(1,num);
    p_old=zeros(1,num);
    Pold=1;
    counter=1;
    u=1;
    iteration=0;
    iter_max=15;
   
end
if isempty(dc)

    %% Initialize the population
    dc=zeros(1,num);
    dc = linspace(0, 1, 8);  
    dc = dc(2:end-1);
    dbest=0;
    dc_old=dc;
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
         p(u)=vpv*ipv; 
    end

    u=u+1;
    if(u<num+1)
        duty=dc(u);
        counter=1;
        return;
    end
    u=1;
    counter=1;
    if iteration>0
        NewFitness=cat(2,p,p_old);
        NewPopulation=cat(2,dc,dc_old);
        [sorted_NewFitness,sorted_indexes]=sort(NewFitness,'descend');
        for i=1:num
            dc(i)=NewPopulation(sorted_indexes(i));
            p(i)= sorted_NewFitness(i);
        end
    end
    iteration=iteration+1;
    p_old=p;
    dc_old=dc;
    [m,j]=sort(p,'descend'); 
    dbest=dc(j(1));
    dc=RIME_UpdateDuty(dbest,dc',iteration,iter_max,num,p); 
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

%% RIME-MPPT
function D=RIME_UpdateDuty(dbest,dc,iteration,iter_max,num,p)
coder.extrinsic('normr');
normalized_rime_rates=zeros(1,num);
D=zeros(1,num);
W = 5;
RimeFactor = (rand-0.5)*2*cos((pi*iteration/(iter_max/10)))*(1-round(iteration*W/iter_max)/W);%Parameters of Eq.(3),(4),(5)
E =(iteration/iter_max)^0.5;%Eq.(6)
newRimepop = dc;%Recording new populations
normalized_rime_rates=normr(p);%Parameters of Eq.(7)
for i=1:num
    r1=rand();
    if r1< E
        newRimepop(i)=dbest+RimeFactor*rand();%Eq.(3)
    end
    %Hard-rime puncture mechanism
    r2=rand();
    if r2<normalized_rime_rates(i)
        newRimepop(i)=dbest;%Eq.(7)
    end
    dup_uplim=(newRimepop(i)-1)<0;
    dup_lowlim=newRimepop(i)>0;
    D(i)=newRimepop(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
end

end