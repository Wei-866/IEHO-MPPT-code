function [duty,iterations,Pos] = SCSO_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num;
if isempty(num)
    num=6; % population quantity
end
if isempty(p)
    p=zeros(1,num);
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
        p(u) = vpv*ipv;
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
    [m,j]=sort(p,'descend'); 
    dbest=dc(j(1));
    dc=SCSO_UpdateDuty(dc',num,iteration,iter_max,dbest); 
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

%% SCSO-MPPT
function D=SCSO_UpdateDuty(dc,num,iteration,iter_max,dbest)
D=zeros(1,num);

p_teta=[1:360];
S=2;                                    %%% S is maximum Sensitivity range
rg=S-((S)*iteration/(iter_max));                %%%% guides R
for i=1:num
    r=rand*rg;
    R=((2*rg)*rand)-rg;                 %%%%   controls to transtion phases
    teta=RouletteWheelSelection(p_teta);
    if((-1<=R)&&(R<=1))              %%%% R value is between -1 and 1
        Rand_position=abs(rand*dbest-dc(i));
        D(i)=dbest-r*Rand_position*cos(teta);
    else
        cp=floor(num*rand()+1);
        CandidatePosition =dc(cp);
        D(i)=r*(CandidatePosition-rand*dc(i));
    end
end


for i=1:num
    dup_uplim=(D(i)-1)<0;
    dup_lowlim=D(i)>0;
    D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
end

end

function j=RouletteWheelSelection(P) 
    r=rand; 
    s=sum(P);
    P=P./s;
    C=cumsum(P); 
    j=find(r<=C,1,'first'); 
end