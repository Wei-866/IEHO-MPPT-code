function [duty,iterations,Pos] = SSA_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num ;
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
    iteration=iteration+1;
    [m,j]=sort(p,'descend'); 
    dbest=dc(j(1));
    dc=SSA_UpdateDuty(dbest,dc',iteration,iter_max,num); 
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

%% SSA-MPPT
function D=SSA_UpdateDuty(dbest,dc,iteration,iter_max,num)

D=zeros(1,num);
c1 = 2*exp(-(4*iteration/iter_max)^2); % Eq. (3.2) in the paper

for i=1:num

    if i<=num/2
        c2=rand();
        c3=rand();
        %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
        if c3<0.5
            D(i)=dbest+c1*c2;
        else
            D(i)=dbest-c1*c2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif i>num/2 && i<num+1
        point1=dc(i-1);
        point2=dc(i);
        D(i)=(point2+point1)/2; % % Eq. (3.4) in the paper
    end
end


for i=1:num
    dup_uplim=(D(i)-1)<0;
    dup_lowlim=D(i)>0;
    D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
end

end