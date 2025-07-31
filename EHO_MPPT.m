function [duty,iterations,Pos] = EHO_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num num_of dc_old p_old  ;

if isempty(num)
    num=6; % population quantity
    MalesRate=0.2;
    num_of = round(num*MalesRate);
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
    dc=EHO_UpdateDuty(dc',num,p,num_of); 
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

%% EHO-MPPT
function D=EHO_UpdateDuty(dc,num,p,num_of)

D=zeros(1,num);
TransposeFitness=zeros(1,num_of);

%Sort the ELK positions
    [sorted_ELKS_fitness,sorted_indexes]=sort(p,'descend');

    % Make a copy of population
    NewElkHerd = dc;
    NewElkHerdFitness = p;
    
    dbest=dc(sorted_indexes(1),:); % ELK with best position
    dbest_fit=sorted_ELKS_fitness(1); % the fitness of best position

    % Number of females for each male
    for i=1:num_of
        TransposeFitness(i)=sorted_ELKS_fitness(i);
    end

    Familes=zeros(1,num);
    for i=(num_of+1):num
        FemaleIndex=sorted_indexes(i); % index of female
        randNumber = rand;
        MaleIndex = 0;
        sum_fitness=0;
               
        for j=1:num_of
            sum_fitness = sum_fitness + (TransposeFitness(j)/sum(TransposeFitness));
            if (sum_fitness > randNumber)
                MaleIndex =  j;
                break;
            end
        end
        Familes(FemaleIndex)=sorted_indexes(MaleIndex);
    end
    
     
    %===================== Reproduction
    for i=1:num
        %Male
        if(Familes(i)==0)
            h=fix(rand*num)+1;
            D(i)=dc(i)+rand*(dc(h)-dc(i));    
        else
            h=fix(rand*num)+1;
            MaleIndex=Familes(i);
            hh = randperm(size(find(Familes==MaleIndex),2));
            h=round(1+(size(hh,2)-1)*rand);
            rd = -2+4*rand;
            D(i)=dc(i)+(dc(Familes(i))-dc(i))+rd*(dc(h)-dc(i));
        end
    end
    for i=1:num
        dup_uplim=(D(i)-1)<0;
        dup_lowlim=D(i)>0;
        D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
    end

end