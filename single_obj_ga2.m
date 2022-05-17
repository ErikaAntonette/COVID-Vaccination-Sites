clear; clc;
D=readmatrix('distances.csv','Range','B2:AQ66'); %distance profile
P=readmatrix('Barangay_Centers_Table.xlsx','Range','A2:B43'); %confirmed cases and population size
names = readtable('AllData.xlsx','Sheet','Sites','Range','C2:C66','ReadVariableNames',false);
names = table2cell(names);
Tc=579; %total no. of confirmed cases
Tp=125252; %total population
TB = 42; %total no. of barangays
TV = 65; %total no. of vaccination sites
L=4; %number of vaccination sites

intcon=1; %the variable for which integer constraint is applied
st=10000; %stopping criterion
combs = nchoosek(1:TV,L); %vector of combinations of vaccination sites
N=length(combs(:,1)); %number of combinations
lb=1; %lowerbound of index
ub=N; %upperbound of index

w1 = zeros(TB,1);
for k = 1 : TB
   w1(k) = P(k,1)/Tc; %constant for confirmed cases
end

w2 = zeros(TB,1);
for k = 1 : TB
   w2(k) = P(k,2)/Tp; %constant for population
end

%stopping criterion for algorithms
%for ga
options = optimoptions('ga','MaxGenerations',st); 
%for surrogate
options1 = optimoptions('surrogateopt','MaxFunctionEvaluations',st,'Display','off','PlotFcn',[]);

%calling ga
NN=20*L*L;
parfor i=1:NN
[xopt,cxopt] = ga(@(x) obfnn(x,D,w1,w2,combs),1,[],[],[],[],lb,ub,[],intcon,options);
min_x(i,:)=xopt;
min_cost(i)=cxopt;
%disp(combs(xopt,:))
%disp(cxopt)
end
[min_C, ~]=min(min_cost);
ind=find(min_cost==min_C);
soln=unique(sort(combs(min_x(ind,:),:),2),'rows')
bestsites = names(soln)
disp(min_C);

%calling surrogate 
% [xopt1,cxopt1] = surrogateopt(@(x) obfnn(x,D,w1,w2,combs),lb,ub,intcon,options1);
% disp(combs(xopt1,:))
% disp(cxopt1)

%cost function of 'true' solution
%J=obfnn(84,D,w1,w2,combs)

%cost function
function J = obfnn(x,D,w1,w2,combs)
    x1=combs(x,:)';
    J = (w1+w2)'*(min(D(x1,:)))';
end
