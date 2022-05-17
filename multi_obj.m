clear; clc; format longG; close all;

for L = 2:4
tic
eps = 3000; %distance in meters for coverage obj fcn

TI = 579; %total infected
TP = 125252; %total population
TB = 42; %total no. of barangays
TV = 65; %total no. of vaccination sites

combs = nchoosek(1:TV,L);

I = readmatrix('AllData.xlsx','Sheet','Barangays','Range','A2:A43');
P = readmatrix('AllData.xlsx','Sheet','Barangays','Range','B2:B43');

names = readtable('AllData.xlsx','Sheet','Sites','Range','C2:C66','ReadVariableNames',false);
names = table2cell(names);

w1 = zeros(TB,1);
for k = 1 : TB
   w1(k) = I(k)/TI;    
end

Distance = readmatrix('AllData.xlsx','Sheet','Distances','Range','B2:AQ66');
Cost1 = zeros(size(combs,1),1);

for j = 1:size(combs,1)
    if L == 1
     Cost1(j) = sum(w1'.*(Distance(combs(j,:)',:)));
    else
        Cost1(j) = sum(w1'.*min(Distance(combs(j,:)',:)));
    end
end

[X1,~]=mink(Cost1,1); %% change 2nd entry to get top N sites
Y1=find(Cost1==X1);
bestcombs1 = combs(Y1,:);
bestsites1 = names(bestcombs1);

Cost2 = zeros(size(combs,1),1);

%%
for j = 1:size(combs,1)
     Index = Distance <= eps;
     alpha = Index(combs(j,:),:);
     if L ~= 1
     alpha = min(sum(alpha),1);
     end
     Cost2(j) = -1*(alpha*P);
end
%%
[X2,~]=mink(Cost2,1); %% change 2nd entry to get top N sites
Y2=find(Cost1==X2);
bestcombs2 = combs(Y2,:);
bestsites2 = names(bestcombs2);

figure;
%scatter(Cost2,Cost1,25)
title(['Optimal Vacc. Sites in San Juan (L = ', num2str(L), ')'])
xlabel('f2')
ylabel('f1')
%xlim([min(Cost2) max(Cost2)]);
%ylim([min(Cost1) max(Cost1)]);
%%
costs = [Cost2,Cost1];
[val,I1] = sortrows([Cost2,Cost1],1);
extract = val(:,2);
bestyval = extract(1);
%%
index = 1;
for i = 2:size(extract,1)
    if extract(i) < bestyval
        bestyval = extract(i);
        index = [index; i];
    end
end

paretoset = val(index,:);
[C,IA,IC] = unique(paretoset(:,1),'last');
paretoset = paretoset(IA,:);

preindex = I1(index);
finalindex = preindex(IA);

sets = combs(ismember(costs,costs(finalindex,:),'rows'),:)
% sets = combs(finalindex,:)
pareto_optimal_set = names(sets)

hold on;
plot(paretoset(:,1),paretoset(:,2),'-mo','LineWidth',2,'MarkerSize',5)
toc
end