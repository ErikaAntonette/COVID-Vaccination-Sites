%input: N - number of villages/baranggays in the considered town
%       M - number of possible vaccination sites
%       Pj P(j,2) - population of the jth village/baranggay
%       Tp - total population of the considered town
%       Tc - total number of COVID-19 cases
%       Cj P(j,1) - number of COVID-19 cases in the jth village/baranggay
%       D - distance profile
%       L - number of optimal vaccination site
%output: x - index array of vaccination site
clear; clc;
D=readmatrix('distances.csv','Range','B2:AQ66');
P=readmatrix('Barangay_Centers_Table.xlsx','Range','A2:B43');
names = readtable('AllData.xlsx','Sheet','Sites','Range','C2:C66','ReadVariableNames',false);
names = table2cell(names);
Tc=579;
Tp=125252;
L=4;
intcon=1:L;
lb=ones(1,L);
ub=65*ones(1,L);
NN=20*L*L;
parfor i=1:NN
[xopt,f_xopt]= ga(@(x) obfn(x,D,P,Tc,Tp,L),L,[],[],[],[],lb,ub,[],intcon);
xbest(i,:)=xopt;
xcost(i)=f_xopt;
end
[val,~]=min(xcost)
ind=find(xcost==val);
%disp(xbest(ind,:));

soln=unique(sort(xbest(ind,:),2),'rows')
bestsites = names(soln)
%disp([soln bla]);
%fun=@(x) @(x) obfn(x,D,P,Tc,Tp,L);
%xopt=particleswarm(@(x) obfn(x,D,P,Tc,Tp,L),length(lb),lb,ub);
%disp(ceil(xopt));

% J=obfn([65,64],D,P,Tc,Tp)

function J = obfn(x,D,P,Tc,Tp,L)
    J=0;
    x=ceil(x);
    for j=1:42
        dum=zeros(L,1);
        for i=1:L
            dum(i) = D(x(i),j);
        end
        [om,~] = sort(dum);
        out = om(1);
        J = J + ((P(j,1)/Tc) + (P(j,2)/Tp))*out;
    end
end
