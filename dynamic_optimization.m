%input: N - number of villages/baranggays in the considered town
%       M - number of possible vaccination sites
%       Pj P(j,2) - population of the jth village/baranggay
%       Tp - total population of the considered town
%       Tc - total number of COVID-19 cases
%       Cj P(j,1) - number of COVID-19 cases in the jth village/baranggay
%       D - distance profile
%       L - number of optimal vaccination site
%       Vmax - maximum daily vaccination per site
%       Vdays - duration of each vaccination schedule
%output: Names of vaccination site for each schedule

clear; 
N=42;
D=readmatrix('distances.csv','Range','B2:AQ66');
P_store=readmatrix('Barangay_Centers_Table.xlsx','Range','A2:B43');
names = readtable('AllData.xlsx','Sheet','Sites','Range','C2:C66','ReadVariableNames',false);
names = table2cell(names);
Tc=579;
Tp=125252;
L=2;
intcon=1:L;
lb=ones(1,L);
ub=65*ones(1,L);
NN=10;
NG = 100*L;
NP = 200*L;
Vmax=200;
Vdays=30;
P=P_store(:,2);

for k=1:floor(0.7*Tp/(L*Vmax*Vdays))
    
parfor i=1:NN
[xopt,f_xopt]= ga(@(x) obfn(x,D,P_store,Tc,Tp,L,P),L,[],[],[],[],lb,ub,[],intcon);%);,options);
xbest(i,:)=xopt;
xcost(i)=f_xopt;
end

  
[val,~]=min(xcost);
ind=find(xcost==val);

soln_temp=unique(sort(xbest(ind,:),2),'rows');
soln(k,:)=soln_temp(1,:);
bestsites = names(soln)

D_assign=zeros(N,1);
ind_store1=[];
ind_store2=[];
for ii=1:N
    if D(soln(k,1),ii)<D(soln(k,2),ii)
        D_assign(ii)=soln(k,1);
        ind_store1=[ind_store1;ii];
    else
        D_assign(ii)=soln(k,2);
        ind_store2=[ind_store2;ii];
    end
end

%sampling
[temp,ind_save1]=sort(D(soln(k,1),ind_store1));
ind_store1=ind_store1(ind_save1)';
Dist1=D(soln(k,1),ind_store1)+1;
Prod1=(1./Dist1);
for_prob1=(Prod1)/sum(Prod1);
for_prob1=[0 cumsum(for_prob1)];

[temp2,ind_save2]=sort(D(soln(k,2),ind_store2));
ind_store2=ind_store2(ind_save2)';
Dist2=D(soln(k,2),ind_store2)+1;
Prod2=(1./Dist2);
for_prob2=(Prod2)/sum(Prod2);
for_prob2=[0 cumsum(for_prob2)];
i2=0;


while i2<Vmax*Vdays
    r1=rand();   
    barangay_selector1=max(find(r1>for_prob1));
    if (P(ind_store1(barangay_selector1))<0.3*P_store(ind_store1(barangay_selector1),2))
        dummy1=find(P(ind_store1)>=0.3*P_store(ind_store1,2));
        if isempty(dummy1)==1
        dummy2=find(P(ind_store2)>=0.3*P_store(ind_store2),2);
        if isempty(dummy2)==1
        i2=Vmax*Vdays;
        disp('70% of population is vaccinated');
        else
        barangay_selector2=dummy2(randi(length(dummy2)));
        P(ind_store2(barangay_selector2))=P(ind_store2(barangay_selector2))-1;
        end
        
        else
        barangay_selector1=dummy1(randi(length(dummy1)));
        P(ind_store1(barangay_selector1))=P(ind_store1(barangay_selector1))-1; 
        end        
    end
    P(ind_store1(barangay_selector1))=P(ind_store1(barangay_selector1))-1;
      
    
    r2=rand();   
    barangay_selector2=max(find(r2>for_prob2));
    if (P(ind_store2(barangay_selector2))<0.3*P_store(ind_store2(barangay_selector2),2))
        dummy2=find(P(ind_store2)>=0.3*P_store(ind_store2,2));
        if isempty(dummy2)==1
        dummy1=find(P(ind_store1)>=0.3*P_store(ind_store1,2));
        if isempty(dummy1)==1
        i2=Vmax*Vdays;
        else
        barangay_selector1=dummy1(randi(length(dummy1)));
        P(ind_store1(barangay_selector1))=P(ind_store1(barangay_selector1))-1;
        end
        
        else
        barangay_selector2=dummy2(randi(length(dummy2)));   
        P(ind_store2(barangay_selector2))=P(ind_store2(barangay_selector2))-1;
        end  
    end
    P(ind_store2(barangay_selector2))=P(ind_store2(barangay_selector2))-1;
    i2=i2+1;
    
end
UNVAX(k,:)=(P./P_store(:,2)*100)';
end

function J = obfn(x,D,P,Tc,Tp,L,Unvax)
    J=0;
    x=ceil(x);
    for j=1:42
        dum=zeros(L,1);
        for i=1:L
            dum(i) = D(x(i),j);
        end
        [om,~] = sort(dum);
        out = om(1);
        if Unvax(j)<0.3*P(j,2)
            eta=0;
        else
            eta=1;
        end
        J = J + eta*((P(j,1)/Tc) + (P(j,2)/Tp))*out;
    end
end
