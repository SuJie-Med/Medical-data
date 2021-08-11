function [GilS,NSilS,S0]=ComputeLoss(G,attribute,GX,k,aa,bb)
%%
%%%%%%%%%%%%%% read the data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=GX;
nn=size(X.N,1);%%%nodes in graph
K=k;%%% the number of elements in a cluster
nr=floor(double(nn)/K);%%%%the number of clusters
n=nr*K;
R=n;
a=aa;b=bb;
ap=a*10;
attr=attribute;
G0=G;
A=G;
%%
%%%%%%%%%%%%%%% record the degree of each node %%%%%%%%%%%%
S0(1,1)=0;
t=0;  %cluster计数器
while(R>0)
if (t==nr)
    break;
end;
t=t+1; 
%% 
%%%%%%%%%%%%%%%% compute the degree of each node %%%%%%%%%%%
D=zeros(nn);
for i=1:nn
    D(i)==length(find(A(i,:)==1));
   
end;

%% 
%%%%%%%%%%%%%%%% find the maxdegree %%%%%%%%%%%%%%%%%%%%%%%%

Maxd=0;k=0;
for i=1:nn
   if(i~=S0)
     if(Maxd<D(i))
       Maxd=D(i)
       k=i;
     end;
   end;
end;
if (k==0)
  for i=1:nn
    if(i~=S0)
      k=i;
       break;
    end;    
  end;
end;
S0(t,1)=k;
%% 
%%%%%%%%%%%%%%%% cluster中其余节点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;%计数器，1:K－1
while (m<K)
%% Computer the generationg information loss
gen=X.N(S0(t,1:m));
for j=1:m
  p=S0(t,j);
  for i=1:nn
  % attibutes c0
    if  strcmp(X.X0(i),X.X0(p))
      S.c0(i,j)=0;
    else
      S.c0(i,j)=1;
     end;
  % attributes c1
    mi=floor(double(X.X1(i))/1000);
    mk=floor(double(X.X1(p))/1000);
    if  strcmp(X.X1(i),X.X1(p))
      S.c1(i,j)=0;
      elseif mi==mk
       S.c1(i,j)=0.5;
    else 
     S.c1(i,j)=1;  
    end;
  % attributes c2
    if strcmp(X.X2(i),X.X2(p))
      S.c2(i,j)=0;
    else
      S.c2(i,j)=1;
    end;
 % attributes c3
    if strcmp(X.X3(i),X.X3(p))
      S.c3(i,j)=0;
    else
      S.c3(i,j)=1;
    end;
  end;

end;
 Nloss=ones(1,nn)*1000;
for i=1:nn    
    if (i~=S0(1:t,1:m))
        geni=gen;
        geni(m+1)=X.N(i);
        Nloss(i)=double(max(geni)-min(geni))/double(max(X.N)-min(X.N));
    end;
end;

%% compute the Gloss
Gloss=ones(1,nn)*1000;
for i=1:nn
   if (i~=S0(1:t,1:m)) 
     Gloss(i)=(m*Nloss(i)+max(S.c0(i,1:m))+max(S.c1(i,1:m))+max(S.c2(i,1:m))+max(S.c3(i,1:m)))/(nn*5) 
   end;
end;
%% Computer DD
%U=n2short(A);% the least distance between any nodes
U=zeros(nn,nn);
pathnumb=zeros(nn,nn);
for i=1:nn
    for j=i+1:nn
     [U(i,j),pathnumb(i,j)]=dijkstra(A,i,j);% the least distance between any nodes
    end;
end;
U=U+U';
pathnumb=pathnumb+pathnumb';
%% cluster内结点到center的距离均值
for i=1:nn
%    if (i~=S0(1:t,1:m))
    if (i~=S0)
        sum=0;
        for dt=1:m   
           
            sum=sum+U(i,S0(t,dt))/(pathnumb(i,S0(t,dt)));           
        end;
       DD(i)=sum/m;
    end;
end;
%% q is the selected node 
q=0;Gl=10000000000  
%for i=1:n
%%%%%%??????????????????????????????????????????????????????????????????????????????????????2020.8.20
for i=1:nn
  if(i~=S0)
      if(Gl>(a*Gloss(i)+b*DD(i)))
         Gl=(a*Gloss(i)+b*DD(i)) 
         q=i;
      end;
  end;
end;
%[Gl,q]=min(a*Gloss+b*DD);
m=m+1;
S0(t,m)=q;
end;
%%
%%%%%%%%%%%%%%%% 从图中删去该族%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(S0,2);
for j=1:m
  % for i=1:nn  %%%%%%%%%%%%%%%%%%%%%???????????????????????????????????????????????????????2020.8.20
  for i=1:nn
      if (i~=S0(1:t,1:m))
        if(S0(t,j)~=0)  
          A(i,S0(t,j))=5000;
          A(S0(t,j),i)=5000;
        end;  
      end;
  end;
end;
R=R-K;

end;
%%
%%%%%%%%%%%%%draw 
Gm=A;
Gm(find(Gm==inf))=0;
Gm(find(Gm==5000))=0;
GG=graph(Gm)
Str=strcat('G(',num2str(K),',',num2str(ap),')');
save(Str,'GG'); 
%plot(GG);
%%
%%%%%%%%%%%%find the additional nodes %%%%%%%%%%%%%%%%%%%%%%%%%
RR=nn-n;
if (RR~=0)
[m1,m2]=size(S0);
 ren=zeros(m2,1);
 S0=[S0' ren];
 S0=S0';
 m1=m1+1;
 j=1;
 for i=1:nn
  if i~=S0(:,:)
     S0(m1,j)=i;
     j=j+1;
  end;
 end;
end;

%%  finish clustering 
%% computer the Gloss for each cluster
[nr,lens]=size(S0);
i=0;j=0;
for i=1:nr
    Max=0;Min=100;
    for j=1:lens
      if(S0(i,j)~=0)
        if Max<X.N(S0(i,j))
           Max=X.N(S0(i,j));
        end;
        if Min>X.N(S0(i,j))
           Min=X.N(S0(i,j));
        end;
      end;    
    end;   
   genloss(i)=double(Max-Min)/double(max(X.N)-min(X.N));     
end;

%% compute C0
for i=1:nr  
  ss=X.X0(S0(i,1));
  C0(i)=0;
  for j=1:lens
      if(S0(i,j)~=0)
        if strcmp(X.X0(S0(i,j)),ss)~=1
          C0(i)=1;   
        end;
      end;
  end;
end;
%% compute C1
for i=1:nr 
    mi=floor(double(X.X1(S0(i,1)))/1000);
    C1(i)=0;
    for j=1:lens
      if(S0(i,j)~=0)
        if((X.X1(S0(i,1))~=X.X1(S0(i,j))))
          mk=floor(double(X.X1(S0(i,j)))/1000);
          if mi~=mk
            C1(i)=1; 
          else 
            C1(i)=0.5;  
          end;
       end;
     end;
  end;
end;
%% compute C2
for i=1:nr  
  ss2=X.X2(S0(i,1));
  C2(i)=0;
  for j=1:lens
      if(S0(i,j)~=0)
        if strcmp(X.X2(S0(i,j)),ss2)~=1
          C2(i)=1;   
        end;
      end;
  end;
end; 
%% compute C3
for i=1:nr  
  ss3=X.X3(S0(i,1));
  C3(i)=0;
  for j=1:lens
      if(S0(i,j)~=0)
        if strcmp(X.X3(S0(i,j)),ss3)~=1
          C3(i)=1;   
        end;
      end;
  end;
end; 
%% Compute GilSA, the loss caused by generated information
for i=1:nr-1
 % Gil(i)=lens(i)*(genloss(i)+C0(i)+C1(i)+C2(i)+C3(i));%%%%%%%????????????????????????????2020.8.20
  Gil(i)=lens*genloss(i)+(C0(i)+C1(i)+C2(i)+C3(i));
end;
if(RR==0)    
  Gil(nr)=lens*genloss(nr)+(C0(nr)+C1(nr)+C2(nr)+C3(nr));
else
  Gil(nr)=RR*genloss(nr)+(C0(nr)+C1(nr)+C2(nr)+C3(nr));
end;
GilSA=0;
for i=1:nr                         
GilSA=GilSA+Gil(i);              
end;                              
GilS=GilSA/(nn*5);%%%%%GILSA
%% compute intasil
for mn=1:nr
  SEN(mn)=0;
  for i=1:lens-1
    for j=i+1:lens
      if (S0(mn,i)~=0)&&(S0(mn,j)~=0)
        if G0(S0(mn,i),S0(mn,j))==1  
         SEN(mn)=SEN(mn)+1;
        end;
      end;
    end;
  end;
end;
for i=1:nr-1
    if lens==2
        C2len=1;
    else   
    C2len=factorial(lens)/(factorial(lens-2)*factorial(2));
    end;   
  intrasil(i)=2*SEN(i)*(1-SEN(i)/C2len);
end;
if (RR==0)
  intrasil(nr)=2*SEN(nr)*(1-SEN(nr)/C2len); 
else
    if (RR>0)&&(RR<=2) 
        C2len=1;
        intrasil(nr)=2*SEN(nr)*(1-SEN(nr)/C2len); 
    else
        C2len=factorial(RR)/(factorial(RR-2)*factorial(2));
        intrasil(nr)=2*SEN(nr)*(1-SEN(nr)/C2len); 
   end;
end;
  
%% compute the intersil
for i=1:nr-1
  for j=i+1:nr-1
      E(i,j)=0;
    for si=1:lens
      for sj=1:lens
        if (S0(i,si)~=0)&&(S0(j,sj)~=0)
          if G0(S0(i,si),S0(j,sj))==1
            E(i,j)=E(i,j)+1; 
          end;
        end;
      end;  
    end;
    intersil(i,j)=2*E(i,j)*(1-E(i,j)/(lens*lens));
  end;  
 
      E(i,nr)=0;
    for si=1:lens
      for sj=1:lens
        if (S0(i,si)~=0)&&(S0(nr,sj)~=0)
          if G0(S0(i,si),S0(nr,sj))==1
            E(i,nr)=E(i,nr)+1; 
          end;
        end;
      end;  
    end;
  if(RR~=0)
    intersil(i,nr)=2*E(i,nr)*(1-E(i,nr)/(lens*RR));
  else
    intersil(i,nr)=2*E(i,nr)*(1-E(i,nr)/(lens*lens));
  end; 
end;

interS=0;
for i=1:nr
   for j=i+1:nr
    interS=interS+intersil(i,j);   
   end;    
end;    
sumintra=0;
for i=1:nr
   sumintra=sumintra+intrasil(i);
end;  
NSilS=(sumintra+interS)/(nn*(nn-1)/4);%%NSIL(G,S)


   


