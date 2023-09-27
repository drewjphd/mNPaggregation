function aggregation_model(Nmax)
tic;
%Default value of Nmax=0.01
if nargin<1
    %Nmax=0.0388; %This value corresponds to true N=0.5 
    Nmax=.5; %This value corresponds to true N=0.5 
    %Nmax=0.1;
    %Nmax=1;
end

%c=1; %Solve diff eq for c=1, then calculate c(N)... (see paper/notebook)
V=linspace(0,3*50,3*2000)'; %Volume is positive definite
N=linspace(0,Nmax,1000); %Number/concentration of linkers positive def
dN=N(2)-N(1); %Need dN "small enough." Check if solution is smooth/makes sense.
dV=V(2)-V(1); %Also needs to be small enough to estimate integrals well

%Set initial condition
%rho0=exp(-((log(V)-1.5).^2)/(2*0.0625))./(V*sqrt(2*pi*0.0625)); rho0(1)=0; %Lognormal dist
%rho0=exp(-((V-5).^2)/(2))/sqrt(2*pi); %Gaussian initial condition with mean 5 and stddev1
%rho0=exp(-((V-0.2).^2)/(2*0.001))/sqrt(2*pi*0.001); %Gaussian initial condition with mean 5 and stddev1
rho0=exp(-((V-5).^2)/(2*0.1))/sqrt(2*pi*0.1); %Gaussian initial condition with mean 5 and stddev rt0.1

[Vx,Ny]=meshgrid(V,N);
rho=zeros(length(N),length(V));
rho(1,:)=rho0;
S=abs(V.^(2/3)); %abs used to make sure V^(2/3) is positive and real
%Also check if \rho vanishes for V->0 and V->infty for every N. Otherwise
%there is some truncation and numerics unreliable

%Integrate diff eq
for j=2:length(N)
    firstterm=(sqrt(12*pi)*S.*rho(j-1,:)');
    secondterm=J(rho(j-1,:),S',dV);
    secondterm=secondterm.*sum(V.*firstterm*dV)./sum(V.*secondterm*dV); %Adjust coefficent to conserve probability
    %c=2*0.0782*N(j)+0.07707; %found by setting c=1, solving dN'=c(N)dN (see paper) Correct up to at least N=0.5
    %c=3*0.1004*N(j)^2+2*(-0.02408)*N(j)+0.06461;
    %c=0.03211+2*0.1054*N(j)-3*0.201*N(j)^2+4*0.1901*N(j)^3;
    %c=3*0.0742*N(j)^2+2*0.02436*N(j)+0.08597;
    c=2*0.0533*N(j)+0.0507;
    drho=c*dN*(-firstterm+secondterm);
    rho(j,:)=rho(j-1,:)+drho';
end;

%Make plot
figure(1); plot(Vx(1,:),rho(1,:),'-','LineWidth',1,'Color',[0.7 0.7 0.7]);
%hold on; plot(Vx(1,:),rho(end/2,:),'k--','LineWidth',1);
hold on; plot(Vx(1,:),rho(round(2*size(rho,1)/6),:),'k--','LineWidth',1);
legend('initial, N=0',['N=',num2str(round(2*Nmax/6,3))]);
xlabel('Hydrodynamic volume, V');
ylabel('Size distribution, \rho(V,N)');
%title(['Change in size distribution when some linkers are added. (c_1=',num2str(c),')']);
set(gca,'FontWeight','bold')
toc;

polymer_distrib=zeros(9,10); %Calculate first 10 fractional numbers of monomers, dimers and so on
%k=round(0.3*size(rho,1)); %First dimension is N. Index of N=0.3.
dVx=Vx(1,2)-Vx(1,1);
for k=1:size(polymer_distrib,1)
    for j=1:size(polymer_distrib,2)
        ind=1+round((k-1)*(size(rho,1)-8)/8);
        polymer_distrib(k,j)=sum(rho(ind,Vx(1,:)>(2.5+(j-1)*5) & Vx(1,:)<(7.5+(j-1)*5))*dVx);
    end
end
%N=1-polymer_distrib(:,1); %True N is just fraction of bound NPs
%which is equal to number of linkers
NN=linspace(0,Nmax,9);
figure; plot(NN,1-polymer_distrib(:,1));
save('agg_model_results.mat','N','polymer_distrib','rho','Vx','Ny');
end

function out=J(rho,S,dV)
out=zeros(length(S),1);
%V=V';
%S=abs(V.^(2/3));
for k=1:length(S)
    out(k)=sum(rho(k:-1:1).*rho(1:k).*S(k:-1:1).*S(1:k));
end
out=out*12*pi*dV;
end