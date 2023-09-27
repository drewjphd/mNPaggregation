%Calculates average magnetization for an aggregated clump of NPs with n
%cores. DJ April 2020
function avmags=av_mag(Nmax)
tic; clear all;
if nargin<1
    Nmax=10;
end
%n=2; %n has to be 2 or greater
samples=10^5; %number of vectors to average over
dims=3; %number of dimensions
avgmag=zeros(Nmax-1,1); %initialize

for n=2:Nmax
    R=2*rand(n,dims,samples)-1; %uniformly distributed random vectors, each normalized to unity
    %R=zeros(n,dims,samples);
    
    %normalize:
    for j=1:samples
        for k=1:n
            r=2*rand(dims,1)-1;
            theta=pi*rand(1);
            phi=2*pi*rand(1);
            R(k,:,j)=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
            %R(k,:,j)/norm(R(k,:,j));
        end
    end
    %R(1,:,:)
    
    dipolesum=sum(R,1);
    %dipolesum
    norms=zeros(samples,1);
    for j=1:samples
        norms(j)=norm(dipolesum(1,:,j));
    end
    %norms
    
    avgmag(n-1)=sum(norms)/length(norms); %Average magnetization of n aggregated cores
end
toc;
figure; plot(1:Nmax,[1 avgmag'],'o-');
avmags=[1 avgmag'];
save('avmags.mat','avmags');
end