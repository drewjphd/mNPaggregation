%Calculate magnetic spectra, ratio of 5th to 3rd harmonics
%from parameters and NP aggregation distribution.
%E.g. if agg_distrib=[1] then it's all monomers
%if agg_distrib=[0.8 0.1] then there is 1 dimer for every 8 monomers
%if agg_distrib=[0.3 0.2 0.1] then there is 1 trimer and 2 dimers
%for every 3 monomers. All three cases have equal number of identical
%individual particles. (This code can only handle identical NPs, i.e.
%identical Rc and Rh). See paper for more explanation.
function [r53,freqs]=calculate_spectra_v2(agg_distrib,freqs, H, temp, visc)
if nargin<1; agg_distrib=1; end
if nargin<2; freqs=linspace(400,1000,5); end %Field frequencies [Hertz]
if nargin<3; H=10*10^-3; end %Field amplitude [Tesla]
if nargin<4; temp=300; end %Room temperature [Kelvin]
if nargin<5; visc=1.0*10^(-3); end %[Pascal-second] Water at r.t.
load('avmags');

%Make sure normalization is correct, distrib needs to add up to one to
%calculate total signal as weighted average of monomers, dimers, trimers etc
S=0; for j=1:length(agg_distrib); S=S+j*agg_distrib(j); end
agg_distrib=agg_distrib/S;

h=zeros(7,length(freqs)); %first seven harmonics at the freqs
%Compute spectra
Rc=10*10^(-9); %[meters] %Micromod specs say BNF starch has Rc~10
Rh=48.5*10^(-9); %[meters]

for j=1:length(freqs)
    Mtot=0; t=0;
    for k=1:length(agg_distrib)
        numk=round(100*agg_distrib(k));
        if numk<1; break; end %If less than one for 100 NPs, ignore contribution.
        for q=1:numk
        %N.B. Dimer volume is twice that of monomer, and so on...
        ran=normrnd(1,0.1); %Gaussian random number, mean 1, std dev 0.2
        ran(ran<0)=0.5; %If the random number ends up being negative...
        [M,t]=effective_field_model(Rc*((ran*avmags(k))^(1/3)),Rh*((ran*k)^(1/3)),freqs(j),H,temp,visc); %Randomly-aligned cores
        %[M,t]=effective_field_model(Rc*(k^(1/3)),Rh*(k^(1/3)),freqs(j),H,temp,visc);
        %Perfectly aligned cores
        Mtot=Mtot+M;
        end
    end
    
    %ind=2*round(length(t)/3);
    %h(:,j)=get_harmonics(t(ind:end),Mtot(ind:end));
    
    %Calculate derivative of magnetization (that's the signal from receive coil E.m.f.=-dB/dt)
    tt=(t(1:end-1)+t(2:end))/2;
    dMtotdT=-(Mtot(1:end-1)-Mtot(2:end))./(t(2)-t(1));
    ind=2*round(length(tt)/3); %Just use last one of the three cycles
    h(:,j)=get_harmonics(tt(ind:end),dMtotdT(ind:end));
end
r53=h(5,:)./h(3,:);
end