%Calculates M as fcn of t in the eff. field. approx. to the F-P eqn?? or
%F-J approx
%Inputs: Default values and units
%t Time [ms]
%Output

function [M,t,H]=effective_field_model(Rc,Rh,freq,H0,temp,visc,ncyc)
%Default values and units
if nargin<1; Rc=10*10^(-9); end % [meters] %See Reeves et al PLOS One paper
if nargin<2; Rh=48.5*10^(-9); end % [meters]
if nargin<3; freq=16*10^2; end % [Hz]
%if nargin<3; freq=7.5*10^2; end % [Hz]
if nargin<4; H0=10*10^(-3); end % [Tesla]
if nargin<5; temp=300; end % [Kelvin] room temperature 20C 
if nargin<6; visc=1*10^(-3); end % [Pascal-second] Water at r.t.
if nargin<7; ncyc=3; end % Number of cycles (need a few to settle down)

%% Constants
kB=1.38*10^(-23); %Boltzmann [Joules per Kelvin]
tauB=(3*visc*(4/3)*pi*Rh^3)/(kB*temp); %Brown. relax. time
%Mu=65.4*5.24*10^3*(4/3*pi*Rc.^3)  % \mu in Joules/Tesla
%Mu=(15/25)*3.2*49*10^3*(4/3*pi*Rc.^3);  % \mu in Joules/Tesla
Mu=(11.9*10^(-18))*(Rc^3)/((10*10^(-9))^3); %Magnetic moment per NP [SI units]
%https://doi.org/10.1016/j.jmmm.2014.09.070
%Coefficient above inherited from older version of code by DBR.
%DJ: Calculated factor below from umod spec sheet, for (rh,rc)=(50nm,10nm)
%Mu=60*5.35*(10^3)*(4/3)*pi*Rc^3
%Mu=(60/((10^6)*2.9*10^12)) %*(Rc/(14*10^(-9)))^3
%Mu=(49*15*10^(-6))/(1.5*10^13)
%Mu=76
xi0=Mu*H0/(kB*temp);

tbar=0; %Non-dimensional time. Cleaner to integrate in non-dimensinal units, 
%and then reintroduce units for final results
delta_t=0; %Initialize integration timestep
if tauB<1/freq
    delta_t=tauB/1000; %Should be smaller than smallest timescale
else
    delta_t=1/(1000*freq);
end
delta_t=0.1*3.8848e-07; %Set it

tbar=linspace(0,ncyc,ncyc*round((1/freq)/delta_t)); %ncyc cycles in non-dimensional time
M=zeros(1,length(tbar));
M(1)=.001; %Initialize. Can't make exactly zero, will cause NaN
delta_tbar=tbar(2)-tbar(1);
%linspace(-1,1,)

%Integrate!
for j=2:length(tbar)
    dMdtbar=-2*M(j-1)*(1-xi0*sin(2*pi*tbar(j))/inv_Langevin(M(j-1)))/(freq*tauB);
    M(j)=M(j-1)+dMdtbar*delta_tbar;
end
t=tbar/freq;
H=H0*sin(2*pi*tbar);

M=M*Mu; %DJ EDIT 5/4/21. Re-introduce dimension. Magnetization is proportional to magnetic moment

%tt=(t(1:end-1)+t(2:end))/2;
%MM=-(M(1:end-1)-M(2:end))./delta_t;
%h=get_harmonics(tt,-MM); h(5)/h(3)

%tt=(t(1:end-1)+t(2:end))/2;
%dMdT=-(M(1:end-1)-M(2:end))./(t(2)-t(1));
%ind=2*round(length(tt)/3); %Just use last one of the three cycles
%h=get_harmonics(t(ind:end),M(ind:end)); h(5)/h(3)

%Show result in a figure
%figure; plot(t,sin(2*pi*tbar),'--'); hold on; plot(t,M,'-');  xlim([t(round(2*length(t)/3)) t(end)]);
%figure; plot(tt,dMdT,'-');%  xlim([t(round(2*length(t)/3)) t(end)]);
%xlabel('Time [s]'); ylabel('Amplitude'); legend('Applied Field (20 mT, 1kHz)','Magnetization');
%title('Effective field model (umod; rc=50nm, rh=10nm, 20C, water)');
end

function alpha=inv_Langevin(L)
%L=abs(L);
%alpha=(0.756*L^5-1.383*L^4+1.5733*L^3-2.9463*L^2+3*L)/(1-L); 
alpha=L*(3-1.00651*L^2-0.962251*L^4+1.47353*L^6-0.48953*L^8)/((1-L)*(1+1.01524*L));
%Best available approximant formula (<0.076% max rel err) (R. Jedynak) (Wikipedia)
end