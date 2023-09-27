%Load mat file solution to aggregation model %Can run again, just
%"agg_model",takes about 5 minutes.

%Calculate NP aggregation polymer distribution
% load('agg_model_results.mat');
% polymer_distrib=zeros(10,1); %Calculate first 10 fractional numbers of monomers, dimers and so on
% k=round(0.3*size(rho,1)); %First dimension is N. Index of N=0.3.
% dVx=Vx(k,2)-Vx(k,1);
% for j=1:length(polymer_distrib)
%     polymer_distrib(j)=sum(rho(k,Vx(k,:)>(2.5+(j-1)*5) & Vx(k,:)<(7.5+(j-1)*5))*dVx);
% end
%
% %Make plots
% [r53,freqs]=calculate_spectra;
% figure(1); hold on; plot(freqs,r53,'o-');
%
% [r53,freqs]=calculate_spectra(polymer_distrib);
% figure(1); hold on; plot(freqs,r53,'x-');
%
% legend('Reference (No biomarker present; 0% NP aggregated)','Biomarker conc. 10% of NP conc. (11% NP bound)', 'Biomarker conc. 20% of NP conc. (13% NP bound)');
% title('Biomarker (e.g. protein) detection via targeted mNP aggregation (increase in avg. relaxation time)');
% %"N.B. (In the caption: Between dashed and dotted line: we doubled the number of biomarker molecules.)
% %Calculate change in average relaxation time

load('agg_model_results.mat');
polymer_distrib=zeros(5,10); %Calculate first 10 fractional numbers of monomers, dimers and so on
%k=round(0.3*size(rho,1)); %First dimension is N. Index of N=0.3.
dVx=Vx(1,2)-Vx(1,1);
numgaps=size(polymer_distrib,1)-1;
for k=1:size(polymer_distrib,1)
    for j=1:size(polymer_distrib,2)
        ind=1+round((k-1)*(size(rho,1)-numgaps)/numgaps);
        polymer_distrib(k,j)=sum(rho(ind,Vx(1,:)>(2.5+(j-1)*5) & Vx(1,:)<(7.5+(j-1)*5))*dVx);
    end
end
%N=1-polymer_distrib(:,1); %True N is just fraction of bound NPs
%which is equal to number of linkers
Nmax=0.5;
NN=linspace(0,Nmax,numgaps+1);
figure; plot(NN,polymer_distrib(:,2),'o-');
hold on; plot(NN,polymer_distrib(:,3),'x-');
hold on; plot(NN,polymer_distrib(:,4),'+-');
hold on; plot(NN,polymer_distrib(:,5),'Marker','square','LineStyle','-');
matrix_linkers=kron(0:size(polymer_distrib,2)-1,ones(size(polymer_distrib,1),1));
numlinkers=diag(matrix_linkers*polymer_distrib');
hold on; plot(NN,numlinkers,'Marker','diamond','LineStyle','-');

matrix_1st=kron(1:size(polymer_distrib,2),ones(size(polymer_distrib,1),1));
firstmoment=diag(matrix_1st*polymer_distrib')./sum(polymer_distrib,2);
secondmoment=diag((matrix_1st.^2)*polymer_distrib')./sum(polymer_distrib,2);
figure; yyaxis left; plot(NN,firstmoment,'ko-','LineWidth',1); ylabel('1st moment (Mean)');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);
set(gca,'YColor','k');
hold on; yyaxis right; plot(NN,secondmoment,'kx-','LineWidth',1); ylabel('2nd moment');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);
set(gca,'YColor','k');
xlabel('N');
legend('1st moment','2nd moment');

[M0,t,H]=effective_field_model;
i1=round(2*length(t)/3); %Get third cycle of three cycles
figure; plot(t(i1:end),H(i1:end)/(max(H(:))),'k:','LineWidth',1); 
hold on; plot(t(i1:end),M0(i1:end)/max(M0(i1:end)),'k-','LineWidth',1);
%xlim([1.999*10^-3 3.001*10^-3]);

Rc=20*10^(-9); %[meters] %Micromod specs say BNF starch has Rc~10
Rh=50*10^(-9); %[meters]
Mtot=0;
load('avmags');
distrib=polymer_distrib(end,:); %Use 50% NP bound
%Make sure normalization is correct, distrib needs to add up to one to
%calculate total signal as weighted average of monomers, dimers, trimers etc
S=0; for j=1:length(distrib); S=S+j*distrib(j); end
distrib=distrib/S;
for k=1:length(distrib)
    %N.B. Dimer volume is twice that of monomer, and so on...
    [M,t]=effective_field_model(Rc*(avmags(k)^(1/3)),Rh*(k^(1/3))); %Randomly-aligned cores
    Mtot=Mtot+distrib(k)*M;
end
%hold on; plot(t(i1:end),Mtot(i1:end)/max(Mtot(i1:end)),'k--','LineWidth',1);
hold on; plot(t(i1:end),Mtot(i1:end),'k--','LineWidth',1); %EDIT: 5/4/21 

% %Make plots
[r53ref,freqs]=calculate_spectra;
figure; plot(freqs,r53ref,'o-');

[r53NP25,freqs]=calculate_spectra(polymer_distrib(3,:));
hold on; plot(freqs,r53NP25,'x-');

[r53NP50,freqs]=calculate_spectra(polymer_distrib(5,:));
hold on; plot(freqs,r53NP50,'x-');

figure; plot(freqs,r53ref,'ko-');
scaling=calc_scaling(freqs,r53ref,r53NP25)
hold on; plot(scaling*freqs,r53NP25,'k+');
scaling=calc_scaling(freqs,r53ref,r53NP50)
hold on; plot(scaling*freqs,r53NP50,'kx');
xlim([400 1000]);

[r53ref,reffreqs]=calculate_spectra;
figure; plot(reffreqs,r53ref,'ko-');
ylim([0 1.1]);
[r53NP50dimonly,freqs]=calculate_spectra([0.5 .25],linspace(400,1000,5));
hold on; plot(freqs,r53NP50dimonly,'kx--');
[r53NP100dimonly,freqs100]=calculate_spectra([0 .5],linspace(400,700,5));
hold on; plot(freqs100,r53NP100dimonly,'k--','Marker','square');
[r53NP50,freqs]=calculate_spectra(polymer_distrib(end,:),linspace(400,1000,5));
hold on; plot(freqs,r53NP50,'k--','Marker','diamond');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',12);

scalings=zeros(5,5);
scalings(1,:)=1;
scalingsdimonly=zeros(5,5);
scalingsdimonly(1,:)=1;
for j=2:5
    [r53NPdimonly,freqs]=calculate_spectra([1-(j-1)/8 (j-1)/16],linspace(400,1000,5));
    scalingsdimonly(j,:)=calc_scaling_v2(reffreqs,r53ref,freqs,r53NPdimonly);
    [r53NP,freqs]=calculate_spectra(polymer_distrib(j,:),linspace(400,1000,5));
    scalings(j,:)=calc_scaling_v2(reffreqs,r53ref,freqs,r53NP);
end
figure; yyaxis left; 
plot(linspace(0,.5,5),mean(scalings'),'ko-','LineWidth',1); xlabel('% NPs bound');
hold on; plot(linspace(0,.5,5),mean(scalingsdimonly'),'ko--','LineWidth',1);
set(gca,'YColor','k'); ylabel('Mean');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);
yyaxis right; plot(linspace(0,.5,5),std(scalings'),'kx-','LineWidth',1);
hold on; plot(linspace(0,.5,5),std(scalingsdimonly'),'kx--','LineWidth',1);
set(gca,'YColor','k'); ylabel('Std. dev.');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);


%% Make a figure with expt data
addpath('C:\Users\Dell\Desktop\NP aggregation\VIsc_Agg_5_21_21');
load('VIsc_Agg_5_21_21\21-May-2021_10h59m_A0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_00.mat')
r42ref=mean(r');
load('21-May-2021_11h9m_A1_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_01.mat')
r42agg(1,:)=mean(r');
load('21-May-2021_11h19m_A2_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_02.mat')
r42agg(2,:)=mean(r');
load('21-May-2021_11h30m_A3_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_03.mat');
r42agg(3,:)=mean(r');
load('21-May-2021_11h40m_A4_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_04.mat');
r42agg(4,:)=mean(r');
clear scalings;
figure;
for j=1:4;
    hold on; plot(Freq,mean(r42agg(j,:)),'o-');
    scalings(j,:)=calc_scaling_v2(Freq,r42ref,Freq,r42agg(j,:));
end
figure; yyaxis right; plot([0 2 4 7 10],[0 std(scalings')]);
yyaxis left; plot([0 2 4 7 10],[1 mean(scalings')]);


%% Take a look at the viscosity data too

load('21-May-2021_11h50m_B0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_05.mat')
r42visc(1,:)=mean(r');
load('21-May-2021_12h41m_C0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_10.mat')
r42visc(2,:)=mean(r');
load('21-May-2021_13h31m_D0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_15.mat');
r42visc(3,:)=mean(r');
load('21-May-2021_14h22m_E0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_20.mat');
r42visc(4,:)=mean(r');

figure;
clear scalings;
for j=1:4;
hold on; plot(Freq,(r42visc(j,:)),'o-');
scalings(j,:)=calc_scaling_v2(Freq,r42ref,Freq,r42visc(j,:));
end

xx=[1 1.13 1.24 1.58 NaN];
stdd=[0 std(scalings')]; stdd(5)=NaN;
meann=[1 mean(scalings')]; meann(5)=NaN;
hold on; yyaxis right; plot(xx,stdd);
yyaxis left; plot(xx,meann);
