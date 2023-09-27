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
%i1=round(2*length(t)/3); %Get third cycle of three cycles
i1=1;
figure; plot(t(i1:end),H(i1:end)/(max(H(:))),'k:','LineWidth',1); 
hold on; plot(t(i1:end),M0(i1:end)/max(M0(i1:end)),'k-','LineWidth',1);
figure; yyaxis right; plot(t(i1:end),H(i1:end)/max(abs(H(i1:end))),'k:','LineWidth',1); 
%figure; yyaxis right; plot(t(i1:end),H(i1:end),'k:','LineWidth',1);
set(gca,'YColor','k');
yyaxis left; hold on; plot(t(i1:end),M0(i1:end)/max(abs(M0(i1:end))),'k-','LineWidth',1);
%yyaxis left; hold on; plot(t(i1:end),M0(i1:end),'k-','LineWidth',1);
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);
set(gca,'YColor','k');
%xlim([1.999*10^-3 3.001*10^-3]);

Rc=10*10^(-9); %[meters] %Micromod specs say BNF starch has Rc~10
Rh=48.5*10^(-9); %[meters]
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
yyaxis left; hold on; plot(t(i1:end),Mtot(i1:end)/max(abs(Mtot(i1:end))),'k--','LineWidth',1);
%yyaxis left; hold on; plot(t(i1:end),Mtot(i1:end),'k--','LineWidth',1);

% %Make plots
[r53ref,freqs]=calculate_spectra;
figure; plot(freqs,r53ref,'o-');

[r53NP25,freqs]=calculate_spectra(polymer_distrib(3,:));
hold on; plot(freqs,r53NP25,'x-');

[r53NP50,freqs]=calculate_spectra(polymer_distrib(5,:));
hold on; plot(freqs,r53NP50,'x-');

% figure; plot(freqs,r53ref,'ko-');
% scaling=calc_scaling(freqs,r53ref,r53NP25)
% hold on; plot(scaling*freqs,r53NP25,'k+');
% scaling=calc_scaling(freqs,r53ref,r53NP50)
% hold on; plot(scaling*freqs,r53NP50,'kx');
% xlim([400 1600]);

[r53ref,reffreqs]=calculate_spectra;
figure; plot(reffreqs,r53ref,'k-','Marker','+','LineWidth',1);
ylim([0 1.1]);
[r53NP50dimonly,freqs]=calculate_spectra([0.5 .25],linspace(400,1000,5));
hold on; plot(freqs,r53NP50dimonly,'ko--','LineWidth',1);
%[r53NP100dimonly,freqs100]=calculate_spectra([0 .5],linspace(400,700,5));
%hold on; plot(freqs100,r53NP100dimonly,'k--','Marker','square');
[r53NP50,freqs]=calculate_spectra(polymer_distrib(end,:),linspace(400,1000,5));
hold on; plot(freqs,r53NP50,'k--','Marker','diamond','LineWidth',1);
[r53visc,freqs]=calculate_spectra(1,linspace(400,1000,5),10*10^-3,300,1.5*10^(-3));
hold on; plot(freqs,r53visc,'k--','Marker','square','LineWidth',1);
set(gca,'FontWeight','Bold');
set(gca,'FontSize',12);

clear scaling;
figure; plot(reffreqs,r53ref,'k-','Marker','+','LineWidth',1);
scaling=calc_scaling_v2(reffreqs,r53ref,freqs',r53NP50dimonly)
hold on; plot(mean(scaling)*freqs',r53NP50dimonly','ko--','LineWidth',1);
%scaling=calc_scaling_v2(reffreqs,r53ref,freqs,r53NP100dimonly)
%hold on; plot(mean(scaling)*freqs',r53NP100dimonly,'k--','Marker','square','LineWidth',1);
scaling=calc_scaling_v2(reffreqs,r53ref,freqs,r53NP50)
hold on; plot(mean(scaling)*freqs',r53NP50,'k--','Marker','diamond','LineWidth',1);
scaling=calc_scaling_v2(reffreqs,r53ref,freqs,r53visc)
hold on; plot(mean(scaling)*freqs',r53visc,'k--','Marker','square','LineWidth',1);
xlim([400 1600]);
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


%Now do viscosity
scalings=zeros(5,5);
scalings(1,:)=1;
for j=2:5
    [r53NP,freqs]=calculate_spectra(1,linspace(400,1000,5),10*10^(-3),300,(1+0.125*(j-1))*10^(-3));
    scalings(j,:)=calc_scaling_v2(reffreqs,r53ref,freqs,r53NP);
end
figure; yyaxis left; 
hold on; plot(linspace(0,.5,5),mean(scalings'),'ko-','LineWidth',1); xlabel('% NPs bound');
set(gca,'YColor','k'); ylabel('Mean');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);
yyaxis right; plot(linspace(0,.5,5),std(scalings'),'kx-','LineWidth',1);
set(gca,'YColor','k'); ylabel('Std. dev.');
set(gca,'FontWeight','Bold'); set(gca,'FontSize',12);

ax_top = axes(); % axis to appear at top
hold(ax_top);
ax_top.XAxisLocation = 'top';
ax_top.YAxisLocation = "right";
ax_top.YTick = [];
% ax_top.XDir = 'reverse';
ax_top.Color = 'none';
xlabel(ax_top, 'Energy($e$V)', 'Interpreter', 'latex', 'FontSize', 14);


%% Make a figure with expt data
%load('21-May-2021_10h49m_Cavity_7frq_24MA_10mTAC_1mTDC_neg99');
addpath('C:\Users\Dell\Desktop\NP aggregation\ConA_Gly_VECTs_2FLDS_10_28_19\10mT_AC_1mT_DC\10mT_Cavity')
load('29-Oct-2019_15h34m_Cavity_31C_7frq_10mTAC_1mTDC')
cavH2=mean(H2,2);
cavH4=mean(H4,2);
cavH2=0; cavH4=0;

%addpath('C:\Users\Dell\Desktop\NP aggregation\VIsc_Agg_5_21_21');
addpath 'C:\Users\Dell\Desktop\NP aggregation\John Weaver - Visc_vs_Binding_ConA_Data_Code\Visc_vs_Binding_ConA_Data_Code\ConA_Gly_VECTs_2FLDS_10_28_19\10mT_AC_1mT_DC\10mT_ConA'
%load('VIsc_Agg_5_21_21\21-May-2021_10h59m_A0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_00.mat')
load('29-Oct-2019_13h18m_ConA_GLY_Pln_NPs_Control_0_7frq_10mTAC_1mTDC_31C_0.mat');
%r42ref=mean(r');
%r42ref=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
r42ref=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
load('29-Oct-2019_13h46m_ConA_Pln_NPs_Expt_1_7frq_10mTAC_1mTDC_31C_1.mat')
%r42agg(1,:)=mean(r');
%r42agg(1,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
r42agg(1,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
load('29-Oct-2019_14h14m_ConA_Pln_NPs_Expt_2_7frq_10mTAC_1mTDC_31C_2')
%r42agg(2,:)=mean(r');
%r42agg(2,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
r42agg(2,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
load('29-Oct-2019_14h48m_ConA_Pln_NPs_Expt_3_7frq_10mTAC_1mTDC_31C_3');
%r42agg(3,:)=mean(r');
%r42agg(3,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
r42agg(3,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
load('29-Oct-2019_15h16m_ConA_Pln_NPs_Expt_4_7frq_10mTAC_1mTDC_31C_4');
%r42agg(4,:)=mean(r');
%r42agg(4,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2));
r42agg(4,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
abs(mean(H2,2))
clear scalings;
figure;
for j=1:4;
    hold on; plot(Freq,mean(r42agg(j,:)),'o-');
    scalings(j,:)=calc_scaling_v2(Freq,r42ref,Freq,r42agg(j,:));
end
figure; yyaxis right; plot([0 1.7 3.4 6.9 NaN]*10^2,[0 std(scalings')]);
yyaxis left; plot([0 1.7 3.4 6.9 NaN]*10^2,[1 mean(scalings')]);


%% Take a look at the viscosity data too
%load 29-Oct-2019_11h33m_Cavity_31C_7frq_10mTAC_1mTDC
addpath('C:\Users\Dell\Desktop\NP aggregation\ConA_Gly_VECTs_2FLDS_10_28_19\10mT_AC_1mT_DC\10mT_Cavity');
load 29-Oct-2019_15h34m_Cavity_31C_7frq_10mTAC_1mTDC
%cavH2=mean(H2,2);
%cavH4=mean(H4,2);
cavH2=0;
cavH4=0;

%load('21-May-2021_11h50m_B0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_05.mat')
addpath('C:\Users\Dell\Desktop\NP aggregation\John Weaver - Visc_vs_Binding_ConA_Data_Code\Visc_vs_Binding_ConA_Data_Code\ConA_Gly_VECTs_2FLDS_10_28_19\10mT_AC_1mT_DC\10mT_Gly')
load 29-Oct-2019_11h53m_GLY_Pln_NPs_Control_0_7frq_10mTAC_1mTDC_31C_0.mat
%r42ref=mean(r');
r42ref=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
%r42visc(1,:)=mean(r');
%r42visc(1,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2));
load 29-Oct-2019_12h7m_GLY_Pln_NPs_Expt_1_7frq_10mTAC_1mTDC_31C_1.mat
%r42visc(1,:)=mean(r');
r42visc(1,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
%load('21-May-2021_12h41m_C0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_10.mat')
%cavH4=1.0*cavH4; cavH2=1.3*cavH2;
%r42visc(2,:)=mean(r');
%r42visc(2,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
load 29-Oct-2019_12h21m_GLY_Pln_NPs_Expt_2_7frq_10mTAC_1mTDC_31C_2.mat
r42visc(2,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
%r42visc(2,:)=mean(r');
%load('21-May-2021_13h31m_D0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_15.mat');
%cavH4=1.0*cavH4; cavH2=1.3*cavH2;
%r42visc(3,:)=mean(r');
%r42visc(3,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2)); 
load 29-Oct-2019_12h50m_GLY_Pln_NPs_Expt_3_7frq_10mTAC_1mTDC_31C_3
%r42visc(3,:)=mean(r');
r42visc(3,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
%load('21-May-2021_14h22m_E0_uMod_S_Nptz_AuFe3O4_B_7frq_24MA_10mTAC_1mTDC_20.mat');
%cavH4=1.0*cavH4; cavH2=1.3*cavH2;
%r42visc(4,:)=mean(r');
%r42visc(4,:)=abs(mean(H4-cavH4,2)./mean(H2-cavH2,2))
%r42visc(4,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
load 29-Oct-2019_13h4m_GLY_Pln_NPs_Expt_4_7frq_10mTAC_1mTDC_31C_4
%r42visc(4,:)=mean(r');
r42visc(4,:)=abs(mean((H4-cavH4)./(H2-cavH2),2)); 
%abs(mean(H2,2))

figure;
clear scalings;
for j=1:4;
%hold on; plot(Freq,(r42visc(j,:)),'o-');
scalings(j,:)=calc_scaling_v2(Freq,r42ref,Freq,r42visc(j,:));
end

%xx=[1 1.075 1.159 1.353 1.593]
xx=[1 1.075 1.159 1.353 NaN]
stdd=[0 std(scalings')]; %stdd(5)=NaN;
meann=[1 mean(scalings')]; %meann(5)=NaN;
hold on; yyaxis right; plot(xx,stdd);
yyaxis left; plot(xx,meann);

%figure; plot(xx,meann);
