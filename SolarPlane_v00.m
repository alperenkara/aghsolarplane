%%AGH Solar Plane
%%PV Cell Shockley diode equation
clear all
clc
Va=0:.01:25;
Suns=.2:.2:1;
TaC=30; % 
Ipv=zeros(size(Va));
Ppv=zeros(size(Va));
for s=1:1:length(Suns)
for i=1:1:length(Va)
k=1.38e-23; % Boltzman’s const
q=1.6e-19; % charge on an electron
A=1.2;% "diode quality" factor, =2 for crystaline, <2 for amorphous
Vg=1.12;% band gap voltage, 1.12eV for xtal Si, ˜1.75 for amorphous Si.
Ns=44;% AGH Solar Plane number of series connected cells (diodes) 

T1=273+25;
Voc_T1=21.06/Ns; % open circuit voltage per cell at temperature T1
Isc_T1=3.80; %short circuit current per cell at temp T1

T2=273+75;
Voc_T2=17.05/Ns; % open circuit voltage per cell at temperature T2
Isc_T2=3.92; % short circuit current per cell at temp T2

TarK=273+TaC;
Tref=273+25; %ref. temperature
%{
Va=0;
Iph_T1=Isc_T1;
%}
Iph_T1=Isc_T1*Suns(s);
a=(Isc_T2-Isc_T1)/Isc_T1*1/(T2-T1);
Iph=Iph_T1*(1+a*(TarK-T1));

Vt_T1=k*T1/q;
Ir_T1=Isc_T1/(exp(Voc_T1/(A*Vt_T1))-1);
Ir_T2=Isc_T2/(exp(Voc_T2/(A*Vt_T1))-1);
b=Vg*q/(A*k);
Ir=Ir_T1*(TarK/T1).^(3/A).*exp(-b.*(1./TarK-1/T1));
X2v=Ir_T1/(A*Vt_T1)*exp(Voc_T1/(A*Vt_T1));
dVdI_Voc=-1.15/Ns/2;
%Sun Power PV Cell Character
Rs=-dVdI_Voc-1/X2v;
%Ia=1:0.01:Iph;
Vt_Ta=A*k*TarK/q;
%{
Ia1=Iph-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1));
solve for I: f(Ia)=Iph-Ia-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1))=0;
Newton raphson Ia2=Ias1-f(Ia1)/f'(Ia1)
%}
Vc=Va(i)/Ns;
Ia=zeros(size(Vc));
%Iav=Ia
for j=1:1:10
    Ia=Ia-(Iph-Ia-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1))./(-1-Ir*(exp((Vc+Ia*Rs)/Vt_Ta)-1).*Rs/Vt_Ta);
end
Ipv(s,i)=Ia;
Ppv(s,i)=Va(i)*Ipv(s,i);
end
end
%% figure 
axes1 = axes('Parent',figure,'OuterPosition',[0 0.5 1 0.5]);
 xlim(axes1,[0 23]);
 ylim(axes1,[0 5]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
title('I-V charateristics at 25 C');
xlabel('V_p_v (V)');
ylabel('I_p_v (A)');
plot1 = plot(Va(1,:),Ipv(:,:),'Parent',axes1,'LineWidth',1.5);
set(plot1(1),'DisplayName','0.2 Sun');
set(plot1(2),'DisplayName','0.4 Sun');
set(plot1(3),'DisplayName','0.6 Sun');
set(plot1(4),'DisplayName','0.8 Sun');
set(plot1(5),'DisplayName','1.0 Sun');


 legend1 = legend(axes1,'show');
 set(legend1,...
     'Position',[0.791450219003363 0.769901214241614 0.0793528505392912 0.151937984496124]);
 legend2 = legend(axes2,'show');
 set(legend2,...
     'Position',[0.189976911157747 0.248883515126569 0.0793528505392912 0.151937984496124]);
