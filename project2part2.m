cbeam_frf=importdata('cbeam_3x2_frfdata_Ts20kHz_hlfadv_20kPts__struc.mat');
cbeam_controller=importdata('cbeam_controller_2x3_Ts20kHz_40pole__struc.mat');

%extracting the controller mtarix's

A_bar=cbeam_controller.a;
B_bar=cbeam_controller.b;
C_bar=cbeam_controller.c;
D_bar=cbeam_controller.d;
ts_bar=cbeam_controller.ts;
controller=ss(A_bar,B_bar,C_bar,D_bar,ts_bar);

%open loop frf

w=cbeam_frf.Frequency;
G=cbeam_frf.ResponseData;
H=freqresp(controller,w);

comp_ol_frf=zeros(size(H,1),size(G,2),numel(w));
for k=1:numel(w)
    comp_ol_frf(:,:,k)=H(:,:,k)*G(:,:,k);
end

%plot of the characterstic loci
char_loci=zeros(size (H,1),numel(w));
for k=1:numel(w)
    char_loci(:,k)=eig(comp_ol_frf(:,:,k));
end
char_loci_mag=20*log10(abs(char_loci));
char_loci_phase = radtodeg(angle(char_loci));
plot(char_loci_phase(1,:),char_loci_mag(1,:),'green*');
hold on
plot(char_loci_phase(2,:),char_loci_mag(2,:),'blue--.');
hold on

ngrid
title('Nichols Chart')
xlabel('Phase(in degrees)')
ylabel('Magnitude(in DB)')
hold off