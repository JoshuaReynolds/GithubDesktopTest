%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Calculator for Chip Baluns %
%           Joshua Reynolds            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment is done with a semi-rigid coax soldered to half of a slotline to form a
% dyson balun. The slotline feeds into the chip balun which is then fed
% into a short length of semi-rigid coax.

clear all; close all; clc;

%%% Universal Constants %%%
u0 = 1.2567e-6; % Vacuum Permeability
c = 299792458; % Speed of light in a vacuum
e0 = 8.854e-12; % Vacuum Permittivity
Zo = 50; % Reference characteristic impedance

%%% Material Properties %%%
ur_copper = 0.999994; % Relative permeability of the coax copper
ur_teflon = 1; % Relative permeability of the coax teflon
er_teflon = 2.1; % Relative permittivity of the coax teflon
er_substrate = 3; % Relative permittivity of the slotline substrate
sig_teflon = 1e-24; % Conductivity of the coax teflon
sig_copper = 5.95948e7; % Conductivity of the coax copper
losstangent_teflon = 0.0125; % Loss tangent of Teflon at 3GHz

%%% Physical Dimensions %%%
l_coax1 = 100e-3; % Length of grounded coax in meters (Actual lenght = 91e-3 meters)
l_coax2 = 65e-3; % Length of ungrounded coax in meters (Actual length = 65e-3 meters)
coax1_inner_rad = .91e-3; % Inner radius of grounded coax
coax1_outer_rad = 3.58e-3; % Outer radius of grounded coax
coax2_inner_rad = coax1_inner_rad; % Inner radius of ungrounded coax
coax2_outer_rad = coax1_outer_rad; % Outer radius of ungrounded coax
l_slot1 = 215e-3; % Length of longer slotline in meters (Actual length = 215e-3 meters)
l_slot2 = 9.29e-3; % Length of shorter slotline in meters
w_slot = 2.76e-3; % Width of both slotline in meters (2.76e-3)
h_slot = 760e-6; % Height of both slotline substrate in meters
h_copper = 17.5e-6; % Height of copper on substrate
wire_l = w_slot;
wire_r = coax1_inner_rad;

%%% Data Import and Conversion %%%
MeasuredData = read(rfdata.data,'cx2156insertion.s2p'); % Reading in S-parameter data from PNA file
Measuredf = MeasuredData.freq; % Extracting frequency data from measurement file
MeasuredS = extract(MeasuredData,'S_PARAMETERS',Zo); % Extracting S parameter data from measurement file
MeasuredABCD = s2abcd(MeasuredS,Zo); % Converting S parameters to ABCD parameters
freq = Measuredf'; % Scaling frequency from data file to proper values
w = 2*pi*freq; % Converting frequency to angular frequency

%%% Pre-allocation %%%
skin_depth = zeros(length(freq));


for n = 1:length(freq)
skin_depth(n) = 1/sqrt(pi*freq(n)*ur_copper*u0*sig_copper);
lambda0(n) = c/freq(n);
    
%%% Coax1 Calculations %%%
coax1_per_unit_R(n) = (1/(2*pi))*sqrt((pi*freq(n)*ur_copper*u0)/sig_copper)*(1/coax1_inner_rad + 1/coax1_outer_rad);
coax1_per_unit_L = (ur_teflon*u0/(2*pi))*log(coax1_outer_rad/coax1_inner_rad);
coax1_per_unit_G(n) = 4*pi^2*freq(n)*e0*er_teflon*losstangent_teflon/log(coax1_outer_rad/coax1_inner_rad);
coax1_per_unit_C = (2*pi*e0*er_teflon)/log(coax1_outer_rad/coax1_inner_rad);
coax1_gamma(n) = sqrt((coax1_per_unit_R(n) + j*coax1_per_unit_L*w(n))*(coax1_per_unit_G(n) + j*coax1_per_unit_C*w(n)));
coax1_Zo(n) = sqrt((coax1_per_unit_R(n) + j*coax1_per_unit_L*w(n))/(coax1_per_unit_G(n) + j*coax1_per_unit_C*w(n)));
coax1_ABCD(:,:,n) = [cosh(coax1_gamma(n)*l_coax1), coax1_Zo(n)*sinh(coax1_gamma(n)*l_coax1); (1/coax1_Zo(n))*sinh(coax1_gamma(n)*l_coax1), cosh(coax1_gamma(n)*l_coax1)];

%%% Bare Wire Calculations %%%
wire_L = 2*wire_l*(log((2*wire_l/wire_r)*(1+sqrt(1+(wire_r/(2*wire_l))^2)))-sqrt(1+(wire_r/(2*wire_l))^2)+ u0/4 + (wire_r/(2*wire_l)));
wire_ABCD = [1, j*w(n)*wire_L; 0, 1];

%%% Coax2 Calculations %%%
coax2_per_unit_R(n) = (1/(2*pi))*sqrt((pi*freq(n)*ur_copper*u0)/sig_copper)*(1/coax2_inner_rad + 1/coax2_outer_rad);
coax2_per_unit_L = (ur_teflon*u0/(2*pi))*log(coax2_outer_rad/coax2_inner_rad);
coax2_per_unit_G(n) = 4*pi^2*freq(n)*e0*er_teflon*losstangent_teflon/log(coax2_outer_rad/coax2_inner_rad);
coax2_per_unit_C = (2*pi*e0*er_teflon)/log(coax2_outer_rad/coax2_inner_rad);
coax2_gamma(n) = sqrt((coax2_per_unit_R(n) + j*coax2_per_unit_L*w(n))*(coax2_per_unit_G(n) + j*coax2_per_unit_C*w(n)));
coax2_Zo(n) = sqrt((coax2_per_unit_R(n) + j*coax2_per_unit_L*w(n))/(coax2_per_unit_G(n) + j*coax2_per_unit_C*w(n)));
coax2_ABCD(:,:,n) = [cosh(coax2_gamma(n)*l_coax2), coax2_Zo(n)*sinh(coax2_gamma(n)*l_coax2); (1/coax2_Zo(n))*sinh(coax2_gamma(n)*l_coax2), cosh(coax2_gamma(n)*l_coax2)];

%%% Slotline Calculations %%%
[slot_Zo(n),slot_lambda(n)] = slotcalcs(freq(n),er_substrate,w_slot,h_slot);
p(n) = lambda0(n)/slot_lambda(n); % p value from Dr. Ruyle's dissertation
d(n) = (w_slot + h_copper)/(20*p(n)^2); % d value from Dr. Ruyle's dissertation
slot_beta(n) = 2*pi/slot_lambda(n); % Phase constant from standard formula
slot_per_unit_L(n) = sqrt(slot_Zo(n)*slot_beta(n)/w(n)^2); % Calculated from Zo = sqrt(L/C) and beta = w*sqrt(LC)
slot_per_unit_C(n) = slot_per_unit_L(n)/slot_Zo(n)^2; % Calculated from Zo = sqrt(L/C) and beta = w*sqrt(LC)
slot_G_rad(n) = 2.63e-3; % Roughly 2.63*10^(-1.98) is the ceiling for this value before IL goes above 0 dB
slot_per_unit_R(n) = 2.63e-3*(w_slot/d(n))*sqrt(freq(n)*10^(-9)); % Calculated from formula in Dr. Ruyle's dissertation
slot_per_unit_G(n) = 2*pi*freq(n)*skin_depth(n)*e0*p(n)^2*4*(h_copper/w_slot) + slot_G_rad(n); % Calculated from formula in Dr. Ruyle's dissertation
slot_gamma(n) = sqrt((slot_per_unit_R(n)+j*w(n)*slot_per_unit_L(n))*(slot_per_unit_G(n)+j*w(n)*slot_per_unit_C(n))); % Propagation constant from standard formula
slot_ABCD(:,:,n) = [cosh(slot_gamma(n)*l_slot1) slot_Zo(n)*sinh(slot_gamma(n)*l_slot1); 1/slot_Zo(n)*sinh(slot_gamma(n)*l_slot1) cosh(slot_gamma(n)*l_slot1)]; %Find ABCD parameters

%%% Chip Calculations %%%
chip_ABCD(:,:,n) = inv(slot_ABCD(:,:,n))*inv(coax1_ABCD(:,:,n))*MeasuredABCD(:,:,n)*inv(coax2_ABCD(:,:,n)); % De-embed the chip by multiplying inverse ABCD matrices
chip_Sparam = abcd2s(chip_ABCD,Zo); % Convert to S parameters from ABCD parameters
chip_S21(n) = chip_Sparam(2,1,n); % Extract S21

%%% Impedance Calculations %%%
S11(n) = MeasuredS(1,1,n);
MeasuredZin(n) = Zo*(1+S11(n))/(1-S11(n));

%%% Slot + Coax %%%
% CoaxSlotCascaded(:,:,n) = coax1_ABCD(:,:,n)*slot_ABCD(:,:,n);
% CalculatedZin(n) = CoaxSlotCascaded(1,2,n)/CoaxSlotCascaded(2,2,n);

%%% Coax Only %%%
CalculatedZin(n) = coax1_ABCD(1,2,n)/coax1_ABCD(2,2,n);
end



figure(2)
plot(freq,imag(MeasuredZin),'r--',freq,imag(CalculatedZin))
legend('Measured','Calculated');
title('Measured Imag(Zin)');

figure(1)
plot(freq,20*log10(abs(chip_S21)))
xlabel('Frequency (Hz)')
ylabel('Insertion Loss (dB)')
title('Insertion Loss')