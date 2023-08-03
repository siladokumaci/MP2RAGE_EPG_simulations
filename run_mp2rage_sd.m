% References for the T1 and T2 values below:
% Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele PF, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. Neuroimage. 2010;49:1271-1281.
% Rooney WD, Johnson G, Li X, et al. Magnetic field and tissue dependencies of human brain longitudinal H2O relaxation in vivo. Magn Reson Med. 2007;57:308-318. 
% Yacoub E, Duong TQ, van De Moortele PF, et al. Spin-echo fMRI in humans using high spatial resolutions and high magnetic fields. Magn Reson Med. 2003;49(4):655-664.

T1_mp2rage_script_arr = [1220 2132 3350];% T1 values for WM/GM/CSF
T2_mp2rage_script_arr = [45.9 55 1000];% T2 values for WM/GM/CSF 
PD_mp2rage_script_arr =[0.69 0.82 1];% PD values for WM/GM/CSF

% 7.9 ms was the TR_GRE value that we used in the last protocol which had 256 partitions (PE steps in the innermost acquisition loop) 
% with a partial Fourier factor of 6/8 which led to 192 steps here.
parfor mp2rage_script_index = 1:3 
   MP2RAGE_sd(7.9,T1_mp2rage_script_arr(mp2rage_script_index),T2_mp2rage_script_arr(mp2rage_script_index),192,5,PD_mp2rage_script_arr(mp2rage_script_index));
end
