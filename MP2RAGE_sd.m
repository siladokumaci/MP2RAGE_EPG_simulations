function [F0_arr,Zn_arr,finit_arr,all_F0,all_Zn,all_finit,PE_steps]=MP2RAGE_sd(TR,T1,T2,PE_steps,number_of_reps,PD)
% EPG Simulations for MP2RAGE 
% Author: Ayse Sila Dokumaci (ayse.dokumaci@kcl.ac.uk) 28.02.2020
% This function requires EPG_GRE2 and EPG_GRE3 functions modified from the original EPG_GRE.m and the EPG-X codes available at https://github.com/mriphysics/EPG-X
% The references for the EPG-X codes above are:
% Malik S. (2017, August 8). mriphysics/EPG-X: First public version (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.840023
% Malik S, Teixeira RPAG, Hajnal JV. Extended phase graph formalism for systems with magnetization transfer and exchange. Magn Reson Med. 2018;80:767-779. 

% Transverse magnetization is not assumed to be zero but is fed into the next run as an input.
% The function can be called by run_mp2rage_sd for WM/GM/CSF.
% The inputs:
% TR is the repetition time of the GRE blocks.
% number_of_reps is the total number of repetitions of the MP2RAGE to reach steady state.
% T1, T2, PD are the relaxation parameters and proton density.
% PE_steps would for instance be 192 for 256 "partitions" (PE-steps in the innermost acquisition loop) with a partial Fourier factor of 6/8.

% There are some lines that might need modification for different simulations which are: 
% line 39 for b1_scale (B1+); 
% line 40 for initial phase in the RF phase cycling; 
% line 42 for k-space centre; 
% lines 44-48 for TR_MP2RAGE/TI2/TI1/FA2/FA1;  
% line 56 or 57 should be commented depending on partial Fourier; 
% lines 148-150-152 if PD values different than 0.69, 0.82 or 1 are used for WM/GM/CSF,
% and line 164 is for the name of the mat file saved from the simulations.

% The outputs: It will produce files ending with the B1+ scale (for instance 1P2 for 120% B1+).
% The outputs all_F0 will be the signals. For instance, for FA1 and FA2 changing from 1 to 10 degrees independently, 
% all_F0 will be a 1 x 100 cell array. For the current case, each cell will be of size 197 x 10 where 
% 197 corresponds to 192 + 5. 192 is resulting from 256 "partitions" * 6/8 with the partial Fourier and 5
% is the result of storing the flip angles scaled by the B1+ and the TI1/TI2/TR_MP2RAGE times in the first 5 columns. 
% Partitions is the number of PE steps in the innermost acquisition loop.
% The signal evolution for the TI2 signal is stored in the last column starting from the 6th element
% (so the k-space centre for this case will be (192/3 + 1) + 5 = 70th element of that column) 
% and the TI1 signal is stored in the penultimate column.

% The code can be written more efficiently but I hope you find this version useful :)
tic
inv_eff_asd = 1;
for b1_scale = 0.5:0.1:1.4 % the loop for B1+ => this one would be 50% to 140% with steps of 10%
        initial_phase = 50;% this was the default value in the sequence implementation
        mp2rage_phase=RF_phase_cycle(PE_steps,initial_phase)';
        kspace_centre =  65; % 65 for the case with 256 partitions and 6/8 partial Fourier => 256*6/8/3 + 1 = 65
        count = 0;% keeping track of the loops
        for MP2RAGE_TR = 4000% TR_MP2RAGE which is set to 4000 ms 
            for TI2 = 2280 % Second Inversion Time
                for TI1 = 650 % First Inversion Time
                    for fa2_scale= 1:10 % Second Flip Angle [degrees]
                        for fa1_scale= 1:10 % First Flip Angle [degrees] 
                            fa1 = (ones(1,PE_steps)*fa1_scale*b1_scale/180*pi)';
                            fa2 = (ones(1,PE_steps)*fa2_scale*b1_scale/180*pi)';
                            if(TI2-TI1>PE_steps*TR)
                                count = count + 1;
                                
                                prep1=struct;
                                prep1.flip=pi*inv_eff_asd;% added for inversion efficiency (eff)
                                %prep1.t_delay=TI1-TR*(PE_steps/2);% TA in the MP2RAGE paper if the k-space centre is sampled in the middle of the PE_steps
                                prep1.t_delay=TI1-TR*(kspace_centre - 1);% TA for partial Fourier acquisitions
                                
                                prep2=struct;
                                prep2.flip=0;
                                prep2.t_delay=TI2-(TI1+PE_steps*TR);% TB in the MP2RAGE paper 
                                
                                % added for diffusion on 22.03.2021
                                d=struct;
                                d.G = [-5.9 4.5 8.1*sqrt(2)]; % mT/m
                                d.tau = [1.4 6.6 3]; %ms
                                d.D = 2.3e-9; %<- diffusion coeff in m^2/s
                                
                                %TC = MP2RAGE_TR-(TI2+TR*PE_steps/2);% TC if the k-space centre is sampled in the middle of the PE_steps 
                                TC = MP2RAGE_TR-(TI2+TR*(PE_steps-(kspace_centre-1)));% TC for partial Fourier acquisitions 
                                
                                % Applying EPG_GRE (the first image in the first run)
                                [F0,~,Zn,F] = EPG_GRE2(fa1,mp2rage_phase,TR,T1,T2,'kmax',inf,'prep',prep1,'zinit',PD,'diff',d);% diffusion was added on 22.03.2021
                                finit = F(:,size(F,2));% The last column of F as an input for the second image
                                
                                F0_arr(1:5,1) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];% recording the FA1, FA2, TI1, TI2, and MP2RAGE_TR in the first 5 columns
                                Zn_arr(1:5,1) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                F0_arr(1:(5+size(F0,1)),1) = cat(1,F0_arr(1:5,1),F0);
                                Zn_arr(1:(5+size(Zn(1,:),2)),1) = cat(1,Zn_arr(1:5,1),Zn(1,:)');
                                
                                % Applying EPG_GRE3 (the second image in the first run)
                                [F0,~,Zn,F] = EPG_GRE3(fa2,mp2rage_phase,TR,T1,T2,'kmax',inf,'prep',prep2,'finit',finit,'zinit',PD,'diff',d);% diffusion was added on 22.03.2021
                                finit = F(:,size(F,2));% The last column of F for relaxation
                                
                                F0_arr(1:5,2) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];% recording the B1+ scaled FA1, B1+ scaled FA2, TI1, TI2, and MP2RAGE_TR in the first 5 columns
                                Zn_arr(1:5,2) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                F0_arr(1:(5+size(F0,1)),2) = cat(1,F0_arr(1:5,2),F0);
                                Zn_arr(1:(5+size(Zn(1,:),2)),2) = cat(1,Zn_arr(1:5,2),Zn(1,:)');
                                
                                ridx = length(finit);
                                E1_3 = exp(-TC/T1);
                                E2_3 = exp(-TC/T2);
                                %%% regrowth
                                b3 = zeros([ridx 1]);
                                b3(3) = PD*(1-E1_3);% just applies to Z0 % PD factor was added by Sila Dokumaci 
                                E_3 = spdiags(repmat([E2_3 E2_3 E1_3],[1 ridx])',0,ridx,ridx);
                                for ii=1:ridx
                                    finit(ii) = E_3(ii,ii)*finit(ii)+b3(ii);
                                end
                                
                                finit_arr(1:5,1) = {fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR};
                                finit_arr(1:(5+size(finit,1)),1) = cat(1,finit_arr(1:5,1),num2cell(finit));
                                
                                for i=1:number_of_reps-1
                                    
                                    finit = cell2mat(finit_arr(6:end,i));% again the first 5 elements are fa1_scale*b1_scale,fa2_scale*b1_scale,TI1, TI2, and MP2RAGE_TR
                                    
                                    % The first image starting from repetition 2 on
                                    [F0,~,Zn,F] = EPG_GRE3(fa1,mp2rage_phase,TR,T1,T2,'kmax',inf,'prep',prep1,'finit',finit,'zinit',PD,'diff',d);% diffusion was added on 22.03.2021
                                    finit = F(:,size(F,2));% The last column of F as an input for the second image
                                    
                                    F0_arr(1:5,2*i+1) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                    Zn_arr(1:5,2*i+1) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                    F0_arr(1:(5+size(F0,1)),2*i+1) = cat(1,F0_arr(1:5,2*i+1),F0);
                                    Zn_arr(1:(5+size(Zn(1,:),2)),2*i+1) = cat(1,Zn_arr(1:5,2*i+1),Zn(1,:)');
                                    
                                    % The second image starting from repetition 2 on
                                    [F0,~,Zn,F] = EPG_GRE3(fa2,mp2rage_phase,TR,T1,T2,'kmax',inf,'prep',prep2,'finit',finit,'zinit',PD,'diff',d);% diffusion was added on 22.03.2021
                                    finit = F(:,size(F,2));% The last column of F for relaxation
                                    
                                    F0_arr(1:5,2*i+2) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                    Zn_arr(1:5,2*i+2) = [fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR];
                                    F0_arr(1:(5+size(F0,1)),2*i+2) = cat(1,F0_arr(1:5,2*i+2),F0);
                                    Zn_arr(1:(5+size(Zn(1,:),2)),2*i+2) = cat(1,Zn_arr(1:5,2*i+2),Zn(1,:)');
                                    
                                    ridx = length(finit);
                                    E1_3 = exp(-TC/T1);
                                    E2_3 = exp(-TC/T2);
                                    %%% regrowth
                                    b3 = zeros([ridx 1]);
                                    b3(3) = PD*(1-E1_3);% just applies to Z0 % PD factor was added by Sila Dokumaci
                                    E_3 = spdiags(repmat([E2_3 E2_3 E1_3],[1 ridx])',0,ridx,ridx);
                                    for ii=1:ridx
                                        finit(ii) = E_3(ii,ii)*finit(ii)+b3(ii);
                                    end
                                    finit_arr(1:5,i+1) = {fa1_scale*b1_scale; fa2_scale*b1_scale; TI1; TI2; MP2RAGE_TR};
                                    finit_arr(1:(5+size(finit,1)),i+1) = cat(1,finit_arr(1:5,i+1),num2cell(finit));
                                end
                                all_finit{count} = finit_arr;
                                all_Zn{count} = Zn_arr;
                                all_F0{count} = F0_arr;
                            end
                        end
                    end
                end
            end
        end
        if(PD == 0.69)
            tissue_type = 'WM';
        elseif(PD == 0.82)
            tissue_type = 'GM';
        elseif(PD == 1)
            tissue_type = 'CSF';
        end
       
        b1_scale_str = num2str(b1_scale);
        if(b1_scale == 1)
            b1_scale_str = '1P0';
        else
            b1_scale_str = strcat(b1_scale_str(1),'P',b1_scale_str(3));
        end
        
        clearvars -except all_F0 all_Zn tissue_type b1_scale_str TR T1 T2 PE_steps number_of_reps PD TI1 TI2 inv_eff_asd
        res_filename = strcat('Sims_650_2280_TR4000_',tissue_type,'_',b1_scale_str);
        save(res_filename,'-v7.3')%,'-nocompression')
end
toc
end