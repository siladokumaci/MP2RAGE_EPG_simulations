# MP2RAGE_EPG_simulations
Extended Phase Graph (EPG) Simulations of the MP2RAGE Sequence including T2 relaxation, diffusion, radio frequency spoiling, and B1+ variability
Transverse magnetization is not assumed to be zero but is fed into the next run as an input.

The main function **MP2RAGE_sd** provides simulations of the MP2RAGE sequence signal behaviour using EPG Simulations.
MP2RAGE sequence is described in this reference: Marques JP, Kober T, Krueger G, van der Zwaag W, Van de Moortele PF, Gruetter R. MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. Neuroimage. 2010;49:1271-1281.  

The MP2RAGE_sd function requires **EPG_GRE2** and **EPG_GRE3** functions modified from the original EPG_GRE function and the **EPG-X codes** which are available at [https://github.com/mriphysics/EPG-X]

The MP2RAGE_sd function can be called by **run_mp2rage_sd** function for WM/GM/CSF.

The **inputs**:
- TR is the repetition time of the GRE blocks.
- number_of_reps is the total number of repetitions of the MP2RAGE to reach steady state.
- T1, T2, PD are the relaxation parameters and proton density.
- PE_steps would for instance be 192 for 256 "partitions" (PE-steps in the innermost acquisition loop) with a partial Fourier factor of 6/8.

There are **some lines that might need modification** for different simulations which are: 
- line 38 for b1_scale (B1+); 
- line 39 for initial phase in the RF phase cycling; 
- line 41 for k-space centre; 
- lines 43-47 for TR_MP2RAGE/TI2/TI1/FA2/FA1;  
- line 55 or 56 should be commented depending on partial Fourier; 
- lines 147-149-151 if PD values different than 0.69, 0.82 or 1 are used for WM/GM/CSF;
- and line 163 is for the name of the mat file saved from the simulations.

The **outputs**:
- It will produce files ending with the B1+ scale (for instance 1P2 for 120% B1+).
- The outputs all_F0 will be the signals. For instance, for FA1 and FA2 changing from 1 to 10 degrees independently, 
- all_F0 will be a 1 x 100 cell array. For the current case, each cell will be of size 197 x 10 where 
197 corresponds to 192 + 5. 192 is resulting from 256 "partitions" * 6/8 with the partial Fourier and 5
is the result of storing the flip angles scaled by the B1+ and the TI1/TI2/TR_MP2RAGE times in the first 5 columns. 
- Partitions is the number of PE steps in the innermost acquisition loop.
- The signal evolution for the TI2 signal is stored in the last column starting from the 6th element 
(so the k-space centre for this case will be (192/3 + 1) + 5 = 70th element of that column) 
and the TI1 signal is stored in the penultimate column.

The code can be written more efficiently but I hope you find it useful :)

**References**
The EPG-X codes are distributed under the MIT license. If you find them useful, please cite the publication below or the code itself:
1. Malik S, Teixeira RPAG, Hajnal JV. Extended phase graph formalism for systems with magnetization transfer and exchange. Magn Reson Med. 2018;80:767-779. 
2. Malik S. (2017, August 8). mriphysics/EPG-X: First public version (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.840023

The references for the EPG theory are:
1. Hennig J. Echoes - how to generate, recognize, use or avoid them in MR-imaging sequences. Part I. Concepts Magn Reson. 1991;3:125-143. 
2. Hennig J. Echoes - how to generate, recognize, use or avoid them in MR-imaging sequences. Part II. Concepts Magn Reson. 1991;3:179-192. 
3. Weigel M. Extended phase graphs: dephasing, RF pulses, and echoes - pure and simple. J Magn Reson Imaging. 2015;41:266-295.

If you find the MP2RAGE simulation code useful, please cite the abstracts below along with the EPG-X code references above. The manuscript has been submitted for publication.  
1. Dokumaci AS. Simultaneous optimisation of MP2RAGE UNI and FLAWS brain images at 7T using Extended Phase Graph (EPG) Simulations. In Proceedings of the 29th Annual Meeting of ISMRM. 2021. Abstract 445.
2. Dokumaci AS. Simultaneous optimisation of MP2RAGE UNI & FLAWS images at 7T with Extended Phase Graph Simulations. In Proceedings of the 27th Annual Meeting of the Organization for Human Brain Mapping. 2021. p. 39.

