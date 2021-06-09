Online supplement for the paper titled "Maximum Expected Entropy Transformed Latin Hypercube Designs"
====

## Requirement: Matlab R2016b (or newer version)
## Usage:

1. MeeBTLHD.m generates Mee-BTLHDs based on the Gaussian correlation function for run size n and dimension d, with the given set of values for correlation parameter theta, the given set of values for nugget and the starting design MmLHD generated from the R package SLHD.
   ```matlab
   D = MeeBTLHD(n,d,theta,nugget,MmLHD);
   ```
2. MeeNTLHD.m generates Mee-NTLHDs based on the Gaussian correlation function for run size n and dimension d, with the given set of values for correlation parameter theta, the given set of values for nugget and the starting design MmLHD generated from the R package SLHD.
   ```matlab
   D = MeeNTLHD(n,d,theta,nugget,MmLHD);
   ```
3. plot_for_Figure3.m reproduces Figure 3.

   Step 1: load('Designs for Gaussian.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Gaussian correlation function and four other designs: MaxPro, MmLHD, uniform and GMLHD

   Step 2: run plot_for_Figure3.m; % reproduces Figure 3

4. plot_for_Figure5.m reproduces Figure 5.

   Step 1: load('Designs for Gaussian.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Gaussian correlation function and four other designs: MaxPro, MmLHD, uniform and GMLHD
   
   Step 2: run plot_for_Figure5.m; % reproduces Figure 5
   
5. plot_for_Figure6.m reproduces Figure 6.

   Step 1: load('Designs for Gaussian.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Gaussian correlation function and four other designs: MaxPro, MmLHD, uniform and GMLHD
   
   Step 2: run plot_for_Figure6.m; % reproduces Figure 6

6. FIRMeeBTLHD.m generates FIRMee-BTLHDs based on the Gaussian correlation function for run size n and dimension d, with the given set of values for correlation parameter theta, the given set of values for nugget and the starting design MmLHD generated from the R package SLHD.
   ```matlab  
   D = FIRMeeBTLHD(n,d,theta,nugget,MmLHD);
   ```
7. plot_for_FigureC2.m reproduces Figure C2.

   Step 1: load('Designs for FIRMee.mat'); % load FIRMee-BTLHD generated based on the Gaussian correlation function, four other designs: MaxPro, MmLHD, uniform and GMLHD, the value of theta used to construct FIRMee-BTLHD 'ctheta' and the value of nugget used to construct FIRMee-BTLHD 'cnugget'

   Step 2: run plot_for_FigureC2.m; % reproduces Figure C2
   
8. MeeBTLHD_matern.m generates Mee-BTLHDs based on the Matern correlation function with smoothness parameter 5/2 for run size n and dimension d, with the given set of values for correlation parameter theta, the given set of values for nugget and the starting design MmLHD generated from the R package SLHD.
   ```matlab  
   D = MeeBTLHD_matern(n,d,theta,nugget,MmLHD);
   ```
9. MeeNTLHD_matern.m generates Mee-NTLHDs based on the Matern correlation function with smoothness parameter 5/2 for run size n and dimension d, with the given set of values for correlation parameter theta, the given set of values for nugget and the starting design MmLHD generated from the R package SLHD.
   ```matlab  
   D = MeeNTLHD_matern(n,d,theta,nugget,MmLHD);
   ```
10. plot_for_FigureD1.m reproduces Figure D1.

    Step 1: load('Designs for Matern.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Matern correlation function with smoothness parameter 5/2 and four other designs: MaxPro, MmLHD, uniform and GMLHD
   
    Step 2: run plot_for_FigureD1.m; % reproduces Figure D1
   
11. plot_for_FigureD2.m reproduces Figure D2.

    Step 1: load('Designs for Matern.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Matern correlation function with smoothness parameter 5/2 and four other designs: MaxPro, MmLHD, uniform and GMLHD
   
    Step 2: run plot_for_FigureD2.m; % reproduces Figure D2
    
12. plot_for_FigureD3.m reproduces Figure D3.

    Step 1: load('Designs for Matern.mat'); % load Mee-BTLHD, Mee-NTLHD generated based on the Matern correlation function with smoothness parameter 5/2 and four other designs: MaxPro, MmLHD, uniform and GMLHD
   
    Step 2: run plot_for_FigureD3.m; % reproduces Figure D3