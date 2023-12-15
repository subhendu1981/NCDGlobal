# Python code to model non-cyanobacterial diazotrophs
# associated with sinking marine particles

# File: parameters.py

import numpy as np
from scipy import interpolate

pi = 3.141592653589793

init_rad = 0.25 # 0.15 # cm # radius of particle

x_B = 5*10**(-8) # mug C/cell # cellular C 

f_PS = 0.238
f_PP = 0.5

# Initial conditions 
C_init = 2.6*10**8 # 
P_init = 1.6*10**8 # 
G_init = 50 # mug/L 
A_init = 5 # mug/L 
B_init = 10**10 # (cells/L)
O2_init = 200 # mumol O2/L
NO3_init = 10*62 # mug NO3/L 
SO4_init = 28*10**3 # mumol SO4/L 

# Boundary conditions 
G_bdry = G_init                                                          
A_bdry = A_init                                                         
O2_bdry = O2_init     
NO3_bdry = NO3_init      
SO4_bdry = SO4_init     

# other parameters 
Q_CN_B = 3.7 # mug C/mug N   

C_m = 0.0412 # mug % Cram et al. 2018
rho_ref = 1.23 # g/ml % Cram et al. 2018
r_ref = 20*10**(-4) # cm % Cram et al. 2018
rho_particle = 1.033 # g/ml % Cram et al. 2018

h_C = 2.25*10**(-5) # Ploug et ak. 1999 % (mug glucose/cell/d) 
h_P = 6.1*10**(-5) # Ploug et ak. 1999 % (mug amino acid/cell/d) (h_P is 2 times than h_c; Steen1 and Arnosti 2013) (19-5232*120.5*10^(-9) mug/cell/d; Grossart and Simon, 1998)
K_h_C = 1.8*10**4 # (mug glucose/L) Huston and Deming (2002) 
A_C = h_C/K_h_C # L/cell/day
K_h_P = 3.6*10**3 # (mug amino acid/L) Huston and Deming (2002) 
A_P = h_P/K_h_P # L/cell/day

J_max_G = 7*10**(-7) # 1.5768e-06; % 26.28*10^(-9)*24/0.4; % mug glucose/cell/d (Azua et al. 2007)       %max glucose uptake rate
K_G = 252.75 # 101.1/f_G_C; %  (4.5-455.4 mug Glu/L; Ayo et al. 2001)
alpha_G = J_max_G/K_G # L/d                     %affinity for glucose
J_max_A = 4.42*10**(-7) # 153*10^(-12)*24*120.5; % mug Amino Acid/cell/d (Ayo et al. 2001)    %max Amino Acid uptake rate 
K_A = 63.62 # (40.36-63.62 mug AA/L; Ayo et al. 2001) 
alpha_A = J_max_A/K_A # L/d                   %affinity for Amino Acid

rB1 = (3/(4*pi)*(x_B*10**9/133.754)**(1/0.438))**(1/3)*10**(-4) # cm      %(radius) (0.29*10^(-6) m)
dB = 2*rB1*10**4 # micro m

M1 = np.loadtxt("Size_vs_NO3_uptake.csv", delimiter =",") # Atkins et al. 2016
diam = M1[:,0] # micro m
uptk = 24*M1[:,1]*x_B*62/12 # mumol NO3/cell/d % 
J_max_NO3_all = interpolate.interp1d(diam, uptk, fill_value = "extrapolate")
J_max_NO3 = J_max_NO3_all(dB) # mumol NO3/cell/d 

M2 = np.loadtxt("Size_vs_halfsaturation_NO3_uptake.csv", delimiter =",") # Atkins et al. 2016
diam1 = M2[:,0] # micro m
hs = M2[:,1]*62 # mumol NO3/L % 
hs1_all = interpolate.interp1d(diam1, hs, fill_value = "extrapolate")
hs1 = hs1_all(dB) # mumol NO3/cell/d 
alpha_NO3 = J_max_NO3/hs1 # L/d/cell 

J_max_SO4 = 5*10**(-10) # mumolN/cell/d, Kondo et al. 2004) % max SO4 uptake rate

M_N2 = 5.77*10**(-8) # Bentzon-Tilia et al. 2015 (1.55*10^(-8); Paerl et al. (2018)); % mug N/cell/d   %(max N fixation rate)

f_G_C = 0.4 # Lopez-Fernandez et al. Biogeosciences 2013
f_A_C = 0.445 # Lee and Cronin 1982 (0.49-Lopez-Fernandez_Biogeosciences_2013)
f_A_N = 0.125 # Lee and Cronin 1982

R_E = 0.6
R_B = 0.05+R_E # 0.05+0.6; % 0.21+1.2; % 1/d % Mislan 2014
R_G = 0.23 # mug C/mug C 
R_A = 0.23 # mug C/mug C % 1*12/14/param.Q_CN_B % Flynn (2005)
R_N2 = 0.4 # mug C/mug C
R_NO3 = 0.4 # 1; % mug C/mug C % 1.7*12/14/param.Q_CN_B % Flynn (2005)
R_SO4 = 0.6 # 2; % mug C/mug C

rho_CO = 10 # mug C/mumol O2 % Ploug et al. 1999
rho_CNO3 = 12.5/62 # mug C/mug NO3 % Paulmier et al. 2009
rho_CSO4 = 20 # mug C/mumol SO4 % Henze et al., 2008


#%-- Diffusive oxygen inflow --%
rB = rB1
vol_cell = 4/3*pi*rB**3*10**(-3) # L (dm^3)
rc = rB-(10+6)*10**(-7) # cm                                            %(thickness of cytoplasm)             
D_O2_water = 2.12*10**(-5) # McCabe and Laurent 1975  (1.38-2.12*10^(-5); % cm^2/s Broecker & Peng 1974)  (2*10^(-5) cm^2/s; Picioreanu et al 1997)                                            %(diff coeff of O2 in water) 
D_O2_particle = 0.95*D_O2_water # cm^2/s (Ploug and Passow, 2007) ((2 orders of lower) Brzezinski et al 1997)            %(diff coeff of O2 in aggregate) 
D_NO3_water = 1.6*10**(-5) #  cm^2/s      (Li and Gregory 1974)     (1.7*10^(-5); Picioreanu et al 1997)             %(diff coeff of NO3 in water) 
D_NO3_particle = 0.95*D_NO3_water # cm^2/s (Ploug and Passow, 2007)          %(diff coeff of NO3 in aggregate) 
D_SO4_water = 8.9*10**(-6) #  cm^2/s      (Li and Gregory 1974)               %(diff coeff of SO4 in water) 
D_SO4_particle = 0.8*D_SO4_water # cm^2/s Williamson and McCarty 1976) (0.95; Ploug and Passow, 2007)          %(diff coeff of SO4 in aggregate) 

eps_m = 7.9*10**(-4) # unitless (Inomura et al. 2017)                       %(diff. of O2 in cytoplasm relative to water)
Lm = 6*10**(-7) # cm                                                     %(thickness of plasma membrane)
K_O2 = D_O2_particle*eps_m*(rc+Lm)/(eps_m*rc+Lm) # cm^2/s                         %(apparent diff of O2 in cell) 
D_O2_cell = 4*pi*rB*10**(-1)*K_O2*10**(-2)*24*60*60 # dm^3/d = L/d/cell             %(oxygen flux in cell)
 
m_B = 0.1 # 0.1; % 0.2; % Miki and Yamamura 2005

D_C = 0 #                                                  %(diff of carbohydrate in aggregate)                
D_P = 0  
D_G = 6.0*10**(-6)*60*60*24 # Stein 1990% (5.62*10^(-6)*60*60*24; Zhang 2005); % 76*60*60*24*10^(-8); % 10*60*60*24*10^(-8); % cm^2/d (Nath et al. 2010)                      %(diff of glucose in aggregate)                
D_A = 6.0*10**(-6)*60*60*24 # Stein 1990% (5.62*10^(-6)*60*60*24; Zhang 2005); % 76*60*60*24*10^(-8); % 10*60*60*24*10^(-8); % cm^2/d (Nath et al. 2010)                      %(diff of amino acid in aggregate)  
D_B = 0 
D_OP = D_O2_particle*60*60*24 # cm^2/d  
D_NP = D_NO3_particle*60*60*24 # cm^2/d  
D_SP = D_SO4_particle*60*60*24 # cm^2/d  
D_O2 = D_O2_water*60*60*24 # cm^2/d

# --- Q10 parameters --- %
Q10_h = 2 # hydrolysis
Q10_u = 1.5 # uptake
Q10_r = 2 # respiration

# ---- reference temperatures ---- %
ref_temp_hC = 17 # Ploug et al. 1999
ref_temp_AC = -1 # Huston and Deming 2002
ref_temp_hP = 17 # Ploug et al. 1999
ref_temp_AP = -1 # Huston and Deming 2002
ref_temp_JmaxG = 20 # Azua et al. 2007 (room temp)
ref_temp_alphaG = 20 # Azua et al. 2007 (room temp)
ref_temp_JmaxA = 13 
ref_temp_alphaA = 13 
ref_temp_JmaxNO3 = 4 # Fouilland et al. 2007
ref_temp_alphaNO3 = 35 # Fagerbakke et al. 1996
ref_temp_JmaxSO4 = 18 # Kondo et al. 2004
ref_temp_MN2 = 20 # Bentzon-Tilia et al. 2015 (room temp)
ref_temp_RB = 20 # Mislan et al. 2014
ref_temp_RE = 20 # Mislan et al. 2014
ref_temp_RG = 20 # assumption
ref_temp_RA = 20 # assumption
ref_temp_RN2 = 28 # Grosskopf and LaRoche 2012
ref_temp_RNO3 = 20 # assumption
ref_temp_RSO4 = 20 # assumption
ref_temp_DO2 = 25 # McCabe and Laurent 1975
ref_temp_DNO3 = 18 # Li and Gregory 2974
ref_temp_DM = 20 # Grigoriev and Meylikhov 1991

K = 273.15

a = 3.5481*10**(-4) # Jackson et al. 1997

alpha = 2.3 # 2.81; % Durkin et al. 2015;

threshold = 10**(-15)

# temperature
temp_current = 15

# viscosity
dataa4 = np.loadtxt("viscocity_temperature_Jumars_1993.csv", delimiter =",") 
vis_temp_data = dataa4[:,0] # C
vis_data = dataa4[:,1] # g/cm/s

ref_vis_DO2_all = interpolate.interp1d(vis_temp_data, vis_data, fill_value = "extrapolate")
ref_vis_DO2 = ref_vis_DO2_all(ref_temp_DO2) 

ref_vis_DNO3_all = interpolate.interp1d(vis_temp_data, vis_data, fill_value = "extrapolate")
ref_vis_DNO3 = ref_vis_DNO3_all(ref_temp_DNO3) 

ref_vis_DM_all = interpolate.interp1d(vis_temp_data, vis_data, fill_value = "extrapolate")
ref_vis_DM = ref_vis_DM_all(ref_temp_DM) 

# sinking speed
eta = 0.26 #         % sinking speed exponent 
cw  = 1.0895*10**4 #         % sinking speed coefficient
