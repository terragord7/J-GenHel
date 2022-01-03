# DESCRIPTION
# This file contains the helicopter parameters.. Always re-run DI_design.jl
# after modifying the aircraft parameters.

# ------------------------- MISCELLANEOUS CONSTANTS ----------------------------

# conversion constants
D2R=pi/180
R2D=180/pi
g=32.17
FtLb_s2Hp=1/550
fps2kts=1/1.69
kts2fps=1.69
# air density [slug/ft³]
rhoSLSTD=0.002378
RHO=rhoSLSTD
# speed of sound [ft/s]
VSOUND=1135.0
# acceleration of gravity [ft/sec^2]
G=32.17

# -------------------------------- MAIN ROTOR ----------------------------------

# radius [ft]
R=26.83
# number of blade segments
NSEG=10
# number of blades
NB=4
# angular speed [rad/s]
OMEGA=27.0
# hinge offset [ft]
HOFFSET=1.25
# spar length (distance from hinge to start of rotor blade) [ft]
SPAR=2.25
# lift slope of airfoil [1/rad]
LIFTSLOPE=6.5
# zero lift angle of attack
ALPHA0=-2*pi/180.0
# profile drag coefficient
CD0=0.008
# hinge offset
e=HOFFSET/R
eprime=SPAR/R
# blade element locations
rseg=zeros(NSEG)
delseg=zeros(NSEG)
rseg[1]=sqrt((1-(e+eprime)^2)/(2*NSEG)+(e+eprime)^2)-e
delseg[1]=sqrt((rseg[1]+e)^2+(1-(e+eprime)^2)/(2*NSEG))-
          sqrt((rseg[1]+e)^2-(1-(e+eprime)^2)/(2*NSEG))
for iseg=2:NSEG
    rseg[iseg]=sqrt((1-(e+eprime)^2)/(NSEG)+(e+rseg[iseg-1])^2)-e
    delseg[iseg]=sqrt((rseg[iseg]+e)^2+(1-(e+eprime)^2)/(2*NSEG))-
                 sqrt((rseg[iseg]+e)^2-(1-(e+eprime)^2)/(2*NSEG))
end
RSEG = rseg
DELSEG = delseg
# dynamic twist mode shape
DynTwistMode = zeros(NSEG)
for iseg=1:NSEG
    DynTwistMode[iseg] = 0.28+(0.72*sin(0.5*pi*rseg[iseg]))
end
TauDT=0.05
# rigid blade prescribed twist
TWISTTABR = collect(0.0:0.05:1.0)
TWISTTABTHET= [0.0,0.0,0.0,0.0,-0.15,-0.95,-1.80,-2.75,-3.55,-4.40,-5.30,-6.15,
               -7.10,-7.90,-8.80,-9.65,-10.30,-10.75,-12.30,-13.10,-10.9]
# first mass moment of blade [sl-ft]
MBETA=86.7
# moment of intertia of blade about flap hinge [sl-ft²]
IBETA=1512.6
# blade weight [lbs]
WBLADE=256.9
# lag damper geometry [in]
ALD=0.227
BLD=3.242
CLD=12.04
DLD=10.0102
RLD=6.898
THETALDGEO=0.0
# shaft incidence
ITHSH=-3.0
IPHSH=0.0
Tshaft=[cos(ITHSH*D2R)                 0              -sin(ITHSH*D2R)
        sin(ITHSH*D2R)*sin(IPHSH*D2R)  cos(IPHSH*D2R)  cos(ITHSH*D2R)*sin(IPHSH*D2R)
        sin(ITHSH*D2R)*cos(IPHSH*D2R) -sin(IPHSH*D2R)  cos(ITHSH*D2R)*cos(IPHSH*D2R)]
# swashplate phasing angle [rad]
DELSP=-9.7*D2R
# δ₃ angle [rad]
DELTA3=0.0
# lag damping parameter [lbs/(ft/s)]
CLAG=1500.0
# lag damper stiffness [ft-lbs/rad]
KLAG=15060.0
# lag angle for 0 spring force [rad]
ZETA0=7*pi/180
# notional flap damper and spring
CFLAP=0.0
KFLAP=0.0
# blade tip chord [ft]
CHORDT=1.73
# blade root chord [ft]
CHORDR=1.73
# chord distribution
SEGOUT = zeros(NSEG)
SEGIN = zeros(NSEG)
CHORD = zeros(NSEG)
for iseg=1:NSEG
    SEGOUT[iseg]=sqrt(iseg*(1-(e + eprime)^2)/NSEG+(e + eprime)^2)
    SEGIN[iseg]=sqrt((iseg-1)*(1-(e + eprime)^2)/NSEG+(e + eprime)^2)
    CHORD[iseg] = (SEGOUT[iseg] + SEGIN[iseg] - 2*(e + eprime))*
                  (0.5/(1-e-eprime))*(CHORDT-CHORDR)+CHORDR;
end
# NACA 0012 airfloil data
fileIn = matopen("naca0012.mat")
ALTAB=read(fileIn, "ALTAB")
RETAB=read(fileIn, "RETAB")
CLTAB=read(fileIn, "CLTAB")
CDTAB=read(fileIn, "CDTAB")
close(fileIn)
MRAOATAB=ALTAB
MRCLTAB=CLTAB[:,5]
MRCDTAB=CDTAB[:,5]
# GenHel rotor airfoil data
fileIn = matopen("GHRotorData.mat")
CDR8UT=read(fileIn, "CDR8UT")
CDR8BT=read(fileIn, "CDR8BT")
CLR0BT=read(fileIn, "CLR0BT")
CLR8UT=read(fileIn, "CLR8UT")
CLR8BT=read(fileIn, "CLR8BT")
CDR0UT=read(fileIn, "CDR0UT")
CDR0BT=read(fileIn, "CDR0BT")
CLR0UT=read(fileIn, "CLR0UT")
close(fileIn)
CLR0UTAB=CLR0UT
AOAUTAB=collect(-180.0:2:180.0)
CLR0BTAB=reshape(CLR0BT,33,11)'
AOABTAB=collect(-32.0:2:32.0)
MACHTAB=collect(0.0:0.1:1.0)
CDR0UTAB=CDR0UT
CDR0BTAB=reshape(CDR0BT,33,11)'
ACL1 = 11.0
ACL2 = 172.0
ACL3 = -5.0
ACL4 = -172.0
DCDMR = 0.002
# wake curvature
KR = -1.0
# schedule from 1.5 in hover to 20, then 0 40 above
KRsched=[1.5, 1.5, 0.0, 0.0]
KRVsched=[0.0, 20.0, 40.0, 160.0]
# main rotor location (station, water line, butt line) [in]
FSMR=341.25
WLMR=315.0
BLMR=0.0

# ----------------------------- MASS PROPERTIES --------------------------------

# CG location (station, water line, butt line) [in]
FSCG=358.0
WLCG=248.2
BLCG=0.0
# aircraft gross weight (rotor blades not included) [lb]
WEIGHT=16000.0
WEIGHTNR = WEIGHT-NB*WBLADE;
MASS=WEIGHTNR/G;
# aircraft moment of inertia [lb-ft²]
IX=4659.0
IY=38512.0
IZ=36796.0
IXZ=1882.0

# --------------------------------- FUSELAGE -----------------------------------

# aerodynamic properties
fileIn = matopen("H60fuselage.mat")
# angle of attack
FUSEAOA=read(fileIn, "FUSEAOA")
# sideslip angle
FUSEBETA=read(fileIn, "FUSEBETA")
# sideslip angle absolute value
FUSEABETA=read(fileIn, "FUSEABETA")
# fuselage drag due to angle of attack [ft²]
FUSEDA=read(fileIn, "FUSEDA")
# fuselage drag due to sideslip [ft²]
FUSEDB=read(fileIn, "FUSEDB")
# fuselage siude force due to sideslip [ft²]
FUSEYB=read(fileIn, "FUSEYB")
# fuselage lift due to angle of attack [ft²]
FUSELA=read(fileIn, "FUSELA")
# fuselage lift due to sideslip [ft²]
FUSELB=read(fileIn, "FUSELB")
# fuselage roll moment due to sideslip [ft³]
FUSERB=read(fileIn, "FUSERB")
# fuselage pitch moment due to angle of attack [ft³]
FUSEMA=read(fileIn, "FUSEMA")
# fuselage pitch moment due to sideslip [ft³]
FUSEMB=read(fileIn, "FUSEMB")
# fuselage yaw moment due to sideslip [ft³]
FUSENB=read(fileIn, "FUSENB")
close(fileIn)
# location of fuselage aerodynamic center
FSfus=345.5
BLfus=0.0
WLfus=234.0
# fuselage CG location (without rotor) (station, water line, butt line) [in]
FSCGB=((WEIGHT*FSCG)-(NB*WBLADE*FSMR))/WEIGHTNR
WLCGB=((WEIGHT*WLCG)-(NB*WBLADE*WLMR))/WEIGHTNR
BLCGB=((WEIGHT*BLCG)-(NB*WBLADE*BLMR))/WEIGHTNR

# -------------------------------- TAIL ROTOR ----------------------------------

# location (station, water line, butt line) [in]
FSTR=732.0
BLTR=-14.0
WLTR=324.7
# cant angle
CantTR=20.0
# nominal RPS
OmegaTR=124.62
# chord
CHRDTR=0.81
# radius
RTR=5.5
# tip loss factor
BTLTR=0.92
# lift slope [1/rad]
a0TR=5.73
# moment of inertia about flapping hinge
IbTR=3.10
# twist
twistTR=-17.2*pi/180
# pitch bias
BIASTR = 5.0
# drag parameters
D0TR =  0.0087
D1TR = -0.0216
D2TR = 0.4
# TR airframe drag coefficient
CDTR=0.0
# tan(δ₃)
TD3TR=0.7002075
# change in coning with thrust
DELTTR = 0.001455
# blockage parameters with vertical tail
BVTTR  = 0.9
BVTTR1 = 1.0
VBVTTR = 30.0
# tail rotor gearing
GEARTR=OmegaTR/OMEGA

# ------------------- HORIZONTAL AND VERTICAL STABILIZERS ----------------------

# HORIZONTAL (idx=1) AND VERTICAL (idx=2) TAIL PROPERTIES
# tail area [ft²]
STAIL=zeros(2)
STAIL[1]=45.0
STAIL[2]=32.3
# cant of the tail (0 for horizontal, -90 for vertical)
PHITAIL=zeros(2)
PHITAIL[1]=0.0
PHITAIL[2]=-90.0
# pitch incidence of tail
ITAIL=zeros(2)
ITAIL[1]=0.0
ITAIL[2]=0.0
# location of tail aero center (station, water line, butt line) [in]
FSTAIL=zeros(2)
FSTAIL[1]=700.0
FSTAIL[2]=695.0
BLTAIL=zeros(2)
BLTAIL[1]=0.0
BLTAIL[2]=0.0
WLTAIL=zeros(2)
WLTAIL[1]=244.0
WLTAIL[2]=273.0
# aerodynamic data
ALTAIL=zeros(2,25)
ALTAIL[1,:]=[collect(-90.0:10:-40.0);collect(-30.0:5:30.0);collect(40.0:10:90.0)]
ALTAIL[2,:]=[collect(-90.0:10:-40.0);collect(-30.0:5:30.0);collect(40.0:10:90.0)]
CLTAIL=zeros(2,25)
CLTAIL[1,:]=[0, -0.294, -0.558, -0.745, -0.847, -0.847, -0.970, -1.030, -1.030,
             -0.930, -0.710, -0.356, 0, 0.356, 0.71, 0.93, 1.03, 1.03, 0.97,
              0.847, 0.847, 0.745, 0.558, 0.294, 0]
CLTAIL[2,:]=[0, -0.12, -0.28, -0.46, -0.66, -0.88, -1.0, -1.0, -0.93, -0.73,
             -0.5, -0.28, -0.06, 0.16, 0.38, 0.61, 0.82, 0.89, 0.89, 0.80, 0.63,
             0.48, 0.32, 0.17, 0.0]
CDTAIL=zeros(2,25)
CDTAIL[1,:]=[1.2000, 1.1610, 1.0500, 0.8880, 0.7020, 0.5310, 0.4300, 0.3700,
             0.3600, 0.1900, 0.0400, 0.0220, 0.01, 0.022, 0.04, 0.19, 0.36,
             0.37, 0.43, 0.531, 0.702, 0.888, 1.05, 1.161, 1.2]
CDTAIL[2,:]=[1.1, 1.025, 0.965, 0.875, 0.745, 0.575, 0.36, 0.265, 0.174, 0.118,
             0.066, 0.033, 0.018, 0.021, 0.044, 0.092, 0.162, 0.248, 0.355,
             0.58, 0.75, 0.875, 0.965, 1.02, 1.08]
# used for elevator control surface (not used here with all moving stab)
CLELEV=zeros(2)
CLELEV[1]=0.0
CLELEV[2]=0.0
# set at non-zero value to fix Horizontal stabilator at fixed position
# with value of 0 will use H60 schedule
STABSET = 0

# ------------------------------ CONTROL MIXING --------------------------------

# conversion from stick percentages to stick inches
perc2in = [0.1 0   0   0
           0   0.1 0   0
           0   0   0.1 0
           0   0   0   5.38/100]
# mixer gain matrix
LNKGAIN = [0.2062 0      0      0
           0      0.2172 0      0
           0      0      0.2025 0
           0      0      0      0.3780]
# convert stick inputs to servo commands
MIXGAIN = [0.6800  0      0.9930   0
           0       1.2020 0.9050  -0.4288
           0      -1.2020 1.2990   0.4288
           0       0      -0.8554  1.6043]
# gain for conversion servo to blade deflections
SWASHGAIN = [11.3413 -5.6706 -5.6706
              0      -5.6706  5.6706
              0       3.5034  3.5034]
TRGAIN = -9.142
# bias for conversion servo to blade deflections
SWASHBIAS = [-7.9738, 9.7280, 10.1111]
TRBIAS = 21.98

# ------------------------- AERODYNAMIC INTERACTION ----------------------------

# EK tables from GenHel
a1f_tab=[-6, 0, 6]
chi_tab=collect(0.0:10.0:100.0)
fileIn = matopen("AeroInteraction.mat")
FKXH1T=read(fileIn, "FKXH1T")
FKXWFT=read(fileIn, "FKXWFT")
FKZH1T=read(fileIn, "FKZH1T")
FKZWFT=read(fileIn, "FKZWFT")
close(fileIn)
ekxf_tab=reshape(FKXWFT,11,3)'
ekzf_tab=reshape(FKZWFT,11,3)'
ekxt_tab=reshape(FKXH1T,11,3)'
ekzt_tab=reshape(FKZH1T,11,3)'
alpha_tab=[-30, -25, -20, 10, 15, 20, 25, 30]*1.0
qlossht_tab=[1.0, 0.875, 0.76, 0.76, 0.82, 0.91, 1.0, 1.0]
psivt_tab=collect(-30.0:5:30.0)
qlossvt_tab=[1.0, 0.88, 0.79, 0.72, 0.66, 0.64, 0.62, 0.64, 0.66, 0.72, 0.79,
             0.88, 1.0]
alphaeps_tab=[collect(-90.0:10:-40.0);collect(-30.0:5:30.0);collect(40.0:10:90.0)]
eps_tab=[0, 0.25, 0.70, 1.2, 1.6, 1.9, 1.8, 1.4, 1.1, 0.8, 0.55, 0.5, 0.45, 0.4,
         0.38, 0.33, 0.19, -0.12, -0.4, -0.7, -0.75, -0.65, -0.45, -0.15, 0]
psisig_tab=[collect(-90.0:10:-40.0);collect(-30.0:5:30.0);collect(40.0:10:90.0)]
sig_tab=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.9, 1.8, 1.5, 0.88, 0.0,
         -0.88, -1.5, -1.8, -1.9, -2.0, -2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# --------------------- LINEARIZATION AND TRIM PARAMETERS ----------------------

# number of fuselage states
NFSTATES = 12
# number of rotor states
NRSTATES = 24
# number of tail rotor states
NTRSTATES = 1
# total number of states
NSTATES = NFSTATES+NRSTATES+NTRSTATES
# number of air frame states
NAFSTATES = NFSTATES+NRSTATES+NTRSTATES
# number of controls
NCTRLS = 4
# fuselage states index vector
IDXF=collect(1:NFSTATES)
# rotor states index vector
IDXR=collect(NFSTATES+1:NFSTATES+NRSTATES)
# tail rotor states index vector
IDXTR=collect(NFSTATES+NRSTATES+1:NFSTATES+NRSTATES+NTRSTATES)
# state scale factors
XSCALE=[1.0; 1.0; 1.0; (pi/180)*ones(6); 1.0; 1.0; 1.0; (pi/180)*ones(16);
        0.2/(OMEGA*R)*ones(3); pi/180; 1.0; 0.01; 0.01; 0.01;
        0.2/(OmegaTR*RTR)]
# perturbations for trim and linearization
DELCLIN=[1.0, 1.0, 1.0, 1.0]
DELXLIN=0.01*XSCALE
# define Trim Variables
TRIMVARS=[collect(1:8);IDXR[[collect(1:19);collect(21:24)]];IDXTR;
          collect(NSTATES+1:NSTATES+NCTRLS)]
TRIMTARG=[collect(1:12);IDXR[[collect(1:19);collect(21:24)]];IDXTR]
# number of azimuth locations to average trim over
NAZTRIM=10
# index of rotor azimuth (needed for azimuthal averaging in trimmer)
IDXAZ=32
