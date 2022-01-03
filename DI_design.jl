# To run:
# include("DI_design.jl")

# DESCRIPTION
# Trims and linearizes the rotorcraft model at incremental flight speeds, and
# generates the gains of the dynamic inversion controller. Always run after
# modifying the aircraft parameters in the AIRCRAFTNAME_constants.jl file.
#
# ------------------------------------------------------------------------------

# include packages
#using DifferentialEquations
using LinearAlgebra
using SparseArrays
using ControlSystems
using MAT
using Interpolations
using ApproxFun
using Printf
using Plots
using JLD

# include functions
include("H60_mixing.jl")
include("table_lookup.jl")
include("H60_StabSched.jl")
include("atan2.jl")
include("eqnmot.jl")
include("interpp1.jl")
include("linearize.jl")
include("interpp2.jl")
include("wrapper.jl")
include("H60_trimmer.jl")
include("ModRed.jl")
include("ModRed8.jl")
include("H60.jl")
#include("H60_DI.jl")
include("rk4.jl")
include("finp.jl")
include("simulate.jl")
# load H60 constants
include("H60_constants.jl")
# vector of trim speeds [kts]
VXTRIM_vec=[0.01,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0]
# initialize DI matrices
CATAB = zeros(3,3,length(VXTRIM_vec))
CBinvTAB = zeros(3,3,length(VXTRIM_vec))
x0_mat = zeros(NSTATES,length(VXTRIM_vec))
u0_mat = zeros(NCTRLS,length(VXTRIM_vec))
# trim at different speeds
for iv=1:length(VXTRIM_vec)
    # load H60 constants
    include("H60_constants.jl")
    # forward, lateral, and velocities in trim [ft/s]
    VXTRIM=VXTRIM_vec[iv]*1.688
    VYTRIM=0.0
    VZTRIM=0.0
    # turn rate in trim [rad/s]
    PSIDTRIM=0.0
    # initial heading [rad]
    PSITRIM = 0.0*D2R
    # initial position [ft]
    ALTTRIM = 0.0
    XNTRIM = 0.0
    YETRIM = 0.0
    # initialize state vectors and initial guess
    # xf = [u v w p q r phi theta psi x y z]'
    xf0=zeros(NFSTATES)
    # xr = [MBC_flap, dMBC_flap, MBC_lag, dMBC_lag, MBC_inflow, azimuth, dyn twist]'
    xr0=zeros(NRSTATES)
    # xtr = TR_inflow
    xtr0=zeros(NTRSTATES)
    # initial guess for controls [%]
    u0=[50.0; 50.0; 50.0; 50.0]
    # set non-zero inflow states
    xr0[17]=0.05
    # initialize dynamic twist blade loads value [lb]
    xr0[21]=WEIGHT/NB
    # tail rotor inflow
    xtr0=0.05
    # initial guess for trim solution
    xf0[1]=VXTRIM
    xf0[2]=VYTRIM
    xf0[3]=VZTRIM
    xf0[9]=PSITRIM
    xf0[10]=XNTRIM
    xf0[11]=YETRIM
    xf0[12]=-ALTTRIM
    x0=[xf0; xr0; xtr0]
    # set trim targets
    xfdot0=zeros(NFSTATES)
    xfdot0[9]=PSIDTRIM
    xfdot0[10:12]=[VXTRIM*cos(PSITRIM)-VYTRIM*cos(PSITRIM); VXTRIM*sin(PSITRIM)+
                   VYTRIM*cos(PSITRIM); VZTRIM]
    xrdot0=zeros(NRSTATES)
    xtrdot0=zeros(NTRSTATES)
    xdot0=[xfdot0; xrdot0; xtrdot0]
    xdot_targ=xdot0[TRIMTARG]
    # Option to trim in zero sideselip or zero bank angle
    if (VXTRIM>=60*1.688)
        if PSIDTRIM==0
            # Trim heading for zero bank angle above 60 kts
            x0[7]=0.0;
            TRIMVARS[7:8]=[8; 9]
        else
            x0[9]=atan2(VYTRIM,VXTRIM)
            TRIMVARS[7:8]=[7; 8]
            # initial guess for phi (using turn coord)
            x0[7]=atan(PSIDTRIM*sqrt(VXTRIM^2+VYTRIM^2+VZTRIM^2)/G)
            # initial guess for r
            x0[6]=G/sqrt(VXTRIM^2+VYTRIM^2+VZTRIM^2)*sin(x0[7])
            # initial guess for q
            x0[5]=x0[6]*tan(x0[7])
        end
    else
        # otherwise set heading to flight path
        x0[9]=atan2(VYTRIM,VXTRIM)
        TRIMVARS[7:8]=[7; 8]
    end
    # time [s]
    t=0.0
    # trim aircraft
    x0, u0 = H60_trimmer(H60!,x0,u0,xdot_targ,t)
    # store trim state and controls
    x0_mat[:,iv] = x0
    u0_mat[:,iv] = u0
    # linearized aircraft dynamics
    A, B = linearize(H60!,x0,u0,xdot0,t)
    # reduce to 8-state model
    sysRed = ModRed8(A,B)
    A8=sysRed[1]
    B8=sysRed[2]
    eigs8=eigvals(A8)
    # dynamic inverse matrices
    C_DI=[1 0 0
          0 1 0
          0 0 1]
    CATAB[:,:,iv]=C_DI*A8[4:6,4:6]
    CBinvTAB[:,:,iv]=inv(C_DI*B8[4:6,[1,2,4]])
end
# command filters time constants
tau_p=1/3.5
tau_q=1/4.5
tau_r=0.5
# washout filter for controls
tau_u=5.0;
# disturbance rejection frequency and damping and integrator poles
wnroll_d=3.5
dmproll_d=1.0
wnpitch_d=4.5
dmppitch_d=1.0
wnyaw_d=2.;
dmpyaw_d=1.;
# feedback gains
kp_roll=2*wnroll_d*dmproll_d
ki_roll=wnroll_d^2
kp_pitch=2*wnpitch_d*dmppitch_d
ki_pitch=wnpitch_d^2
kp_yaw=2*dmpyaw_d*wnyaw_d
ki_yaw=wnyaw_d^2;
# pilot gains
k_lat=pi/100
k_lon=pi/100
k_ped=0.644
# save DI matrices and gains into a file
save("DI.jld","CATAB",CATAB,"CBinvTAB",CBinvTAB,"kp_roll",kp_roll,"ki_roll",
    ki_roll,"kp_pitch",kp_pitch,"ki_pitch",ki_pitch,"kp_yaw",kp_yaw,"ki_yaw",
    ki_yaw,"k_lat",k_lat,"k_lon",k_lon,"k_ped",k_ped,"tau_p",tau_p,
    "tau_q",tau_q,"tau_r",tau_r,"tau_u",tau_u)
