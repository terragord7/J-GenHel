# TO RUN:
# include("main.jl")
#
# AUTHOR
# Dr. Umberto Saetti, Assistant Professor, Department of Aerospace Engineering,
# Auburn university
# saetti@auburn.edu
#
# LAST UPDATED
# 1/3/2022
#
# DESCRIPTION
# Complex helicopter model based on the General Helicopter (GenHel) flight
# dynamics simulation model (Ref. 1). This model is representative of a utility
# helicopter similar to a UH-60. The model contains a 6-DoF rigid-body dynamic
# model of the fuselage, nonlinear aerodynamic lookup tables for the fuselage,
# rotor blades, and empennage, rigid flap and lead-lag rotor blade dynamics,
# a three-state Pitt-Peters inflow model, and a one-state tail rotor model.
# The trim, linearization routines, and overall architecture of the simulation
# follows the teachings of Dr. Joe Horn at Penn State. If using the code,
# please cite Ref. 2.
#
# REFERENCES
# 1) Howlett, J. J., “UH-60A Black Hawk Engineering Simulation Program.
#    Volume 1: Mathematical Model,” Tech. rep. NASA-CR-166309, 1980.
# 2) Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia
#    Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
#    DOI: https://arc.aiaa.org/doi/10.2514/6.2022-2354.
#
# STATES             | NINDICES    | DESCRIPTION
# ______________________________________________________________________________
# u v w              |  1  2  3    | body velocities [ft/s]
# p q r              |  4  5  6    | angula rates [rad]
# phi theta psi      |  7  8  9    | Euler angles [rad]
# x  y  z            | 10 11 12    | position [ft]
# b0 b1c b1s b0d     | 13 14 15 16 | flapping angles [rad]
# db0 db1c db1s db0d | 17 18 19 20 | flapping angles derivatives [rad/s]
# z0 z1c z1s z0d     | 21 22 23 24 | lead-lag angles [rad]
# dz0 dz1c dz1s dz0d | 25 26 27 28 | lead-lag angles derivative [rad/s]
# lam0 lam1c lam1s   | 29 30 31    | inflow angles [rad]
# psi                | 32          | rotor azimuth [rad]
# dynamic twist      | 33          | dynamic twist force? [?]
# lam0T              | 34          | tail rotor inflow [rad]
#
# CONTROL INPUTS     | INDICED     | DESCRIPTION
# ______________________________________________________________________________
# delta_lat          | 1           | lateral stick
# delta_lon          | 2           | longitudinal stick
# delta_col          | 3           | collective stick
# delta_ped          | 4           | pedals
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
#ENV["MPLBACKEND"]="tkagg"
using Plots
using JLD
#using PyPlot
#using GR

@printf("\n                  J-GenHel           \n")
@printf("\n          Author: Dr. Umberto Saetti\n")
@printf("\n-----------------------------------------------\n")

# include functions
include("control_mixing.jl")
include("table_lookup.jl")
include("stabSched.jl")
include("atan2.jl")
include("eqnmot.jl")
include("interpp1.jl")
include("linearize.jl")
include("interpp2.jl")
include("wrapper.jl")
include("trimmer.jl")
include("ModRed.jl")
include("ModRed8.jl")
include("ModRed10.jl")
include("GenHel.jl")
include("GenHel_DI.jl")
include("rk4.jl")
include("finp.jl")
include("simulate.jl")
include("simulate_DI.jl")
# load H60 constants
include("H60_constants.jl")

# ------------------------ TRIM FLIGHT DYNAMICS MODEL --------------------------

# forward, lateral, and velocities in trim [ft/s]
VXTRIM=80*1.688
VYTRIM=0.0
VZTRIM=.0
# turn rate in trim [rad/s]
PSIDTRIM=0.0
# initial heading [rad]
PSITRIM = 0.0*D2R
# initial position [ft]
ALTTRIM = 0.0
XNTRIM = 0.0
YETRIM = 0.0
# initialize trim variables
include("trimInit.jl")
# trim aircraft
x0, u0 = trimmer(GenHel!,x0,u0,xdot_targ,t)
# linearized aircraft dynamics
A, B = linearize(H60!,x0,u0,xdot0,t)

# ----------------------------- SPECTRAL ANALYSIS ------------------------------

# eigenvalues
eigs=eigvals(A)
# reduce to 8-state model
sysRed = ModRed8(A,B)
A8=sysRed[1]
B8=sysRed[2]
eigs8=eigvals(A8)
# reduce to 10-state model
sysRed = ModRed10(A,B)
A10=sysRed[1]
B10=sysRed[2]
eigs10=eigvals(A10)
# store coefficient matrices and eigenvalues
#save("AB_J-GenHel.jld","A",A,"B",B,"A8",A8,"B8",B8)
#save("eigs_J-GenHel.jld","eigs",eigs,"eigs8",eigs8)

# plot eigenvalues
gr(size=(800,600))
plt_eigs=plot(real(eigs),imag(eigs),seriestype=:scatter,
      label="Full-Order", legendfontsize=12, legend=:topleft,
      marker = (:circle, :royalblue, 6), markerstrokecolor = :royalblue)
plot!(real(eigs10),imag(eigs10),seriestype=:scatter,
      label="Reduced-Order (10-state)", legendfontsize=12,
      marker = (:utriangle, :brown3, 6), markerstrokecolor = :brown3,
      reuse = false, color=:purple)
plot!(real(eigs8),imag(eigs8),seriestype=:scatter,
      label="Reduced-Order (8-state)", legendfontsize=12,
      marker = (:star5, :forestgreen, 6), markerstrokecolor = :forestgreen)
display(plt_eigs)
xaxis!("Real",xguidefontsize=14)
yaxis!("Imag",yguidefontsize=14)

# ---------------------------- FREQUENCY RESPONSES -----------------------------

# transfer functions (full-order model)
sys=ss(A,B,Matrix{Float64}(I, 37, 37),zeros(37,4))
pdlat=tf(sys[4,1])
# transfer functions (10-state model)
sys10=ss(A10,B10,Matrix{Float64}(I, 10, 10),zeros(10,4))
pdlat10=tf(sys10[4,1])
# transfer functions (8-state model)
sys8=ss(A8,B8,Matrix{Float64}(I, 8, 8),zeros(8,4))
pdlat8=tf(sys8[4,1])
# freq response (pdlat)
mag, phase, w = bode(pdlat, 10.0.^(range(-1,stop=2,length=1000)))
mag10, phase10, w10 = bode(pdlat10, 10.0.^(range(-1,stop=2,length=1000)))
mag8, phase8, w8 = bode(pdlat8, 10.0.^(range(-1,stop=2,length=1000)))

# plot (pdlat)
gr(size=(800,600))
p1_mag=plot(vec(w), 20*log10.(vec(mag)/perc2in[1,1]),
      label="Full Order", reuse = false, line=(3, :royalblue, :solid),
      legendfontsize=12)
plot!(vec(w10), 20*log10.(vec(mag10)/perc2in[1,1]), label="Reduced Order (10-state)",
                          reuse = false, line=(3, :brown3, :dash))
plot!(vec(w8), 20*log10.(vec(mag8)/perc2in[1,1]), label="Reduced Order (8-state)",
                          reuse = false, line=(3, :forestgreen, :dashdot))
xaxis!(:log10)
yaxis!("Magnitude [dB]",yguidefontsize=14)
xlims!((0.1,30))
ylims!((-50,0))
p1_phase=plot(vec(w), wrapper(vec(phase)).-360, label="", line=(3, :royalblue, :solid))
plot!(vec(w10), wrapper(vec(phase10)).-360, label="", line=(3, :brown3, :dash))
plot!(vec(w8), wrapper(vec(phase8)).-360, label="", line=(3, :forestgreen, :dashdot))
yaxis!("Phase [deg]",yguidefontsize=14)
xaxis!(:log10)
xlims!((0.1,30))
p1=plot(p1_mag, p1_phase, layout=(2,1), leftmargin=3Plots.mm)
display(p1)

# ------------------------------ TIME SIMULATION -------------------------------

# time step [s]
dt=0.01
# length of simulation [s]
Tsim=10
# run open-loop simulation
state_OL, time_OL = simulate(GenHel!,finp,Tsim,dt,x0,xdot0,u0)
# number of dynamic inverse states
NDISTATES=9
# initial state vector of closed-loop simulation
x0_DI=[x0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; u0[1]; u0[2]; u0[4]]
# initial state derivative vector of closed-loop simulation
xdot0_DI=[xdot0; zeros(NDISTATES,1)]
# run closed-loop simulation
state_CL, time_CL = simulate_DI(GenHel_DI!,finp,Tsim,dt,x0,xdot0,u0)
# store open- and closed-loop simulations
save("resp_GenHel.jld","time_OL",time_OL,"time_CL",time_CL,
    "state_OL",state_OL,"state_CL",state_CL)

# plot attitude
gr(size=(800,600))
i=7
p1=plot(time_OL,state_OL[i,:]*180/pi,label="Open Loop",legend=:topright,
    legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="Closed Loop",legend=:topright,
    legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("φ [deg]",yguidefontsize=14)
xlims!((0,10))
i=8
p2=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
yaxis!("θ [deg]",yguidefontsize=14)
xlims!((0,10))
i=9
p3=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("ψ [deg]",yguidefontsize=14)
xlims!((0,10))
p4=plot(p1, p2, p3, layout=(3,1), leftmargin=3Plots.mm)
display(p4)
#savefig(".\\Plots\\att_J-GenHel.svg")

# plot angular rates
gr(size=(800,600))
i=4
p5=plot(time_OL,state_OL[i,:],label="Open Loop",legend=:bottomleft,
    legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="Closed Loop",legend=:bottomleft,
    legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("p [rad/s]",yguidefontsize=14)
xlims!((0,10))
i=5
p6=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
yaxis!("q [rad/s]",yguidefontsize=14)
xlims!((0,10))
i=6
p7=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
xaxis!("Time [s]",xguidefontsize=14)
yaxis!("r [rad/s]",yguidefontsize=14)
xlims!((0,10))
p8=plot(p5, p6, p7, layout=(3,1), leftmargin=3Plots.mm)
display(p8)
#savefig(".\\Plots\\ang_J-GenHel.svg")
