# DESCRIPTION
# Trim initialization.
#
# ------------------------------------------------------------------------------

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
