# DESCRIPTION
# Rigid-body equations of motion.
#
# INPUT
# - xf: fuselage states
# - Ftot: net forces acting at the CG
# - Mtot: net moments acting on the CG
#
# OUTPUT
# - xfdot: rigid-body dynamics 
#
# ------------------------------------------------------------------------------

function eqnmot(xf,Ftot,Mtot)
    # unpack forces and moments
    X=Ftot[1]
    Y=Ftot[2]
    Z=Ftot[3]
    L=Mtot[1]
    M=Mtot[2]
    N=Mtot[3]
    # unpack states
    u=xf[1]
    v=xf[2]
    w=xf[3]
    p=xf[4]
    q=xf[5]
    r=xf[6]
    phi=xf[7]
    theta=xf[8]
    psi=xf[9]
    # trig functions
    cphi=cos(phi)
    sphi=sin(phi)
    cthe=cos(theta)
    sthe=sin(theta)
    cpsi=cos(psi)
    spsi=sin(psi)
    # state derivative vector
    xf_dot=zeros(12)
    # equations of motion
    # velocities
    xfdot=zeros(NFSTATES)
    xfdot[1]=X/MASS-G*sthe-q*w+r*v
    xfdot[2]=Y/MASS+G*cthe*sphi-r*u+p*w
    xfdot[3]=Z/MASS+G*cthe*cphi-p*v+q*u
    # angular rates
    gam=IX*IZ-IXZ^2;
    xfdot[4]=(IZ*L+IXZ*N+IXZ*(IX-IY+IZ)*p*q-(IZ^2-IY*IZ+IXZ^2)*q*r)/gam
    xfdot[5]=(M+(IZ-IX)*p*r-IXZ*(p^2-r^2))/IY
    xfdot[6]=(IX*N+IXZ*L-IXZ*(IX-IY+IZ)*q*r+(IX^2-IX*IY+IXZ^2)*p*q)/gam
    # attitude
    xfdot[7]=p+q*sphi*sthe/cthe+r*cphi*sthe/cthe
    xfdot[8]=q*cphi-r*sphi
    xfdot[9]=q*sphi/cthe+r*cphi/cthe
    # position
    xfdot[10]=u*cthe*cpsi+v*(sphi*sthe*cpsi-cphi*spsi)+w*(cphi*sthe*cpsi+sphi*spsi)
    xfdot[11]=u*cthe*spsi+v*(sphi*sthe*spsi+cphi*cpsi)+w*(cphi*sthe*spsi-sphi*cpsi)
    xfdot[12]=-u*sthe+v*sphi*cthe+w*cphi*cthe
    # return state derivative vector
    return xfdot
end
