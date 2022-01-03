# DESCRIPTION
# Computes the linearized dynamics of a nonlinear system (acfun), given the
# trim state and control vectors (x0 and u0).
#
# INPUT
# - acfun: name of system dynamics function
# - x0: trim state vector
# - u0: trim control vector
# - xdot_in: state derivative
# - t: time
#
# OUTPUT
# - F, G, M, A, B: coefficient matrices of linearized dynamics 
#
# ------------------------------------------------------------------------------

function linearize(acfun,x0,u0,xdot_in,t)
    # init state-space matrices
    F=zeros(NSTATES,NSTATES)
    G=zeros(NSTATES,NCTRLS)
    M=Matrix{Float64}(I, NSTATES, NSTATES)
    A=zeros(NSTATES,NSTATES)
    B=zeros(NSTATES,NCTRLS)
    # F matrix
    x_p=zeros(NSTATES)
    xd_p=zeros(NSTATES)
    u_p=zeros(NCTRLS)
    for k=1:NSTATES
        #x_p=x0
        x_p[:]=x0[:]
        x_p[k]=x_p[k]+DELXLIN[k]
        xdot_p1=acfun(xdot_in,x_p,u0,t)
        #xdot_p1=H60!(xdot_in,x_p,u0,t)
        x_p[k]=x_p[k]-2*DELXLIN[k]
        #xdot_p2=eval(acfun,x_p,u0,xdot_in)
        xdot_p2=acfun(xdot_in,x_p,u0,t)
        F[:,k]=(xdot_p1-xdot_p2)/(2*DELXLIN[k])
    end
    # M matrix
    xd_p=zeros(NSTATES)
    for k=1:NSTATES
        #xd_p=xdot_in
        xd_p[:]=xdot_in[:]
        xd_p[k]=xd_p[k]+DELXLIN[k]
        xdot_p1=acfun(xd_p,x0,u0,t)
        xd_p[k]=xd_p[k]-2*DELXLIN[k]
        xdot_p2=acfun(xd_p,x0,u0,t)
        M[:,k]=M[:,k]-(xdot_p1-xdot_p2)/(2*DELXLIN[k])
    end
    # G matrix
    u_p=zeros(NCTRLS)
    for k=1:NCTRLS
        #u_p=u0
        u_p[:]=u0[:]
        u_p[k]=u_p[k]+DELCLIN[k]
        xdot_p1=acfun(xdot_in,x0,u_p,t)
        u_p[k]=u_p[k]-2*DELCLIN[k]
        xdot_p2=acfun(xdot_in,x0,u_p,t)
        G[:,k]=(xdot_p1-xdot_p2)/(2*DELCLIN[k])
    end
    A=M\F
    B=M\G
    return F, G, M, A, B
end
