# DESCRIPTION
# Trims the aircraft dynamics at an arbitrary flight condition.
#
# INPUT
# - acfun: name of system dynamics function
# - x0: initial guess for trim state vector
# - u0: initial guess for trim control vector
# - targ_des: vector of trim targets
# - t: time
#
# OUTPUT
# - x0trim: trim state vector
# - u0trim: trim control vector
# - itrim: number of trim iterations
#
# ------------------------------------------------------------------------------

function trimmer(acfun,x0,u0,targ_des,t)
    # trim state and controls
    x0trim=x0
    u0trim=u0
    #
    it=0
    err=100.0
    trim_tol=5e-4
    itmax=80
    #
    xdot_in=zeros(NSTATES)
    idx_targ_state=findall(x -> x<=NSTATES, TRIMTARG)
    xdot_in[TRIMTARG[idx_targ_state]]=targ_des[1:length(idx_targ_state)]
    daztrim=2*pi/NB/NAZTRIM
    #
    @printf("trimming the aircraft\n")
    #
    while ((it<itmax) && (err>trim_tol))
        it=it+1
        # zimuthal average of xdot
        xdot0=zeros(NSTATES)
        for iaz=1:NAZTRIM
            x0trim[IDXAZ]=(iaz-1)*daztrim
            xdotaz=acfun(xdot_in,x0trim,u0trim,t)
            xdot0=xdot0+xdotaz/NAZTRIM
        end
        x0trim[IDXAZ]=0.0
        # or: do not average over one rotor revolution
        # xdot0=H60!(xdot_in,x0trim,u0trim,t)
        targvec=xdot0[TRIMTARG]
        targ_err=targvec-targ_des
        err=maximum(abs.(targ_err)./XSCALE[TRIMTARG])
        @printf("%d          %f\n",it,err)
        if (err>trim_tol)
            F, G, M=linearize(acfun,x0trim,u0trim,xdot_in,t)
            # don't use M in trimmer, as it is prone to diverge
            Jac=[F G]
            Jac2 = Jac[TRIMTARG,TRIMVARS]
            trimvec=[x0trim;u0trim]
            trimvec[TRIMVARS]=trimvec[TRIMVARS]-0.5*pinv(Jac2)*targ_err
            x0trim=trimvec[1:NSTATES]
            u0trim=trimvec[NSTATES+1:NSTATES+NCTRLS]
        end
    end
        #
        if (err>trim_tol)
            @printf("warning: trim not achieved\n\n")
            itrim=0
        else
            @printf("successful trim\n\n")
            itrim=1
        end
        #
        return x0trim, u0trim, itrim
end
