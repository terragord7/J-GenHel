# DESCRIPTION
# Closed-loop aircraft dynamics.
#
# INPUT
# - dx: state derivative
# - x: current state vector
# - u_curr: current control vector
# - t: time
#
# OUTPUT
# - xdot: system dynamics
#
# ------------------------------------------------------------------------------

function GenHel_DI!(dx,x,u_curr,t)
    # unpack states
    u = x[1]
    v = x[2]
    w = x[3]
    p = x[4]
    q = x[5]
    r = x[6]
    ϕ = x[7]
    # unpack command filter states
    p_cmd = x[NSTATES+1]
    intperr_cmd = x[NSTATES+2]
    q_cmd = x[NSTATES+3]
    intqerr_cmd = x[NSTATES+4]
    r_cmd = x[NSTATES+5]
    intrerr_cmd = x[NSTATES+6]
    δ_lat = x[NSTATES+7]
    δ_lon = x[NSTATES+8]
    δ_ped = x[NSTATES+9]
    # control perturbations
    Δδ_lat = u_curr[1]-δ_lat
    Δδ_lon = u_curr[2]-δ_lon
    Δδ_ped = u_curr[4]-δ_ped
    # feedforward
    pdot_cmd = -1/τ_p*p_cmd + k_lat/τ_p*Δδ_lat
    qdot_cmd = -1/τ_q*q_cmd + k_lon/τ_q*Δδ_lon
    rdot_cmd = -1/τ_r*r_cmd + k_ped/τ_r*Δδ_ped
    # feedback compensation
    νp = pdot_cmd + kp_roll*(p_cmd-p) + ki_roll*(intperr_cmd)
    νq = qdot_cmd + kp_pitch*(q_cmd-q) + ki_pitch*(intqerr_cmd)
    νr = rdot_cmd + kp_yaw*(r_cmd-r) + ki_yaw*(intrerr_cmd)
    # pseudo control
    ν=[νp, νq, νr]
    # total velocity
    V=sqrt(u^2+v^2+w^2)/1.688
    # matrix interpolation
    idx1=Int(min(max(floor(V/20)+1,1),7))
    idx2=idx1+1;
    intscale=min((V-(idx1-1)*20)/20,1.)
    CA=CATAB[:,:,idx1]+intscale*(CATAB[:,:,idx2]-CATAB[:,:,idx1])
    CBinv=CBinvTAB[:,:,idx1]+intscale*(CBinvTAB[:,:,idx2]-CBinvTAB[:,:,idx1])
    # dynamic inversion
    Δu=CBinv*(ν-CA*[p, q, r])
    # control inputs
    u=[Δu[1]+δ_lat; Δu[2]+δ_lon; u_curr[3]; Δu[3]+δ_ped]
    # state derivative
    der = GenHel!(dx[1:NSTATES],x[1:NSTATES],u,t)
    # rigid-body state derivative
    xdot=zeros(NSTATES+NDISTATES,1)
    xdot[1:NSTATES]=der
    # command filters state derivative
    xdot[NSTATES+1] = pdot_cmd
    xdot[NSTATES+2] = p_cmd-p
    # turn compensation
    if V<40
        xdot[NSTATES+3] = qdot_cmd
    elseif V<60
        xdot[NSTATES+3] = qdot_cmd+r*sin(ϕ)*((V-40)/(60-40))
    else
        xdot[NSTATES+3] = qdot_cmd+r*sin(ϕ)
    end
    xdot[NSTATES+4] = q_cmd-q
    # turn coordination
    if V<40
        xdot[NSTATES+5] = rdot_cmd
    elseif V<60
        xdot[NSTATES+5] = rdot_cmd+g/(V*1.688)*ϕ*((V-40)/(60-40))
    else
        xdot[NSTATES+5] = rdot_cmd+g/(V*1.688)*ϕ
    end
    xdot[NSTATES+6] = r_cmd-r
    xdot[NSTATES+7] = (u_curr[1]-δ_lat)/τ_u
    xdot[NSTATES+8] = (u_curr[2]-δ_lon)/τ_u
    xdot[NSTATES+9] = (u_curr[4]-δ_ped)/τ_u
    #
    return xdot
end
