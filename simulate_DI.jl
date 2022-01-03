# DESCRIPTION
# Simulates closed-loop dynamics using a dynamic inverse controller.
#
# INPUT
# - acfun: name of system dynamics function
# - finp: name of control input function
# - Tsim: total simulation time
# - dt: time step
# - x0: initial state vector
# - xdot0: initial state derivative vector
# - u0: initial control vector
#
# OUTPUT
# - state: state vector history
# - time: time vector history
#
# ------------------------------------------------------------------------------

function simulate_DI(acfun,finp,Tsim,dt,x0,xdot0,u0)
    # set up simulation
    # load DI parameters
    global k_lat = load("DI.jld","k_lat")
    global k_lon = load("DI.jld","k_lon")
    global k_ped = load("DI.jld","k_ped")
    global kp_roll = load("DI.jld","kp_roll")
    global ki_roll = load("DI.jld","ki_roll")
    global kp_pitch = load("DI.jld","kp_pitch")
    global ki_pitch = load("DI.jld","ki_pitch")
    global kp_yaw = load("DI.jld","kp_yaw")
    global ki_yaw = load("DI.jld","ki_yaw")
    global τ_p = load("DI.jld","tau_p")
    global τ_q = load("DI.jld","tau_q")
    global τ_r = load("DI.jld","tau_r")
    global τ_u = load("DI.jld","tau_u")
    global CATAB = load("DI.jld","CATAB")
    global CBinvTAB = load("DI.jld","CBinvTAB")
    # vector of times
    time=collect(0:dt:Tsim)
    # initialize state time history
    state=zeros(NSTATES+NDISTATES,length(time))
    dstate=zeros(NSTATES+NDISTATES,length(time))
    # set initial condition
    state[:,1]=x0_DI
    dstate[:,1]=xdot0_DI
    # simulate closed-loop system
    for i=1:length(time)-1
        # integrate
        sol=rk4(acfun,finp,time[i],dt,dstate[:,i],state[:,i],u0)
        state[:,i+1]=sol[1]
        dstate[:,i+1]=sol[2]
    end
    #
    return state, time
end
