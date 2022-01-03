# DESCRIPTION
# Simulates open-loop dynamics using a dynamic inverse controller.
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

function simulate(acfun,finp,Tsim,dt,x0,xdot0,u0)
    # set up simulation
    # vector of times
    time=collect(0:dt:Tsim)
    # initialize state time history
    state=zeros(NSTATES,length(time))
    dstate=zeros(NSTATES,length(time))
    # set initial condition
    state[:,1]=x0
    dstate[:,1]=xdot0
    # simulate open-loop system
    for i=1:length(time)-1
        # integrate
        sol=rk4(H60!,finp,time[i],dt,dstate[:,i],state[:,i],u0)
        state[:,i+1]=sol[1]
        dstate[:,i+1]=sol[2]
    end
    #
    return state, time
end
