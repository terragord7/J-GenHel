# DESCRIPTION
# Runge-Kutta 4 integration scheme.
#
# INPUT
# - fdyn: name of system dynamics function
# - finp: name of control input function
# - t: time
# - dt: time step
# - dx_curr: current state derivative vector
# - x_curr: current state vector
# - u0: initial control vector
#
# OUTPUT
# - x_new: new state vector after integration
# - dx_new: new state derivative vector after integration
#
# ------------------------------------------------------------------------------

function rk4(fdyn,finp,t,dt,dx_curr,x_curr,u0)
    #
    u=finp(t)+u0
    xd=fdyn(dx_curr,x_curr,u,t)
    dx_new=xd
    k1=dt*xd
    #
    u=finp(t+0.5*dt)+u0
    xd=fdyn(dx_curr,x_curr+0.5*k1,u,t+0.5*dt)
    k2=dt*xd
    #
    xd=fdyn(dx_curr,x_curr+0.5*k2,u,t+0.5*dt)
    k3=dt*xd
    #
    u=finp(t+dt)+u0
    xd=fdyn(dx_curr,x_curr+k3,u,t+dt)
    k4=dt*xd
    # weighted average
    x_new=x_curr+k1/6+k2/3+k3/3+k4/6
    #
    return x_new, dx_new
end
