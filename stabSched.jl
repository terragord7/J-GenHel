# DESCRIPTION
# Horizontal stabilizer scheduling with airspeed and pitch rate.
#
# ------------------------------------------------------------------------------

function stabSched(Coll,Veq,ayg,q,StabSet)
    if (StabSet==0)
        aylim=max(min(ayg,0.13),-0.13)
        stab_ay=-0.0433*min((Veq-30.0),30.0)*aylim
        stab_q=0.16*180/pi*q
        VSchedule=-0.1119*min(max(Veq-80.0,0.0),67.0)
        PR50=min(max(Veq-30.0,0.0),50.0)
        C1=min(max(70.0-Coll,0.),20.0)
        C2=min(max(50.0-Coll,0.),50.0)
        CollSchedule=-1.0*PR50*(0.003286*(C1+C2)+0.67)
        IHT0=42.0+CollSchedule+VSchedule+stab_q+stab_ay
        # stabilator incidence [deg]
        StabInc=min(max(IHT0,-8.0),39.0)
    else
        StabInc=StabSet
    end
end
