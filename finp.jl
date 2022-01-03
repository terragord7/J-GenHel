# DESCRIPTION
# Control input time history.
#
# ------------------------------------------------------------------------------

function finp(t)
    # amplitude of perturbation from trim controls
    amp=5
    # controls history
    if t < 1
        u = [0.0; 0.0; 0.0; 0.0]
    elseif t < 3
        u = [amp; 0.0; 0.0; 0.0]
    elseif t < 5
        u = [-amp; 0.0; 0.0; 0.0]
    else
        u = [0.0; 0.0; 0.0; 0.0]
    end
end
