# DESCRIPTION
# Converts pilot stick percentages to swashplate inputs.
#
# ------------------------------------------------------------------------------

function control_mixing(u)
    # convert stick percentages to inches
    uin = perc2in*u
    # input to mixer
    umix = LNKGAIN*uin
    # servo inputs
    servos=MIXGAIN*umix
    # servo inputs to blade deflections (MR cyclic, MR collective, TR collective)
    control = zeros(4)
    control[1:3] = SWASHGAIN*servos[1:3]+SWASHBIAS
    # MR cyclic convection: β_1c = -A1, β_1s = -A2
    control[1:2] = -control[1:2]
    # add tail rotor collective
    control[4] = TRGAIN*servos[4]+TRBIAS
    #
    return control
end
