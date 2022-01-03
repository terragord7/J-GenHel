# DESCRIPTION
# Wraps a phase signal between 0 and 360 deg.
#
# ------------------------------------------------------------------------------

function wrapper(input)
    phase=zeros(length(input),1)
    phase[:]=input[:]
    for i=1:length(phase)
        if phase[i]>0
            n=floor(phase[i]/360)
            phase[i]=phase[i]-n*360
        else
            n=ceil(abs(phase[i])/360)
            phase[i]=phase[i]+n*360
        end
    end
    return phase
end
