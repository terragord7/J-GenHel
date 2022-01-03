# DESCRIPTION
# Model order reduction based on residualization. Keep only rigid-body states.
#
# ------------------------------------------------------------------------------

function ModRed8(A,B)
    idxss=collect(1:8)
    idxfs=[collect(13:31); collect(33:23)]
    Ared=A[idxss,idxss]-A[idxss,idxfs]/A[idxfs,idxfs]*A[idxfs,idxss]
    Bred=B[idxss,1:4]-A[idxss,idxfs]/A[idxfs,idxfs]*B[idxfs,1:4]
    return Ared, Bred
end
