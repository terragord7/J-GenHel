# DESCRIPTION
# Model order reduction based on residualization. Keep rigid-body states and
# long and lat flapping angles.
#
# ------------------------------------------------------------------------------

function ModRed10(A,B)
    idxss=[collect(1:8); collect(15:16)]
    idxfs=[collect(13:14); collect(17:31); collect(33:23)]
    Ared=A[idxss,idxss]-A[idxss,idxfs]/A[idxfs,idxfs]*A[idxfs,idxss]
    Bred=B[idxss,1:4]-A[idxss,idxfs]/A[idxfs,idxfs]*B[idxfs,1:4]
    return Ared, Bred
end
