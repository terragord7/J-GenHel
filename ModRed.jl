# DESCRIPTION
# Model order reduction based on residualization. 
#
# INPUT
# - A: system matrix
# - B: control matrix
# - idxss: indices for slow states
# - idxfs: indices fro fast states
#
# OUTPUT
# - Ared: reduced-order system matrix
# - Bred: reduced-order control matrix

function ModRed(A,B,idxss,idxfs)
    Ared=A[idxss,idxss]-A[idxss,idxfs]/A[idxfs,idxfs]*A[idxfs,idxss]
    Bred=B[idxss,1:4]-A[idxss,idxfs]/A[idxfs,idxfs]*B[idxfs,1:4]
    return Ared, Bred
end
