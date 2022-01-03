function table_lookup(datin,datout,in)
    out=datout[1]
    if in<=datin[1]
        done=1
    else
        done=0
    end
    k=1
    while done==0
        k=k+1
        if in<=datin[k]
            out=datout[k-1]+(datout[k]-datout[k-1])/(datin[k]-datin[k-1])*(in-datin[k-1])
            done=1
        elseif k>=length(datin)
            out=datout[k]
            done=1
        end
    end
    return out
end
