function atan2(y,x)
    if x>0
        return atan(y/x)
    elseif (x<0 && y>=0)
        return atan(y/x)+pi
    elseif (x<0 && y<0)
        return atan(y/x)-pi
    elseif (x==0 && y>0)
        return pi/2
    elseif (x==0 && y<0)
        return -pi/2
    else
        return NaN
    end
end
