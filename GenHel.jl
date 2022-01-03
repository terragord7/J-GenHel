# DESCRIPTION
# Open-loop aircraft dynamics.
#
# INPUT
# - dx: state derivative
# - x: current state vector
# - p: current control vector
# - t: time
#
# OUTPUT
# - xdot: system dynamics
#
# ------------------------------------------------------------------------------

function GenHel!(dx,x,p,t)

    # control system mixing
    controls = control_mixing(p)

    # PARTITION STATE VECTOR
    # fuselage states
    xf=x[1:NFSTATES]
    # rotor states
    xr=x[NFSTATES+1:NFSTATES+NRSTATES]
    # tail rotor states
    xtr=x[NFSTATES+NRSTATES+1:NFSTATES+NRSTATES+NTRSTATES]

    # ROTOR INTERFERENCE
    # body to shaft coordinate transform
    Tshaft=[cos(ITHSH*D2R)                 0.0            -sin(ITHSH*D2R)
            sin(ITHSH*D2R)*sin(IPHSH*D2R)  cos(IPHSH*D2R)  cos(ITHSH*D2R)*sin(IPHSH*D2R)
            sin(ITHSH*D2R)*cos(IPHSH*D2R) -sin(IPHSH*D2R)  cos(ITHSH*D2R)*cos(IPHSH*D2R)]
    # body-frame velocities
    u=x[1]
    v=x[2]
    w=x[3]
    # longitudinal flap angle
    beta1c=x[15]
    # convert longitudinal flap angle to different sign convention and deg
    A1F=-beta1c*R2D
    # average harmonic inflow
    lambda0=x[29]
    cg=0.0
    Gif=1.0-cg
    temp=Tshaft*[u;v;w]
    us=temp[1]
    vs=temp[2]
    ws=temp[3]
    muxs=us/(OMEGA*R)
    muzs=ws/(OMEGA*R)
    chi=R2D*atan2(muxs,abs(lambda0-muzs))+A1F
    alpha=atan2(w,abs(u))*R2D
    beta=atan2(v,sqrt(u^2+w^2))*R2D
    psif=-beta
    # check these two lines
    chitab=min(max(chi,0.0),100.0)
    A1Ftab=min(max(A1F,-6.0),6.0)
    # interpolate
    itp1=interpolate((chi_tab, a1f_tab), ekxf_tab', Gridded(Linear()))
    ekxf=itp1(chitab,A1Ftab)
    itp2=interpolate((chi_tab, a1f_tab), ekzf_tab', Gridded(Linear()))
    ekzf=itp2(chitab,A1Ftab)
    itp3=interpolate((chi_tab, a1f_tab), ekxt_tab', Gridded(Linear()))
    ekxt=itp3(chitab,A1Ftab)
    itp4=interpolate((chi_tab, a1f_tab), ekzt_tab', Gridded(Linear()))
    ekzt=itp4(chitab,A1Ftab)
    #
    rotor_if=zeros(8)
    rotor_if[1]=-Gif*ekzf*lambda0*OMEGA*R
    rotor_if[2]=Gif*ekxf*lambda0*OMEGA*R
    rotor_if[3]=-Gif*ekzt*lambda0*OMEGA*R
    rotor_if[4]=Gif*ekxt*lambda0*OMEGA*R
    rotor_if[5]=table_lookup(alpha_tab,qlossht_tab,alpha)
    rotor_if[6]=table_lookup(psivt_tab,qlossvt_tab,psif)
    rotor_if[7]=table_lookup(alphaeps_tab,eps_tab,alpha)
    rotor_if[8]=table_lookup(psisig_tab,sig_tab,psif)

    # FUSELAGE
    # external velocities zeroed out in trim and linearization
    VextF=zeros(3)
    VextTR=zeros(3)
    VextHT=zeros(3)
    VextVT=zeros(3)
    # fudelage geometry [ft]
    FWT = (FSCG-FSfus)/12
    WWT = (WLCG-WLfus)/12
    BWT = (BLCG-BLfus)/12
    # gust velocities
    VXGWF = 0.0
    VYGWF = 0.0
    VZGWF = 0.0
    # NED to Body transformation
    sphi=sin(xf[7])
    cphi=cos(xf[7])
    sthe=sin(xf[8])
    cthe=cos(xf[8])
    spsi=sin(xf[9])
    cpsi=cos(xf[9])
    TNED2body=[cthe*cpsi                  cthe*spsi                  -sthe
               (sphi*sthe*cpsi-cphi*spsi) (sphi*sthe*spsi+cphi*cpsi)  sphi*cthe
               (cphi*sthe*cpsi+sphi*spsi) (cphi*sthe*spsi-sphi*cpsi)  cphi*cthe]
    # external velocities in body coordinates
    Vextb=TNED2body*VextF
    # interference velocities
    VXIWF = rotor_if[2]+Vextb[1]
    VYIWF = 0.0+Vextb[2];
    VZIWF = rotor_if[1]+Vextb[3]
    # air velocity inputs
    VXB=xf[1]
    VYB=xf[2]
    VZB=xf[3]
    # fuselage velocities in body axes
    VXWF = VXB+VXGWF+VXIWF
    VYWF = VYB+VYGWF+VYIWF
    VZWF = VZB+VZGWF+VZIWF
    # angle of attack and sideslip
    AVXWF=abs(VXWF);
    if AVXWF < 0.000001
        AVXWF = 0.000001
    end
    ALFWFR = atan(VZWF/AVXWF)
    ALFWF = R2D*ALFWFR
    FVTERM = VXWF^2+VZWF^2
    #
    if FVTERM < 0.000001
        FVTERM = 0.000001
    end
    RVTERM = 1/sqrt(FVTERM)
    CALFWF = VXWF*RVTERM
    SALFWF = VZWF*RVTERM
    #
    BETWFR = atan(VYWF*RVTERM)
    BETWF  = R2D*BETWFR
    SBETWF = sin(BETWFR)
    CBETWF = cos(BETWFR)
    PSIWF = -BETWF
    # dynamic pressure
    QWF = 0.5*RHO*(VYWF^2 + FVTERM)
    # fuselage aerodyanmic loading coefficients
    AALFWF = abs(ALFWF)
    APSIWF = abs(PSIWF)
    SGNPSI = sign(PSIWF)
    # functions of angle of attack
    CDA = table_lookup(FUSEAOA,FUSEDA,ALFWF)
    CLA = table_lookup(FUSEAOA,FUSELA,ALFWF)
    CMA = table_lookup(FUSEAOA,FUSEMA,ALFWF)
    # functions of sideslip
    CLB = table_lookup(FUSEBETA,FUSELB,PSIWF)
    CMB = table_lookup(FUSEABETA,FUSEMB,APSIWF)
    CDB = table_lookup(FUSEABETA,FUSEDB,APSIWF)
    CYB = table_lookup(FUSEBETA,FUSEYB,PSIWF)
    CRB = table_lookup(FUSEBETA,FUSERB,PSIWF)
    CNB = table_lookup(FUSEBETA,FUSENB,PSIWF)
    # total coefficients
    CDTOT = CDA + CDB
    CLTOT = CLA + CLB
    CYTOT = CYB
    CRTOT = CRB
    CMTOT = CMA + CMB
    CNTOT = CNB
    # transform forces and moments from wind axes to body axes
    Tbw=[CALFWF*CBETWF  CALFWF*SBETWF -SALFWF
         SBETWF        -CBETWF         0.0
         SALFWF*CBETWF  SALFWF*SBETWF  CALFWF]
    Ffw=-QWF*[CDTOT;CYTOT;CLTOT]
    Ff=Tbw*Ffw
    # check this line
    Mf=Tbw*QWF*[CRTOT;-CMTOT;CNTOT]+[0.0 -WWT  BWT; WWT 0.0 -FWT; -BWT FWT 0.0]*Ff
    pfusNED=xf[10:12]+TNED2body'*[FWT;BWT;WWT]

    # TAIL ROTOR
    # absolute speed
    VKT = sqrt(u^2+w^2)*fps2kts
    # tail rotor twist [rad]
    TWSTTR = twistTR/D2R
    # ratio between main and tail rotor angular speeds
    OMGRAT = 1 # OMEGA/OMEGAT
    # tail rotor collective
    THETTRC = controls[4]
    # constants
    C12TR = 0.5
    C13TR = 0.333333333
    C14TR = 0.25
    C15TR = 0.2
    C16TR = 0.166666666
    C23TR = 0.666666666
    C25TR = 0.4
    C43TR = 1.333333333
    C54TR = 1.25
    C58TR = 0.625
    C83TR = 2.666666666
    C89TR = 0.888888888
    # dynamic inflow constant
    XKINF = 4.0/(3.0*pi)
    # angular rates
    PB=xf[4]
    QB=xf[5]
    RB=xf[6]
    # external velocities in body frame
    Vextb=TNED2body*VextTR
    # tail rotor geometry [ft]
    XTR = -(FSTR-FSCGB)/12.0
    ZTR = -(WLTR-WLCGB)/12.0
    YTR = (BLTR-BLCGB)/12.0
    # body to TR transform
    cgam=cos((90-CantTR)*D2R)
    sgam=sin((90-CantTR)*D2R)
    Tcant=[1.0  0.0  0.0
           0.0  cgam sgam
           0.0 -sgam cgam]
    #
    GTR = a0TR*NB*CHRDTR/(2.0*pi*RTR)
    #
    S0 = a0TR*CHRDTR*RTR^4/IbTR
    S1 = BTLTR/2.0
    S2 = S1*BTLTR
    S3 = S2^2
    S4 = S2/2
    S5 = (BTLTR^3)/3.0
    S6 = DELTTR*TD3TR;
    CBLK = 1.0/(VBVTTR^2)
    #
    B2 = BTLTR^2
    B3 = BTLTR*B2
    B4 = B2^2
    B5 = B2*B3
    B6 = B3^2
    B7 = B3*B4
    B8 = B4^2
    B9 = B4*B5
    B10 = B5^2
    B11 = B5*B6
    B12 = B6^2
    #
    S10 = B2/1296.0
    S11 = C83TR*BTLTR
    S12 = B9/864.0
    S13 = 2.0*B2
    S14 = B10/1080.0
    S15 = C89TR*B2
    S16 = B10/2304.0
    S17 = C43TR*B3
    S18 = B11/1440.0
    S19 = B4/2.0
    S20 = B12/3600.0
    S21 = -B5/108.0
    S22 = -B6/144.0
    S23 = -B7/180.0
    S24 = B2/36.0
    S25 = -C14TR+(1.0/B2)+(C12TR/B4)
    S26 = (B4/162.0)-(B5/81.0)+(B6/144.0)
    S27 = C43TR/BTLTR*(1.0+(1.0/B2))
    S28 = (B5/108.0)-(B6/54.0)+(B7/96.0)
    S29 = 1.0+(1.0/B2)
    S30 = (B6/135.0)*(1.0-(2.0*BTLTR))+(B8/120.0)
    S31 = C14TR+(C89TR/B2)
    S32 = (B6/288.0)-(B7/144.0)+(B8/256.0)
    S33 = C13TR+(C43TR/BTLTR)
    S34 = (B7/180.0)-(B8/90.0)+(B9/160.0)
    S35 = (B8/450.0)-(B9/225.0)+(B10/400.0)
    #
    BLDSTR = NB
    SIGTR = BLDSTR*CHRDTR/(pi*RTR)
    THT1TR = TWSTTR*D2R
    #
    S7 = pi*((OmegaTR*OMGRAT*RTR^2)^2)
    S8 = 1.0/(OmegaTR*RTR*OMGRAT)
    S9 = S7*RTR
    # dynamic pressure ratio (same as for vertical tail)
    idxtail = 1
    xtail=-(FSTR[idxtail]-FSCG)/12.0
    ytail=(BLTR[idxtail]-BLCG)/12.0
    ztail=-(WLTR[idxtail]-WLCG)/12.0
    qlossfac=sqrt(rotor_if[6])
    ublocal=xf[1]*qlossfac+Vextb[1]+rotor_if[2]+(QB)*ztail-(RB)*ytail
    vblocal=(xf[2]-D2R*xf[1]*rotor_if[8])*qlossfac+Vextb[2]-(PB)*ztail+(RB)*xtail;
    wblocal=(xf[3]-D2R*xf[1]*rotor_if[7])*qlossfac+Vextb[3]+rotor_if[1]+(PB)*ytail-(QB)*xtail
    # ransform velocities to orientation of tail
    temp=Tcant*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]
    # advance ratios
    XMUXTR = ut*S8
    XMUYTR = vt*S8
    XMUZTR = wt*S8
    XMU2 = ((XMUXTR^2) + (XMUYTR^2))
    # Bailey coefficents
    BT31 = S2+(0.25*XMU2)
    BT32 = S5+(S1*XMU2)
    BT33 = S3+(S4*XMU2)
    # additional Bailey coefficients for torque
    XGAM = RHO*S0
    XGAM2 = XGAM^2
    BT41  = S2+(C54TR+(S10*XGAM2))*XMU2
    BT42  = S5+(S11+(S12*XGAM2))*XMU2
    BT43  = S3+(S13+(S14*XGAM2))*XMU2
    BT44  = (S15+(S16*XGAM2))*XMU2
    BT45  = (S17+(S18*XGAM2))*XMU2
    BT46  = (S19+(S20*XGAM2))*XMU2
    BT47  = S21*XMU2
    BT48  = S22*XMU2
    BT49  = S23*XMU2
    BT410 = S24*XMU2
    BT51  = C14TR*(1.0+XMU2)
    BT52  = C13TR
    BT53  = BT51
    BT54  = C15TR+(C16TR*XMU2)
    BT55  = C12TR+((S25+(S26*XGAM2))*XMU2)
    BT56  = C23TR+((S27+(S28*XGAM2))*XMU2)
    BT57  = C12TR+((S29+(S30*XGAM2))*XMU2)
    BT58  = C14TR+((S31+(S32*XGAM2))*XMU2)
    BT59  = C25TR+((S33+(S34*XGAM2))*XMU2)
    BT510 = C16TR+((C58TR+(S35*XGAM2))*XMU2)
    # vertical tail blockage factor
    BLKTR = BVTTR1
    if VKT < VBVTTR
            BLKTR = ((1.0-BVTTR)*CBLK*(VKT^2))+BVTTR
    end
    TTR = 200
    DWTRSS = xtr[1]
    XLAMTR=XMUZTR-xtr[1]
    THETTR = D2R*(THETTRC)
    iter=0
    err=10
    # iterate for steady-state induced velocity and thrust
    while abs(err)>1e-6 && iter<200
          # tail rotor blade pitch [rad]
          THETTR = D2R*(THETTRC-(TTR*S6)+BIASTR)
          # tail rotor downwash and thrust
          DWTRSS = GTR*((XMUZTR*BT31)+(THETTR*BT32)+(D2R*TWSTTR*BT33))/
                   (2.0*sqrt(XMU2+((XMUZTR-DWTRSS)^2))+(GTR*BT31))
          XLAMTR = XMUZTR-DWTRSS
          TTR_new = 2.0*DWTRSS*sqrt(XMU2+(XLAMTR^2))*RHO*S7*BLKTR
          err = TTR_new-TTR
          # this thrust is steady-state version, it will get overwritten for current inflow
          TTR = TTR+0.5*err
          iter=iter+1
    end
    if iter>=199
        @printf("warning: TR not converged")
    end
    # inflow dynamics (no time delay)
    VT = sqrt(XMU2+(XLAMTR^2))
    CDWTR = XKINF/(OmegaTR*OMGRAT*VT)
    # inflow ratio nflow ratio
    XLAMTR = XMUZTR-xtr[1]
    # tail rotor state derivative (check this line: 10^-6 vs 10^-9)
    xtrdot=(DWTRSS-xtr[1])/(CDWTR)
    # tail rotor thrust
    CTHTR = GTR*((XLAMTR*BT31)+(THETTR*BT32)+(D2R*TWSTTR*BT33))
    TTR = CTHTR*RHO*S7*BLKTR
    # tail rotor torque
    ATR = a0TR
    CQTR = (((D2TR*BT55)-(ATR*BT41))*(XLAMTR^2)+((D2TR*BT56)-
           (ATR*BT42))*THETTR*XLAMTR+((D2TR*BT57)-(ATR*BT43))*THT1TR*XLAMTR+
           ((D2TR*BT58)-(ATR*BT44))*(THETTR^2)+((D2TR*BT59)-
           (ATR*BT45))*THT1TR*THETTR+((D2TR*BT510)-(ATR*BT46))*(THT1TR^2)+
           (D0TR*BT51 )+D1TR*((BT52*XLAMTR)+(BT53*THETTR)+
           (BT54*THT1TR)))*0.5*SIGTR
    Qtr = CQTR*RHO*S9
    # tail rotor drag
    DTR=0.5*CDTR*RHO*(ut^2)
    # tail rotor horsepower [hp]
    HPTR = Qtr*OmegaTR*OMGRAT*FtLb_s2Hp
    # transform tail rotor thrust, drag, torque to body frame
    Ftr=Tcant'*[-DTR;0.0;-TTR]
    Mtr_local=Tcant'*[0.0;0.0;-Qtr]
    Mtr=[0.0  -ZTR  YTR
         ZTR   0.0 -XTR
         -YTR  XTR  0.0]*Ftr+Mtr_local
    ptrNED = xf[10:12]+TNED2body'*[XTR;YTR;ZTR]

    # HORIZONTAL TAIL
    # horizontal tail index
    idxtail=1
    # equivalent airspeed used in stabilator schedule
    Veq=sqrt( (x[1])^2 + (x[2])^2 + (x[3])^2 )/1.688*sqrt(RHO/rhoSLSTD)
    ayg=(dx[1]-G*sin(x[7])*cos(x[8])+x[6]*x[1]-x[4]*x[3])/G
    StabInc=H60_StabSched(p[3],Veq,ayg,x[5],STABSET)
    # geometry [ft]
    xtail=-(FSTAIL[idxtail]-FSCG)/12.0
    ytail=(BLTAIL[idxtail]-BLCG)/12.0
    ztail=-(WLTAIL[idxtail]-WLCG)/12.0
    # external velocities in body frame
    Vextb=TNED2body*VextHT
    #
    qlossfac=sqrt(rotor_if[5])
    epsilon=rotor_if[7]
    sigma=0.0
    # local velocities
    uHT=(xf[1])*qlossfac+Vextb[1]+rotor_if[4]
    vHT=(xf[2]-D2R*xf[1]*sigma)+Vextb[2]
    wHT=(xf[3]-D2R*xf[1]*epsilon)*qlossfac+Vextb[3]+rotor_if[3]
    ublocal=uHT+(QB)*ztail-(RB)*ytail
    vblocal=vHT-(PB)*ztail+(RB)*xtail
    wblocal=wHT+(PB)*ytail-(QB)*xtail
    #
    if (CLELEV[idxtail] == 0.0)
        cthe=cos((ITAIL[idxtail]+StabInc)*D2R)
        sthe=sin((ITAIL[idxtail]+StabInc)*D2R)
    else
        cthe=cos(ITAIL[idxtail]*D2R)
        sthe=sin(ITAIL[idxtail]*D2R)
    end
    cphi=cos(PHITAIL[idxtail]*D2R)
    sphi=sin(PHITAIL[idxtail]*D2R)
    T=[cthe  sthe*sphi -sthe*cphi
       0     cphi       sphi
       sthe -sphi*cthe  cphi*cthe]
    # transfrom veclocities to orientation of tail
    temp=T*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]
    Vtot2=ut^2+vt^2+wt^2
    Vtot=sqrt(Vtot2)
    Vxz=sqrt(ut^2+wt^2)
    qt=0.5*RHO*Vtot2
    # tail angle of attack and sideslip
    alphat=atan2(wt,ut)*R2D
    sa=wt/Vxz
    ca=ut/Vxz
    sb=vt/Vtot
    cb=Vxz/Vtot
    # lift and drag coefficients
    CLt=table_lookup(ALTAIL[idxtail,:],CLTAIL[idxtail,:],alphat)+CLELEV[idxtail]*StabInc
    CDt=table_lookup(ALTAIL[idxtail,:],CDTAIL[idxtail,:],alphat)
    Lift=qt*STAIL[idxtail]*CLt
    Drag=qt*STAIL[idxtail]*CDt
    #
    Tw2t=[cb*ca -sb  -cb*sa
          sb*ca  cb  -sb*sa
          sa     0.0  ca]
    #
    temp=T'*Tw2t*[-Drag;0.0;-Lift]
    # horizontal tail forces and moments
    Fht=temp
    Mht=[ 0.0   -ztail  ytail
         ztail  0.0   -xtail
        -ytail  xtail  0.0]*Fht
    phtNED=xf[10:12]+TNED2body'*[xtail;ytail;ztail]

    # VERTICAL TAIL
    # vertical tail index
    idxtail=2
    # geometry [ft]
    xtail=-(FSTAIL[idxtail]-FSCG)/12.0
    ytail=(BLTAIL[idxtail]-BLCG)/12.0
    ztail=-(WLTAIL[idxtail]-WLCG)/12.0
    # external velocities in body frame
    Vextb=TNED2body*VextVT
    #
    qlossfac=sqrt(rotor_if[6])
    epsilon=0.0
    sigma=rotor_if[8]
    # local velocities
    uHT=(xf[1])*qlossfac+Vextb[1]+rotor_if[4]
    vHT=(xf[2]-D2R*xf[1]*sigma)+Vextb[2]
    wHT=(xf[3]-D2R*xf[1]*epsilon)*qlossfac+Vextb[3]+rotor_if[3]
    ublocal=uHT+(QB)*ztail-(RB)*ytail
    vblocal=vHT-(PB)*ztail+(RB)*xtail
    wblocal=wHT+(PB)*ytail-(QB)*xtail
    #
    if (CLELEV[idxtail] == 0.0)
        # set StabInc to 0 for vertical tail
        cthe=cos((ITAIL[idxtail]+0.0)*D2R)
        sthe=sin((ITAIL[idxtail]+0.0)*D2R)
    else
        cthe=cos(ITAIL[idxtail]*D2R)
        sthe=sin(ITAIL[idxtail]*D2R)
    end
    cphi=cos(PHITAIL[idxtail]*D2R)
    sphi=sin(PHITAIL[idxtail]*D2R)
    T=[cthe  sthe*sphi -sthe*cphi
       0     cphi       sphi
       sthe -sphi*cthe  cphi*cthe]
    # transfrom veclocities to orientation of tail
    temp=T*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]
    Vtot2=ut^2+vt^2+wt^2
    Vtot=sqrt(Vtot2)
    Vxz=sqrt(ut^2+wt^2)
    qt=0.5*RHO*Vtot2
    # tail angle of attack and sideslip
    alphat=atan2(wt,ut)*R2D
    sa=wt/Vxz
    ca=ut/Vxz
    sb=vt/Vtot
    cb=Vxz/Vtot
    # lift and drag coefficients (set incidence of vertical tail to zero)
    CLt=table_lookup(ALTAIL[idxtail,:],CLTAIL[idxtail,:],alphat)+CLELEV[idxtail]*0.0
    CDt=table_lookup(ALTAIL[idxtail,:],CDTAIL[idxtail,:],alphat)
    Lift=qt*STAIL[idxtail]*CLt
    Drag=qt*STAIL[idxtail]*CDt
    #
    Tw2t=[cb*ca -sb  -cb*sa
          sb*ca  cb  -sb*sa
          sa     0.0  ca]
    #
    temp=T'*Tw2t*[-Drag;0.0;-Lift]
    # horizontal tail forces and moments
    Fvt=temp
    Mvt=[ 0.0   -ztail  ytail
         ztail  0.0   -xtail
        -ytail  xtail  0.0]*Fvt
    pvtNED=xf[10:12]+TNED2body'*[xtail;ytail;ztail]

    # MAIN ROTOR
    # define inputs to module
    ur=zeros(23+NB*NSEG*3+1)
    ur[1]=controls[1]*pi/180
    ur[2]=controls[2]*pi/180
    ur[3]=controls[3]*pi/180
    ur[4:15]=xf[1:12]
    ur[16:21]=dx[1:6]
    # nominal rotor airspeed
    ur[22]=OMEGA
    ur[23]=0 #dx[IDXP[1]]
    # zero out external flow velocities in trim and linearization
    ur[24:23+NB*NSEG*3+1]=zeros(NB*NSEG*3+1,1)
    if ur[23]!=0.0
        soptval=1
    end
    # define blade twist at blade elements
    # twist=interp1(TWISTTABR,TWISTTABTHET,RSEG)
    itp = interpolate((TWISTTABR,),TWISTTABTHET, Gridded(Linear()))
    twist = itp(RSEG)
    # map to local variables names
    beta1c=xr[3]
    lambda0=xr[17]
    lambda1s=xr[18]
    lambda1c=xr[19]
    psi1=xr[20]
    # blade azimuth angles
    psi=psi1*ones(4)+collect(0:(NB-1))*(2*pi/NB)
    for i=1:NB
        psi[i]=mod(psi[i]',2*pi)
    end
    # trig functions
    cpsi = zeros(4)
    spsi = zeros(4)
    for i=1:NB
        cpsi[i]=cos(psi[i])
        spsi[i]=sin(psi[i])
    end
    # convert from IBC to MBC
    LMBC2IBC=[ones(4) (-1).^collect(1:NB) cpsi spsi]
    LMBC2IBCdot=OMEGA*[zeros(4,2) -spsi cpsi]
    beta=LMBC2IBC*xr[1:NB]
    beta_dot=LMBC2IBC*xr[NB+1:2*NB]+LMBC2IBCdot*xr[1:NB]
    zeta=LMBC2IBC*xr[2*NB+1:3*NB]
    zeta_dot=LMBC2IBC*xr[3*NB+1:4*NB]+LMBC2IBCdot*xr[2*NB+1:3*NB]
    #
    LIBC2MBC=1/NB*[ones(4)'
                   ((-1).^collect(1:NB))'
                   2*cpsi'
                   2*spsi']
    LIBC2MBCdot=OMEGA/NB*[zeros(2,4)
                           -2*spsi'
                           2*cpsi']
    LIBC2MBCddot=(OMEGA*OMEGA)/NB*[zeros(2,4)
                                   -2*cpsi'
                                   -2*spsi']
    # trig functions
    cbeta = zeros(4)
    sbeta = zeros(4)
    czeta = zeros(4)
    szeta = zeros(4)
    cpsizeta = zeros(4)
    spsizeta = zeros(4)
    for i=1:NB
        cbeta[i]=cos(beta[i])
        sbeta[i]=sin(beta[i])
        czeta[i]=cos(zeta[i])
        szeta[i]=sin(zeta[i])
        cpsizeta[i]=cos(psi[i]+zeta[i])
        spsizeta[i]=sin(psi[i]+zeta[i])
    end
    # control inputs
    theta1c=ur[1]
    theta1s=ur[2]
    theta0=ur[3]
    # air velocity inputs
    Vxb=ur[4]
    Vyb=ur[5]
    Vzb=ur[6]
    # angular rates
    p=ur[7]
    q=ur[8]
    r=ur[9]
    # Euler Angles
    sphi=sin(ur[10])
    cphi=cos(ur[10])
    sthe=sin(ur[11])
    cthe=cos(ur[11])
    spsiE=sin(ur[12])
    cpsiE=cos(ur[12])
    Tbody2NED=[cthe*cpsiE (sphi*sthe*cpsiE-cphi*spsiE) (cphi*sthe*cpsiE+sphi*spsiE)
               cthe*spsiE (sphi*sthe*spsiE+cphi*cpsiE) (cphi*sthe*spsiE-sphi*cpsiE)
               -sthe                  sphi*cthe                 cphi*cthe]
    # CG location
    xcg=ur[13]
    ycg=ur[14]
    zcg=ur[15]
    # linear accelerations
    Vxbdot=ur[16]
    Vybdot=ur[17]
    Vzbdot=ur[118]
    # angular acceleration inputs
    pdot=ur[19]
    qdot=ur[20]
    rdot=ur[21]
    # rotor acceleration
    Omega_dot=ur[23]
    # wind gusts
    VgNED=reshape(ur[24:23+NB*NSEG*3],3,NB*NSEG)
    # coupling gain to fade out inflow model
    InflowFade=ur[23+NB*NSEG*3+1]
    # convert gusts to shaft axes
    Vgs=Tshaft*Tbody2NED'*VgNED;
    Vgxs=reshape((Vgs[1,:])',NSEG,NB)'
    Vgys=reshape((Vgs[2,:])',NSEG,NB)'
    Vgzs=reshape((Vgs[3,:])',NSEG,NB)'
    # accelerations at Rotor Hub
    gx=G*sthe
    gy=-G*sphi*cthe
    gz=-G*cphi*cthe
    # moment arms from fuselage CG to main rotor hub
    Xh=(FSCGB-FSMR)/12.0
    Yh=(BLCGB-BLMR)/12.0
    Zh=(WLCGB-WLMR)/12.0
    # acceleration of the hub
    Vxhdot=Vxbdot-r*Vyb+q*Vzb-Xh*(q^2+r^2)+Yh*(p*q-rdot)+Zh*(p*r+qdot)+gx
    Vyhdot=Vybdot-p*Vzb+r*Vxb+Xh*(p*q+rdot)-Yh*(p^2+r^2)+Zh*(q*r-pdot)+gy
    Vzhdot=Vzbdot+p*Vyb-q*Vxb+Xh*(p*r-qdot)+Yh*(q*r+pdot)-Zh*(p^2+q^2)+gz
    # non-dimensional advance ratios
    muxh=(Vxb+q*Zh-r*Yh)/(OMEGA*R)
    muyh=(Vyb+r*Xh-p*Zh)/(OMEGA*R)
    muzh=(Vzb-q*Xh+p*Yh)/(OMEGA*R)
    # non-dimensional advance ratios derivatives
    muxhdot=(Vxbdot+qdot*Zh-rdot*Yh)/(OMEGA*R)
    muyhdot=(Vybdot+rdot*Xh-pdot*Zh)/(OMEGA*R)
    muzhdot=(Vzbdot-qdot*Xh+pdot*Yh)/(OMEGA*R)
    # shaft axes transformations
    temp=Tshaft*[Vxhdot;Vyhdot;Vzhdot]
    Vxsdot=temp[1]
    Vysdot=temp[2]
    Vzsdot=temp[3]
    #
    temp=Tshaft*[pdot;qdot;rdot]
    psdot=temp[1]
    qsdot=temp[2]
    rsdot=temp[3]
    rsdotmomd=rsdot-Omega_dot
    #
    temp=Tshaft*[muxh;muyh;muzh]
    muxs=temp[1]
    muys=temp[2]
    muzs=temp[3]
    mu=sqrt(muxs^2+muys^2)
    betawind = atan2(muys, muxs)
    #
    temp=Tshaft*[muxhdot;muyhdot;muzhdot]
    muxsdot=temp[1]
    muysdot=temp[2]
    muzsdot=temp[3]
    betawind_dot=(muxs*muysdot-muys*muxsdot)/(mu^2)
    #
    temp=Tshaft*[p;q;r]
    ps=temp[1]
    qs=temp[2]
    rs=temp[3]
    rminusO=rs-OMEGA
    # blade pitch
    thetab = zeros(4)
    thetab_dot = zeros(4)
    for i=1:NB
        thetab[i]=theta0+theta1c*cos(psi[i]+DELSP)+theta1s*sin(psi[i]+DELSP)-
                  tan(DELTA3)*beta[i]
        thetab_dot[i]=-theta1c*OMEGA*sin(psi[i]+DELSP)+theta1s*OMEGA*
                      cos(psi[i]+DELSP)
    end
    # dynamic Twist Component
    Veq=sqrt(RHO/rhoSLSTD)*sqrt(Vxb^2+Vyb^2+Vzb^2)/1.688
    KVDT=min(max(1.4e-4-(4.4e-6*Veq),-0.00052),-0.0003)
    theta_DTtip=KVDT*xr[21]
    theta_DT=vec(ones(4,1))*vec(DynTwistMode)'*theta_DTtip*pi/180.0
    # segment pitch
    theta=vec(thetab)*ones(1,NSEG)+ones(NB,1)*vec(twist)'*pi/180+theta_DT
    # gust velocities, transform to blade system
    Utg=1/(OMEGA*R)*( (spsizeta*ones(1,NSEG)).*Vgxs + (cpsizeta*ones(1,NSEG)).*Vgys)
    Urg=1/(OMEGA*R)*( ((-cpsizeta.*cbeta)*ones(1,NSEG)).*Vgxs +
        ((spsizeta.*cbeta)*ones(1,NSEG)).*Vgys - (sbeta*ones(1,NSEG)).*Vgzs )
    Upg=1/(OMEGA*R)*( ((-cpsizeta.*sbeta)*ones(1,NSEG)).*Vgxs +
        ((spsizeta.*sbeta)*ones(1,NSEG)).*Vgys + (cbeta*ones(1,NSEG)).*Vgzs )
    # for Pitt-Peters Model
    Upd_pp=(-lambda0*vec(cbeta)*ones(1,NSEG)-lambda1c*(e*vec(cbeta).*vec(cpsi)*ones(1,NSEG)+
           vec(cbeta).*vec(cpsizeta)*vec(RSEG)')-lambda1s*(e*vec(cbeta).*vec(spsi)*ones(1,NSEG)+vec(cbeta).*
           vec(spsizeta)*vec(RSEG)') )*(1-InflowFade)
    Urd_pp=(-lambda0*vec(sbeta)*ones(1,NSEG)-lambda1c*(e*vec(sbeta).*vec(cpsi)*ones(1,NSEG)+
           vec(sbeta).*vec(cpsizeta)*vec(RSEG)')-lambda1s*(e*vec(sbeta).*vec(spsi)*ones(1,NSEG)+vec(sbeta).*
           vec(spsizeta)*vec(RSEG)') )*(1-InflowFade)
    Upd_ph = zeros(NB,NSEG)
    Urd_ph = zeros(NB,NSEG)
    # blade segment velocities
    Up=(-muxs*vec(sbeta).*vec(cpsizeta)+muys*vec(sbeta).*vec(spsizeta)+muzs*vec(cbeta))*ones(1,NSEG)+
       (e/OMEGA)*(vec(cbeta).*(qs*vec(cpsi)+ps*vec(spsi))-vec(sbeta).*vec(szeta)*rminusO)*ones(1,NSEG)+
       (-vec(beta_dot)+qs*vec(cpsizeta)+ps*vec(spsizeta))*vec(RSEG)'/OMEGA+Upd_pp+Upg+Upd_ph
    Ut=(muxs*vec(spsizeta)+muys*vec(cpsizeta)-(e/OMEGA)*vec(czeta)*rminusO)*ones(1,NSEG)+
       ((vec(zeta_dot)-ones(4)*rminusO).*vec(cbeta)+(ps*vec(cpsizeta)-qs*vec(spsizeta)).*vec(sbeta))*vec(RSEG)'/OMEGA+
       Utg
    Ur=(muxs*vec(cbeta).*vec(cpsizeta)-muys*vec(cbeta).*vec(spsizeta)+muzs*vec(sbeta))*ones(1,NSEG)+
       (e/OMEGA)*(vec(sbeta).*(qs*vec(cpsi)+ps*vec(spsi))+vec(cbeta).*vec(szeta)*rminusO)*ones(1,NSEG)+
       Urd_pp+Urg+Urd_ph
    Utot=(Up.^2+Ut.^2+Ur.^2).^(1/2)
    # local angles of attack, Mach, skew
    cosgam=abs.(Ut)./((Ut.^2+Ur.^2).^(1/2))
    acosgam=abs.(cosgam)
    alphay=(180.0/pi)*atan2.((Ut.*tan.(theta)+Up).*acosgam,(Ut-Up.*tan.(theta).*
           cosgam.*cosgam))
    Mach=((Ut.^2+Up.^2).^(1/2))*OMEGA*R/VSOUND
    Machtab=max.(min.(Mach,1.0),0.0)
    # GenHel rotor model
    atrans0=cosgam.*alphay
    atranslow=(alphay+ones(4,10)*180.0).*cosgam-ones(4,10)*180.0
    atranshigh=ones(4,10)*180.0+(alphay-ones(4,10)*180.0).*cosgam
    #
    atrans1=((1.0-ACL1/90.0)*atrans0-ACL1*(ones(4,10)-cosgam))./(-ACL1/90.0*ones(4,10)+cosgam)
    atrans2=((1.0-ACL2/90.0)*atranshigh+ACL2*(ones(4,10)-cosgam))./(2.0*ones(4,10)-ACL2/90.0*ones(4,10)-cosgam)
    atrans3=((1.0+ACL3/90.0)*atrans0-ACL3*(ones(4,10)-cosgam))./(ACL3/90.0*ones(4,10)+cosgam)
    atrans4=((1.0+ACL4/90.0)*atranslow+ACL4*(ones(4,10)-cosgam))./(2.0*ones(4,10)+ACL4/90.0*ones(4,10)-cosgam)
    #
    atrans=atrans0.*(atrans0.<=ACL1).*(atrans0.>=ACL3).*(abs.(alphay).<90.0)+
           atrans1.*(atrans0.>ACL1).*(alphay.<90.0)+(atrans2.*(atranshigh.<ACL2)+
           atranshigh.*(atranshigh.>=ACL2)).*(alphay.>=90.0)+atrans3.*
           (atrans0.<ACL3).*(alphay.>=-90.0)+(atrans4.*(atranslow.>ACL4)+
           atranslow.*(atranslow.<=ACL4)).*(alphay.<-90.0)
    atransBtab=max.(min.(atrans,32.0),-32.0)
    alphayBtab=max.(min.(alphay,32.0),-32.0)
    #
    CL=interpp1(vec(AOAUTAB),vec(CLR0UTAB),atrans).*(abs.(atrans).>32.0)+
       interpp2(vec(AOABTAB),vec(MACHTAB),CLR0BTAB,atransBtab,Machtab).*
       (abs.(atrans).<=32.0)
    CD=interpp1(vec(AOAUTAB),vec(CDR0UTAB),alphay).*(abs.(alphay).>32.0)+
       interpp2(vec(AOABTAB),vec(MACHTAB),CDR0BTAB,alphayBtab,Machtab).*
       (abs.(alphay).<=32.0)+ones(4,10)*DCDMR
    # lift array specific to Peters-He inflow form
    BLSEGLIFT=((0.5*(ones(NB,1)*vec(CHORD)'))./R).*CL.*sign.(Ut).*Utot.*
              sqrt.(Ut.^2+Ur.^2)
    # blade segment forces
    Fp=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.*
       (CL.*Ut./acosgam+CD.*Up).*(ones(NB,1)*vec(delseg)')
    Ft=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.*
       (CD.*Ut-CL.*Up.*acosgam).*(ones(NB,1)*vec(delseg)')
    Fr=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.*
       (CD-CL.*Up.*acosgam./Ut).*Ur.*(ones(NB,1)*vec(delseg)')
    # total shear force at each blade hinge
    Fpb=zeros(4,1)
    Ftb=zeros(4,1)
    Frb=zeros(4,1)
    for i=1:NB
        Fpb[i]=sum(Fp[i,:])
        Ftb[i]=sum(Ft[i,:])
        Frb[i]=sum(Fr[i,:])
    end
    #Fpb=sum(Fp,2)
    #Ftb=sum(Ft,2)
    #Frb=sum(Fr,2)
    # sum total blade forces for dynamic twist model
    Fp_DT=sum(sqrt.(Fpb.^2+Ftb.^2))/NB
    # total aero moment about each hinge
    Mf_aero=zeros(4,1)
    Ml_aero=zeros(4,1)
    for i=1:NB
        Mf_aero[i]=R*sum(rseg.*Fp[i,:])
        Ml_aero[i]=R*sum(rseg.*Ft[i,:])
    end
    #Mf_aero=R*sum((ones(NB,1)*vec(rseg)').*Fp,2)
    #Ml_aero=R*sum((ones(NB,1)*vec(rseg)').*Ft,2)
    # aero moments
    Lha=-sum(Mf_aero.*spsizeta)   # maybe a mistake here?
    Mha=-sum(Mf_aero.*cpsizeta)
    # aero shear forces at each hinge
    Fxa=Frb.*cbeta.*szeta-Ftb.*czeta-Fpb.*sbeta.*szeta
    Fya=Frb.*cbeta.*czeta+Ftb.*szeta-Fpb.*sbeta.*czeta
    Fza=-Frb.*sbeta-Fpb.*cbeta
    # total aero thrust
    Tha=-sum(Fza)
    # aero force and momemnt coefficients (used by Pitt Peters model)
    CTa=Tha/(pi*RHO*OMEGA^2*R^4)
    CLa=Lha/(pi*RHO*OMEGA^2*R^5)
    CMa=Mha/(pi*RHO*OMEGA^2*R^5)
    # lag damper kinematics
    thetald=thetab-THETALDGEO*ones(4)
    xld=ALD*sbeta+ones(4)*CLD+(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*
        sin.(zeta+ones(4)*ZETA0)).*cbeta
    yld=-RLD*cos.(thetald)-BLD*sin.(zeta+ones(4)*ZETA0)+DLD*
        cos.(zeta+ones(4)*ZETA0)
    zld=ALD*cbeta-RLD*sin.(thetald)-(BLD*cos.(zeta+ones(4)*ZETA0)+
        DLD*sin.(zeta+ones(4)*ZETA0)).*sbeta
    daxld=sqrt.(xld.^2+yld.^2+zld.^2)
    #
    xld_dot = ALD*beta_dot.*cbeta-beta_dot.*sbeta.*(BLD*cos.(zeta+ones(4)
              *ZETA0)+DLD*sin.(zeta+ones(4)*ZETA0))+(-BLD*sin.(zeta+ones(4)*
              ZETA0)+DLD*cos.(zeta+ones(4)*ZETA0)).*zeta_dot.*cbeta
    yld_dot = RLD*sin.(thetald).*thetab_dot-(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*
              cos.(zeta+ones(4)*ZETA0)).*zeta_dot
    zld_dot = -ALD*sbeta.*beta_dot-RLD*cos.(thetald).*thetab_dot -cbeta.*
              beta_dot.*(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*sin.(zeta+ones(4)*
              ZETA0))-(-BLD*sin.(zeta+ones(4)*ZETA0)+DLD*cos.(zeta+ones(4)*
              ZETA0)).*zeta_dot.*sbeta
    raxld = (xld.*xld_dot+yld.*yld_dot+zld.*zld_dot)./daxld
    # lag damper force (linear lag damper)
    Fld=CLAG*raxld
    # flap and lag moments due to lag damper force
    Mf_ld=-Fld.*((zld*CLD)+(xld.*sin.(thetald)*RLD))./(12.0*daxld)
    Ml_ld=-Fld.*(RLD*cos.(thetald).*(xld.*cbeta-zld.*sbeta)+yld.*
          (CLD*cbeta+RLD*sin.(thetald).*sbeta))./(12.0*daxld)-
          KLAG*(zeta+ones(4)*ZETA0)
    # notional flap sping /damper
    Mf_fd=-KFLAP*beta-CFLAP*beta_dot
    # flap and lag acclerations
    beta_ddot = zeros(NB)
    zeta_ddot = zeros(NB)
    for i=1:NB
        beta_ddot[i]=(MBETA/IBETA)*(cbeta[i].*(Vzsdot+HOFFSET*(2*OMEGA*(ps*cpsi[i]-qs*spsi[i])+
              psdot*spsi[i]+qsdot*cpsi[i]))+sbeta[i].*czeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
              HOFFSET*rminusO^2)+sbeta[i].*szeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
              HOFFSET*rsdotmomd))+(cbeta[i].^2).*(czeta[i].*(psdot*spsi[i]+qsdot*
              cpsi[i]-2*(zeta_dot[i]+OMEGA).*(qs*spsi[i]-ps*cpsi[i]))-szeta[i].*(2*(OMEGA+
              zeta_dot[i]).*(ps*spsi[i]+qs*cpsi[i])+(qsdot*spsi[i]-psdot*cpsi[i])))+
              cbeta[i].*sbeta[i].*(2*zeta_dot[i]*rminusO-rminusO^2-zeta_dot[i].^2)+
              (Mf_aero[i]+Mf_ld[i]+Mf_fd[i])/IBETA
        zeta_ddot[i]=(MBETA./(IBETA*cbeta[i])).*(szeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
                  HOFFSET*rminusO^2)-czeta[i].*(Vxsdot*spsi[i]+Vysdot*cpsi[i]))+rsdotmomd+
                  (sbeta[i]./cbeta[i]).*(2*beta_dot[i].*(OMEGA+zeta_dot[i]-rs)+qsdot*spsizeta[i]-
                  psdot*cpsizeta[i])+2*beta_dot[i].*(czeta[i].*(qs*spsi[i]-ps*cpsi[i])+szeta[i].*
                  (ps*spsi[i]+qs*cpsi[i]))+(Ml_ld[i])./(IBETA*cbeta[i].^2)-
                  Ml_aero[i]./(IBETA*cbeta[i])
    end
    # inertial shear at each hinge
    Fxi = zeros(NB)
    Fyi = zeros(NB)
    Fzi = zeros(NB)
    for i=1:NB
        Fxi[i]=MBETA*(cbeta[i].*czeta[i].*(rsdotmomd-zeta_ddot[i])+2*sbeta[i].*czeta[i].*(zeta_dot[i].*
            beta_dot[i]-rminusO*beta_dot[i])+cbeta[i].*szeta[i].*(zeta_dot[i].^2+beta_dot[i].^2-2*
            rminusO*zeta_dot[i]+rminusO^2)+2*beta_dot[i].*cbeta[i].*(ps*cpsi[i]-qs*spsi[i])+
            sbeta[i].*szeta[i].*beta_ddot[i])-WBLADE/G*(Vxsdot*spsi[i]+Vysdot*cpsi[i])
        Fyi[i]=MBETA*(cbeta[i].*czeta[i].*(zeta_dot[i].^2+beta_dot[i].^2-2*rminusO*zeta_dot[i]+
            rminusO^2)+sbeta[i].*czeta[i].*beta_ddot[i]+cbeta[i].*szeta[i].*zeta_ddot[i]-2*beta_dot[i].*
            cbeta[i].*(ps*spsi[i]+qs*cpsi[i])+WBLADE*HOFFSET/(G*MBETA)*rminusO^2)+WBLADE/G*
            (Vxsdot*cpsi[i]-Vysdot*spsi[i])
        Fzi[i]=MBETA*(beta_ddot[i].*cbeta[i]-(beta_dot[i].^2).*sbeta[i] + 2*sbeta[i].*czeta[i].*
            beta_dot[i].*(ps*spsi[i]+qs*cpsi[i])+cbeta[i].*szeta[i].*(2*(OMEGA+zeta_dot[i]).*
            (ps*spsi[i]+qs*cpsi[i])+qsdot*spsi[i]-psdot*cpsi[i])-cbeta[i].*czeta[i].*(2*(OMEGA+
            zeta_dot[i]).*(ps*cpsi[i]-qs*spsi[i])+psdot*spsi[i]+qsdot*cpsi[i])+WBLADE*HOFFSET/
            (G*MBETA)*((2*OMEGA*(ps*cpsi[i]-qs*spsi[i])+psdot*spsi[i]+qsdot*cpsi[i])))-
            WBLADE/G*Vzsdot
    end
    # total shear force
    Fxt=Fxi+Fxa
    Fyt=Fyi+Fya
    Fzt=Fzi+Fza
    # total rotor forces and moments
    Thrust=-sum(Fzt)
    Hforce=sum(Fyt.*cpsi-Fxt.*spsi)
    Jforce=-sum(Fxt.*cpsi+Fyt.*spsi)
    Lhub=sum(HOFFSET*Fzt.*spsi+Mf_ld.*spsizeta)
    Mhub=sum(HOFFSET*Fzt.*cpsi+Mf_ld.*cpsizeta)
    Qhub=-sum(HOFFSET*Fxt-Ml_ld)
    # iteration to solve for quasi-steady lambda0
    lambda0_qs=lambda0
    delta_lambda=10.0
    iter=0
    while ((abs(delta_lambda)>1e-9) && (iter<70))
        lambda0_qsn=CTa/(2*sqrt(mu^2+(lambda0_qs-muzs)^2))
        delta_lambda=lambda0_qsn-lambda0_qs
        lambda0_qs=lambda0_qs+0.5*delta_lambda
        iter=iter+1
    end
    # angle of attack and mass flow parameters
    mutot=sqrt(mu^2+(muzs-lambda0_qs)^2)
    Vmpp=(mu^2+(lambda0_qs-muzs)*(2*lambda0_qs-muzs))/mutot
    # wind to shaft transformation
    cbetawind=muxs/mu
    sbetawind=muys/mu
    Ts2w=[1  0         0
          0  cbetawind sbetawind
          0 -sbetawind cbetawind]
    Tw2s=Ts2w'
    Tw2sdot=[0  0          0
             0 -sbetawind -cbetawind
             0 -cbetawind  sbetawind]*betawind_dot
    Ts2wdot=Tw2sdot'
    #
    lambda_w=Ts2w*[lambda0;lambda1s;lambda1c]
    #
    omega_w=Ts2w*[rs;ps;qs]/OMEGA
    beta1c_dot=xr[NB+3]
    beta1s_dot=xr[NB+4]
    beta0_dot=xr[NB+1]
    betastar_w=Ts2w*[beta0_dot;beta1s_dot;beta1c_dot]/OMEGA
    # Pitt-Peters with quasi-steady wake curvature effects.  Based on Zhao 2004
    Xpp=xr[4*NB+6]
    kcpp=xr[4*NB+7]
    kspp=xr[4*NB+8]
    # limit the denominator to be positive
    kcqs=(omega_w[3]-betastar_w[3])/max((lambda0_qs-muzs),0.001)
    ksqs=(omega_w[2]-betastar_w[2])/max((lambda0_qs-muzs),0.001)
    chi=atan2(mu,(lambda0_qs-muzs))
    #
    Xqs=tan(0.5*chi)
    #
    Xpp_dot=OMEGA*(15.0*pi*Vmpp)/32.0*(Xqs-Xpp)     # a bit off
    kcpp_dot=OMEGA*(15.0*pi*mutot)/32.0*(kcqs-kcpp)
    kspp_dot=OMEGA*(15.0*pi*mutot)/32.0*(ksqs-kspp)
    #
    Vpp=[mutot 0.0  0.0
         0.0   Vmpp 0.0
         0.0   0.0  Vmpp]
    Lpp=[0.5          0.0         -15*pi/64*Xpp
         0.0          2*(1+Xpp^2)  0
         15*pi/64*Xpp 0.0          2*(1-Xpp^2)]
    # use constant value of KR>=0, otherwise use airspeed schedule
    #if KR<0
    #    ciao=interp1(vec(KRVsched),vec(KRsched),Veq)
    #end
    # assume no wake curvature for now
    delL1=zeros(3,3)
    delL1[2,1]=0.5*KR*kspp
    delL1[3,1]=0.5*KR*kcpp
    delL1[1,2]=0.5*KR*kspp
    delL1[1,3]=0.5*KR*kcpp
    delL2=zeros(3,3)
    delL2[2,1]=0.75*KR*kspp*Xpp^2
    delL2[3,1]=-0.75*KR*kcpp*Xpp^2
    delL3=zeros(3,3)
    delL3[2,1]=1.25*mu*kcpp*Xpp*KR
    delL3[2,2]=(-2.5*kcpp*Xpp-1.5*mu*kspp*(1.0+1.5*Xpp^2))*KR
    delL3[2,3]=-2.5*kspp*Xpp*KR
    delL3[3,1]=1.25*mu*kspp*Xpp*KR
    delL3[3,2]=(-2.5*kspp*Xpp-1.5*mu*kcpp*(1.0-1.5*Xpp^2))*KR
    delL3[3,3]=-0.3*kcpp*Xpp*KR
    #
    ciao=OMEGA*[(75*pi)/128 0.0        0.0
                0.0         (45*pi)/16 0.0
                0.0         0.0        (45*pi)/16]
    ciao2=Ts2w*[CTa;-CLa;-CMa]-(inv(Lpp+delL1+delL2+delL3)*Vpp)*lambda_w
    lambda_dot_w=ciao*ciao2
    #lambda_dot_w=ciao*(Ts2w*[CTa;-CLa;-CMa]-Vpp/(Lpp+delL1+delL2+delL3)*lambda_w)
    lambda_dot=Tw2s*lambda_dot_w+Tw2sdot*lambda_w
    # build rotor derivatives vector
    xrdot=zeros(NRSTATES)
    xrdot[1:NB]=LIBC2MBC*beta_dot+LIBC2MBCdot*beta
    xrdot[NB+1:2*NB]=LIBC2MBC*beta_ddot+2*LIBC2MBCdot*beta_dot+LIBC2MBCddot*beta
    xrdot[2*NB+1:3*NB]=LIBC2MBC*zeta_dot+LIBC2MBCdot*zeta
    xrdot[3*NB+1:4*NB]=LIBC2MBC*zeta_ddot+2*LIBC2MBCdot*zeta_dot+LIBC2MBCddot*zeta
    xrdot[4*NB+1:4*NB+3]=lambda_dot
    xrdot[4*NB+4]=OMEGA
    # added state for filtered dynamic twist force
    xrdot[4*NB+5]=(1.0/TauDT)*(Fp_DT-xr[4*NB+5])
    xrdot[4*NB+6]=Xpp_dot
    xrdot[4*NB+7]=kcpp_dot
    xrdot[4*NB+8]=kspp_dot
    # rotor forces and moments in body frame at aircraft CG
    Fmrb=Tshaft'*[-Hforce;-Jforce;-Thrust];
    Mmrb=Tshaft'*[Lhub;Mhub;Qhub]+[0 -Zh Yh;Zh 0 -Xh;-Yh Xh 0]*Fmrb
    # main rotor forces and moments
    Fr=Fmrb
    Mr=Mmrb
    Qr=Qhub

    # EQUATIONS OF MOTION
    # total forces and moments
    Ftot=Fr+Ff+Fht+Fvt+Ftr
    Mtot=Mr+Mf+Mht+Mvt+Mtr
    Qreq=Qr+Qtr*GEARTR
    # equations of motion
    xfdot=eqnmot(xf,Ftot,Mtot)

    # STATE DERIVATIVE
    xdot=zeros(NSTATES)
    xdot[1:NFSTATES]=xfdot
    xdot[NFSTATES+1:NFSTATES+NRSTATES]=xrdot
    #xdot[NFSTATES+NRSTATES+1:NFSTATES+NRSTATES+NTRSTATES]=xtrdot[1]
    xdot[NSTATES]=xtrdot
    return xdot
end
