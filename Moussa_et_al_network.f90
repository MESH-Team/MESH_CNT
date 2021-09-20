subroutine diffusive_CNT(j,ntim,repeatInterval)

    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    !doubleprecision, intent(in) :: t0, t, tfin, saveInterval

    integer, intent(in) :: j, ntim, repeatInterval

    real,allocatable :: E_cnt(:,:), F_cnt(:,:) !,ini_E(:),ini_F(:)!, celerity(:,:), diffusivity(:,:) !, dx(:,:), ini_q(:), qp(:,:,:)
    !real :: USBoundary(2,100)

    integer :: i, n,kk!,repeat_interval

    real :: hi, gi, ki, pj, qj, rj, pj_p, qj_p, rj_p, sj_p, mi, ni, qp_ghost, r_interpol_time, qp_ghost_1
    real :: alpha_1, alpha_2, t_prime_1, t_prime_2
!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!

    ncomp = nx1(j)

    allocate(E_cnt(ncomp,ntim))
    allocate(F_cnt(ncomp,ntim))

    ! qp_ghost_1 is the estimated q value at the i = i-1 and n = ntim + 1
    ! at the first node of a reach, this is taken as equal to the same as the q of previous timestep and same location
    qp_ghost_1 = qp(1,ntim,j)

    E_cnt(1:ncomp,1) = ini_E(1:ncomp,j)
    F_cnt(1:ncomp,1) = ini_F(1:ncomp,j)

    !do n = 1, ntim
    !    i = 1
    !    if (qp(i,n,j) .lt. min_Q) then
    !        added_Q(i,n,j) = min_Q - qp(i,n,j)
    !        qp(i,n,j) = max(qp(i,n,j),min_Q)
    !    end if
    !end do

    do i = 2,ncomp

        do n = 2, ntim

            hi = dx(i-1,j) / (dtini * celerity(i-1,j))
            gi = diffusivity(i-1,j) * dx(i-1,j) / ((celerity(i-1,j)**3.0)*(dtini**2.0))
            ki = (gi ** 2.0)/(hi ** 2.0)
            mi = diffusivity(i-1,j) / (celerity(i-1,j)**2.0) * dx(i-1,j) / dtini
            ni = dx(i-1,j)

            if (n .eq. ntim) then
                ! applying ghost node !
                if (i .eq. 2) then
                    t_prime_1 = dx(i-1,j) / celerity(i-1,j)
                else
                    t_prime_1 = dx(i-2,j) / celerity(i-1,j)
                end if

                t_prime_2 = dx(i-1,j) / celerity(i-1,j)
                alpha_1 = t_prime_1 / dtini
                alpha_2 = t_prime_2 / dtini


                !print*, j, i, alpha_1, alpha_2


                !!! set 1 equation !!!
                pj = - 1.0/(2.0 * (1.0 + alpha_2)) * hi - gi / 2.0 - 2.0 * ki
                qj = 1.0 + (alpha_2 + 1.0)/(2.0 * alpha_2) * gi + (2.0*(alpha_2+1.0)/alpha_2) * ki
                rj = 1.0/(2.0*(alpha_2 + 1.0)) * hi - gi / (2.0 * alpha_2) - 2.0 / alpha_2 * ki

                pj_p = 1.0/(2.0 * (1.0 + alpha_1)) * hi + gi / 2.0 - 2.0 * ki
                qj_p = 1.0 - (alpha_1 + 1.0)/(2.0 * alpha_1) * gi + (2.0*(alpha_1+1.0)/alpha_1) * ki
                rj_p = - 1.0/(2.0*(alpha_1 + 1.0)) * hi + gi / (2.0 * alpha_1) - 2.0 / alpha_1 * ki

                !!! set 1 equation (with possible correction) !!!
                !pj = - 1.0/(2.0 * (1.0 + alpha_2)) * hi - gi / 2.0 /alpha_2 - 2.0 * ki / alpha_2
                !qj = 1.0 + (alpha_2 + 1.0)/(2.0 * alpha_2 * alpha_2) * gi + 2.0*(alpha_2+1.0)/(alpha_2*alpha_2) * ki
                !rj = 1.0/(2.0*(alpha_2 + 1.0)) * hi - gi / (2.0 * alpha_2*alpha_2) - 2.0 / (alpha_2*alpha_2) * ki

                !pj_p = 1.0/(2.0 * (1.0 + alpha_1)) * hi + gi / 2.0 / alpha_1 - 2.0 * ki / alpha_1
                !qj_p = 1.0 - (alpha_1 + 1.0)/(2.0 * alpha_1 * alpha_1) * gi + 2.0*(alpha_1+1.0)/(alpha_1 * alpha_1) * ki
                !rj_p = - 1.0/(2.0*(alpha_1 + 1.0)) * hi + gi / (2.0 * alpha_1*alpha_1) - 2.0 / alpha_1/alpha_1 * ki



                !!! set 2 equation !!!
                ! Equation set 2 is not working so far !
                !pj = -1.0/2.0*alpha_2/(1.0+alpha_2) * hi - 1.0/(1.0+alpha_2)*gi - 4.0/(1.0+alpha_2)*ki
                !qj = 1.0 + 1.0/2.0*(alpha_2 - 1.0)/alpha_2 * hi + 1.0/alpha_2 * gi + 4.0/alpha_2 * ki
                !rj = 1.0/2.0/alpha_2/(alpha_2+1.0) * hi - 1.0/alpha_2/(alpha_2+1.0) * gi - 4.0/alpha_2/(alpha_2+1.0)*ki

                !pj_p = 1.0/2.0*alpha_1/(1.0+alpha_1)*hi + 1.0/(1.0+alpha_1)*gi -4.0/(1.0+alpha_1)*ki
                !qj_p = 1.0 - 1.0/2.0*(alpha_1-1)/alpha_1 * hi - 1.0/alpha_1 * gi + 4.0 / alpha_1 * ki
                !rj_p = - 1.0/2.0/alpha_1/(1.0+alpha_1)*hi + 1.0/(1.0+alpha_1)*gi - 4.0/alpha_1/(1.0+alpha_1)*ki


                sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp_ghost_1 &
                !+ mi / 2.0 * (lateralFlow(i-1,n,j) - lateralFlow(i-1,n-1,j) ) &
                + ni * lateralFlow(i-1,n,j)

            else
                ! for any node other than ghost node !
                pj = - hi / 4.0 - gi / 2.0 - 2.0 * ki
                qj = 1.0 + gi + 4.0 * ki
                rj = hi / 4.0 - gi / 2.0 - 2.0 * ki

                pj_p = hi / 4.0 + gi / 2.0 - 2.0 * ki
                qj_p = 1.0 - gi + 4.0 * ki
                rj_p = -hi / 4.0 + gi / 2.0 - 2.0 * ki

                sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp(i-1,n+1,j) &
                    + mi / 2.0 * (lateralFlow(i-1,n+1,j) - lateralFlow(i-1,n-1,j) ) &
                    + ni * lateralFlow(i-1,n,j)
            end if

            !print*, pj, qj, rj, pj_p, qj_p, rj_p

            E_cnt(i,n) = -1.0 * rj / (pj * E_cnt(i,n-1) + qj)
            F_cnt(i,n) = ( sj_p - pj * F_cnt(i,n-1) ) / ( pj * E_cnt(i,n-1) + qj )

        end do

        !pause

        qp_ghost = qp(i-1,ntim,j)+lateralFlow(i-1,ntim,j)*dx(i-1,j)

        qp_ghost_1 = qp_ghost

        ! boundary at t=ntim, calculated from the ghost node at t=ntim+1
        qp(i,ntim,j) = E_cnt(i,ntim) * qp_ghost + F_cnt(i,ntim)

        do n = ntim-1,1,-1
            !print*, n, E_cnt(i,n), qp(i,n+1,j), F_cnt(i,n)
            qp(i,n,j) = E_cnt(i,n) * qp(i,n+1,j) + F_cnt(i,n)

            if (qp(i,n,j) .lt. min_Q) then
                added_Q(i,n,j) = min_Q - qp(i,n,j)
                qp(i,n,j) = max(qp(i,n,j),min_Q)
            end if

            ! change 20210715: provision of negative flow
            !if (abs(qp(i,n,j)) .lt. min_Q) then
            !    if (qp(i,n,j) .lt. 0.) then
            !        added_Q(i,n,j) = - min_Q - qp(i,n,j)
            !        qp(i,n,j) = - min_Q
            !    else
            !        added_Q(i,n,j) = min_Q - qp(i,n,j)
            !        qp(i,n,j) = min_Q
            !    end if
            !end if

        end do

    end do



    !do i = 1,ncomp
    !    ini_q_repeat(i,j) = max(ini_q_repeat(i,j),min_Q)
    !end do


    ! replacing with the initial value at all nodes
    qp(1:ncomp,1,j) = ini_q_repeat(1:ncomp,j)




    ! taking the new initial value for the next cycle
    ini_E(1:ncomp,j) = E_cnt(1:ncomp,repeatInterval+1)
    ini_F(1:ncomp,j) = F_cnt(1:ncomp,repeatInterval+1)
    ini_q_repeat(1:ncomp,j) = qp(1:ncomp,repeatInterval+1,j)


    deallocate (E_cnt, F_cnt)

end subroutine diffusive_CNT


subroutine mesh_diffusive_backward(j,ntim,repeatInterval)


    use constants_module
    use arrays_module
    use matrix_module
    use var_module
    use arrays_section_module
    use xsec_attribute_module
    use subtools

    implicit none

    integer, intent(in) :: j,ntim,repeatInterval

    real :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
    real :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope

    real :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, temp, dkdh
    real :: D_lim1, D_lim2, y_norm, y_crit, area_n, area_c, chnWidth, vel,y_norm1

    real :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ, stg1, stg2
    integer :: tableLength, jj, iii


    real :: elevTable_1(nel),areaTable_1(nel),rediTable_1(nel),convTable_1(nel),topwTable_1(nel),currentSquaredDepth_1(nel)
    real :: dKdATable_1(nel)
    real :: pereTable_1(nel),depthYi,tempDepthi_1,tempCo_1,temp_q_sk_multi_1,tempY_1,tempArea_1,tempRadi_1,tempbo_1,ffy
    real :: ffy_1, ffy_2, ffy_3, tempCo_2, tempCo_3, tempsfi_2, tempsfi_3
    real :: ffprime,tempDepthi_1_new,tempsfi_1,toll, dkda, tempPere_1, tempdKdA_1, tempY_2, tempY_3, temp_v_1, tempDepthi_2
    real :: skkkTable_1(nel), tempsk_1
    real :: usFroud, dsFroud, dUdt, eHds, y_alt, area_alt, y_cnj, area_cnj, y_crit_test, y_norm_test
    integer :: depthCalOk(ncomp), wlCalcMethod

    integer :: i, pp

!!++++++++++++++++++++ Diffusive wave Backward sweep starts +++++++++++++++++++!!

       !print*, j, 'normal depth', normalDepth(j)
        D_lim1 = -10.
        D_lim2 = -15.
        wlCalcMethod = 3

        ! wlCalcMethod = 1 : Mid-point bisection method
        ! wlCalcMethod = 2 : simple iterative method method with contribution from dU/dX
        ! wlCalcMethod = 3 : Newton Raphson method
        ! wlCalcMethod = 6 : Only Normal depth

        ncomp = nx1(j)

        depthCalOk(ncomp) = 1

        do i=ncomp,1,-1

            !print*, j, i, newY(i,j),newY(i,j)-z(i,j)

            currentQ = qp(i,repeatInterval+1,j)
            call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

    !      Calculating : read all attributes from tab file
            elevTable = xsec_tab(1,:,i,j)
            convTable = xsec_tab(5,:,i,j)
            areaTable = xsec_tab(2,:,i,j)
            pereTable = xsec_tab(3,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
            skkkTable = xsec_tab(11,:,i,j)
    !     interpolate the cross section attributes based on water elevation
            xt=newY(i,j)

            currentSquareDepth=(elevTable-z(i,j))**2.

            call r_interpol(currentSquareDepth,convTable,nel,(newY(i,j)-z(i,j))**2.0,co(i))
            if (co(i) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of conveyence was not possible, wl', &
                newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j), 'previous z',z(i+1,j)
                stop
            end if
			co(i) =q_sk_multi * co(i)

            call r_interpol(elevTable,areaTable,nel,xt,newArea(i,j))
            if (newArea(i,j) .eq. -9999) then
                print*, 'At j = ',j,', i = ',i, 'interpolation of newArea was not possible'
                stop
            end if
            call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
            call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
            call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))

            sfi = qp(i,repeatInterval+1,j) * abs(qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
            if (abs(qp(i,repeatInterval+1,j)) .lt. min_Q) then
                sfi = min_Q ** 2.0 * sign(1.0,qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
                ! at some head water basin, the Q boundary becomes 0.0 m3/s at some time. If Q = 0, then
                ! sfi = 0. This leads to diffusivity = NaN. To avoid NaN diffusivity, those sfi is calculated using the min_Q
                ! for those cases.
                ! print*, 'here'
            end if



            !print*, j, i, newY(i,j), sfi



            ! Celerity actual equation:
            !!! new for dkdh
            !do pp = 2,nel
            !    if (newY(i,j) .le. elevTable(pp)) then
            !        dkdh  =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
            !        EXIT
            !    endif
            !    if (pp .eq. nel) dkdh =(convTable(pp)-convTable(pp-1))/(elevTable(pp)-elevTable(pp-1))
            !end do
            !if (i .gt. 1) then
            !    dbdx(i)=(bo(i,j)-bo(i-1,j))/dx(i-1,j)
            !else
            !    dbdx(i)=dbdx(i+1)
            !end if

            ! Diffusivity actual equation:
            !diffusivity2(i) = co(i) * co(i) / 2.0 / qp(i,repeatInterval+1,j) / bo(i,j)

            !celerity2(i) = qp(i,repeatInterval+1,j) / co(i) * dkdh / bo(i,j) !+ dbdx(i) * diffusivity2(i) / bo(i,j)


            ! Calculating empirical celerity:
            !celerity2(i)=5.0 / 3.0 * sfi ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
            !celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6

             !min(rightBank(i,j)-leftBank(i,j), bo(i,j))
            chnWidth = rightBank(i,j)-leftBank(i,j)
            chnWidth = min(chnWidth,bo(i,j))
            !celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** 0.4 / chnWidth ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6


            if (depthCalOk(i) .eq. 1) then
                celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** &
                    0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6

                diffusivity2(i) = abs(qp(i,repeatInterval+1,j)) / 2.0 / bo(i,j) / abs(sfi)
                vel = qp(i,repeatInterval+1,j)/newArea(i,j)

                velocity(i,j) = vel






                ! Applying the upper limit of celerity
                if (celerity2(i) .gt. 3.0*abs(vel) ) celerity2(i) = abs(vel)*3.0
            else
                if (abs(qp(i,repeatInterval+1,j)) .lt. 1) then
                    celerity2(i)=0.5
                else
                    celerity2(i)=1.0
                end if

                diffusivity2(i)=diffusivity(i,j)
            end if


            ! Calculating Depth at the upstream node
            if (i .gt. 1) then
           !     print*, 'backward j', j,'i', i, 'routingNotChanged(i-1,j)',routingNotChanged(i-1,j), &
           !     'currentRoutingNormal(i-1,j)',currentRoutingNormal(i-1,j),'dimensionless_D(i-1,j)',dimensionless_D(i-1,j)

                slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                depthYi = newY(i,j) - z(i,j)

                tempDepthi_1 = oldY(i-1,j)-z(i-1,j)

                elevTable_1 = xsec_tab(1,:,i-1,j)
                areaTable_1 = xsec_tab(2,:,i-1,j)
                pereTable_1 = xsec_tab(3,:,i-1,j)
                rediTable_1 = xsec_tab(4,:,i-1,j)
                convTable_1 = xsec_tab(5,:,i-1,j)
                topwTable_1 = xsec_tab(6,:,i-1,j)
                dKdATable_1 = xsec_tab(9,:,i-1,j)


                !tempDepthi_1 = elevTable_1(100) - z(i-1,j)

                currentSquaredDepth_1=(elevTable_1-z(i-1,j))**2.

                toll = 1.0
                iii = 0

                dsFroud = abs(qp(i,repeatInterval+1,j))/sqrt(grav*newArea(i,j)**3.0/bo(i,j))

                !print*, qp(i,repeatInterval+1,j), newArea(i,j), bo(i,j), dsFroud*dsFroud; pause 111

                ! Applying mid point bisection
                if (wlCalcMethod .eq. 1) then
                tempY_1 = elevTable_1(2)
                tempY_2 = depthYi * 3. + z(i-1,j)
                tempY_3 = (tempY_1 + tempY_2) / 2.
                do while ( abs(toll) .gt. 0.001)
                    iii = iii +1

                    call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_1-z(i-1,j))**2.0,tempCo_1)
                    call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_2-z(i-1,j))**2.0,tempCo_2)
                    call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_3-z(i-1,j))**2.0,tempCo_3)

                    call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
                    tempCo_1 = tempCo_1 * temp_q_sk_multi_1
                    tempCo_2 = tempCo_2 * temp_q_sk_multi_1
                    tempCo_3 = tempCo_3 * temp_q_sk_multi_1

                    tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )
                    tempsfi_2 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_2** 2.0 )
                    tempsfi_3 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_3** 2.0 )

                    ffy_1 = (tempY_1-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)
                    ffy_2 = (tempY_2-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_2)
                    ffy_3 = (tempY_3-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_3)

                    if ((ffy_1 * ffy_2) .gt. 0.) then
                        tempY_2 = (tempY_2 - z(i-1,j)) * 2.0 + z(i-1,j)
                    elseif ((ffy_1 * ffy_3) .le. 0.) then
                        tempY_2 = tempY_3
                    elseif ((ffy_2 * ffy_3) .le. 0.) then
                        tempY_1 = tempY_3
                    end if
                    tempY_3 = (tempY_1 + tempY_2) / 2.0
                    toll = tempY_2 - tempY_1
                    tempDepthi_1 = tempY_3 - z(i-1,j)
                    depthCalOk(i-1) = 1


                end do
                end if

                ! Applying simple iteration method with inclusion of dU/dX
                if (wlCalcMethod .eq. 2) then
                vel = qp(i,repeatInterval+1,j)/newArea(i,j)
                toll = 1.0

                !qp(i,repeatInterval+1,j)=10.;depthYi = 1; vel = 0.5; sfi = 0.000309142; dx(i-1,j)=600;&
                !tempDepthi_1=1.;slope=0.001;qp(i-1,repeatInterval+1,j)=10.


                do while ( abs(toll) .gt. 0.001)
                    iii = iii +1
                    tempY_1 = tempDepthi_1 + z(i-1,j)

                    call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)

                    !print*, currentSquaredDepth_1
                    !print*, convTable_1
                    !print*, tempDepthi_1,(tempDepthi_1)**2.0
                    !print*, tempCo_1

                    call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
                    tempCo_1 = tempCo_1 * temp_q_sk_multi_1

                    call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)

                    !call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
                    !call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
                    !call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
                    !call r_interpol(elevTable_1,skkkTable_1,nel,tempY_1,tempsk_1)

                    !tempArea_1 = tempDepthi_1*20
                    !temp = tempArea_1/(20+2*tempDepthi_1)
                    !tempCo_1 = tempArea_1*(temp**(2./3.))/0.033


                    !print*, tempDepthi_1, tempArea_1,temp, tempCo_1; pause

                    temp_v_1 = qp(i-1,repeatInterval+1,j) / tempArea_1

                    tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

                    !tempDepthi_2 = depthYi - dx(i-1,j)*slope + 0.5*dx(i-1,j)*(sfi + tempsfi_1)
                    tempDepthi_2 = depthYi - dx(i-1,j)*slope + 0.5*dx(i-1,j)*(sfi + tempsfi_1)&
                        +0.5*(vel+temp_v_1)*(vel-temp_v_1)/grav


                   ! print*, depthYi, dx(i-1,j),slope,dx(i-1,j),(sfi + tempsfi_1),&
                   !     +(vel+temp_v_1),(vel-temp_v_1); pause
                    !print*, iii, tempDepthi_2

                    toll = tempDepthi_2 - tempDepthi_1
                    !print*, slope
                    !print*, depthYi, dx(i-1,j), slope, 0.5*(sfi + tempsfi_1); pause
                    !print*, iii,  depthYi, dx(i-1,j), slope, sfi, tempCo_1, tempsfi_1, vel, temp_v_1,tempDepthi_2; pause
                    !print*, iii,(sfi+tempsfi_1)/2.,depthYi,dx(i-1,j),slope,sfi,tempCo_1,tempsfi_1,vel,temp_v_1,tempDepthi_2
                    !pause

                    tempDepthi_1 = tempDepthi_2

                    print*, i, iii, tempCo_1, tempDepthi_1!; pause



                end do
                depthCalOk(i-1) = 1

                !print*, i, iii, tempDepthi_1

                end if

                ! Applying Newton–Raphson method corrected
                if (wlCalcMethod .eq. 3) then
                vel = qp(i,repeatInterval+1,j)/newArea(i,j)
                !print*, j, i, qp(i,repeatInterval+1,j),newArea(i,j), vel; pause 1000

                do while ( abs(toll) .gt. 0.001)
                    iii = iii +1
                    tempY_1 = tempDepthi_1 + z(i-1,j)

                    call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)
                    call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
                    tempCo_1 = tempCo_1 * temp_q_sk_multi_1

                    !print*, 'NewRap 1', j, i, iii

                    call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
                    call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
                    call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
                    call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
                    call r_interpol(elevTable_1,dKdATable_1,nel,tempY_1,tempdKdA_1)!    Change 20210520

                    !print*, 'NewRap 2', j, i, iii, tempArea_1, tempPere_1, tempRadi_1,tempbo_1,tempdKdA_1

                    tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

                    temp_v_1 = qp(i-1,repeatInterval+1,j)/tempArea_1

                    !dUdt = 0.5*(temp_v_1-oldQ(i-1,j)/oldArea(i-1,j))/dtini &
                    !      + 0.5*(vel-oldQ(i,j)/oldArea(i,j))/dtini

                    ffy = tempDepthi_1- depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1) !&
                        !+dUdt/grav*dx(i-1,j)  *0.0 &
                        !+0.5*(vel+temp_v_1)*(vel-temp_v_1)/grav  *0.0

                    !dkda=tempsk_1*((5.0/3.0*tempArea_1**(2.0/3.0)*tempPere_1)- &
                    !    (tempArea_1**(5.0/3.0)*2.0/tempbo_1))/(tempPere_1**(5.0/3.0))

                    dkda = tempdKdA_1

                    ffprime = 1 + dx(i-1,j) * tempbo_1 *  qp(i-1,repeatInterval+1,j) * &
                        abs(qp(i-1,repeatInterval+1,j)) / (tempCo_1 ** 3.0) * dkda

                    if (ffprime .eq. 0.) then
                        print*, 'ffprime = 0'
                        print*, j, i, qp(i-1,repeatInterval+1,j), tempCo_1
                        pause
                    end if

                    !print*, tempDepthi_1,depthYi,dx(i-1,j),slope,dx(i-1,j),(sfi + tempsfi_1),ffy
                    !print*, tempDepthi_1- depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1); pause

                    !print*, tempY_1,dx(i-1,j),newArea(i,j),tempsfi_1,sfi,ffy,ffprime,dkda,tempDepthi_1,depthYi; pause

                    tempDepthi_1_new = tempDepthi_1 - ffy / ffprime

                    tempDepthi_1_new = max(tempDepthi_1_new,0.005)

                    toll = abs(tempDepthi_1_new - tempDepthi_1)

                    if(iii .gt. 30)then
                        print*, 'Warning: Depth iteration reached maximum trial at j=', j, 'i=', i-1 , 'and',i, &
                        'depths are',tempDepthi_1, tempDepthi_1_new, 'slope=', slope, 'tempbo_1=',tempbo_1, 'dx=',&
                        dx(i-1,j), 'Q=', qp(i-1,repeatInterval+1,j)
                        print*, 'depth at d/s',depthYi
                        depthCalOk(i-1) = 0
                        EXIT
                    endif
                    tempDepthi_1 = tempDepthi_1_new
                    depthCalOk(i-1) = 1
                end do


                end if      ! calculation of Wl based on d/s wl ends

                !end if

                !pause 2020

                usFroud = abs(qp(i-1,repeatInterval+1,j))/sqrt(grav*tempArea_1**3.0/tempbo_1)

                newY(i-1,j) = tempDepthi_1 + z(i-1,j)

                ! calculating the normal depth at i=i-1 using the sfi as the slope
                !call normal_crit_y(i-1, j, q_sk_multi, sfi, max(abs(qp(i-1,repeatInterval+1,j)),min_Q), y_norm1, y_crit, &
                !    area_norm, area_crit)
                !print*, j, i-1,qp(i-1,repeatInterval+1,j)

                ! calculating the normal depth at i=i-1 using the slope between i and i-1
                ! Calculating the normal depth and critical depth at the river reach upstream as an output
                ! very small slope is creating an issue. Applying a min limit of bed slope.
                ! for Normal depth calculation, applying the absolute value of the Q
                call normal_crit_y(i-1, j, q_sk_multi, max(slope,0.0001), max(abs(qp(i-1,repeatInterval+1,j)),min_Q),&
                    y_norm, y_crit, area_norm, area_crit)

                ! to avoid sudden jump in water level calculation
                ! if the WL jumps abruptly, the tempsfi_1 would be very small compared to sfi.
                ! we are replacing the water level using normal water level
                if (sfi / tempsfi_1 .gt. 1e4)  newY(i-1,j) = y_norm



                !print*, j, i-1, newY(i-1,j),newY(i,j), y_norm

                !print*, j, i-1, tempDepthi_1, newY(i-1,j); pause

                !print*, 'Diffu', j, i, qp(i-1,repeatInterval+1,j)

                !print*, 'Normal', j, i

                ! no checking of conjugate depth. Applied after confirming with Goodwin creek.
                !if (y_norm .lt. y_crit) then
                !    call conjugate_depth(i-1, j, qp(i-1,repeatInterval+1,j), y_norm, y_cnj, area_cnj)

                !    !print*, i-1, newY(i-1,j)-z(i-1,j), y_norm-z(i-1,j), y_crit-z(i-1,j), y_cnj-z(i-1,j); pause 3030
                !    if (newY(i-1,j) .lt. y_cnj) then
                !        !newY(i-1,j) = y_norm                cmt 5
                !        !newArea(i-1,j) = area_norm              cmt 6
                !    end if
                !    !print*, i-1, newY(i-1,j), y_norm, y_crit, y_alt, z(i-1,j); pause
                !else
                !    if (usFroud .gt. 1.0) then      ! To avoid any wrong calculation of very shallow water depth when the channel is not steep
                !        !newY(i-1,j) = y_norm              cmt 7
                !        !newArea(i-1,j) = area_norm              cmt 8
                !    end if
                !end if

                !!!Test
                ! running Amite: temporarily turning off the auto switching
                !if (dimensionless_D(i-1,j) .lt. 1.0) newY(i-1,j) = y_norm
                !if (dimensionless_D(i-1,j) .lt. 1.0) newArea(i-1,j) = area_norm
                !if (qp(i-1,repeatInterval+1,j) .lt. 0.1) newY(i-1,j) = y_norm        !change 9
                !if (qp(i-1,repeatInterval+1,j) .lt. 0.1) newArea(i-1,j) = area_norm  ! change 10

                ! limiting very high value of WL
                if (newY(i-1,j) .gt. 10.0**5.) newY(i-1,j) = 10.0**5.
                !if (newY(i-1,j)-z(i-1,j) .gt. 10.0**2.) then
                !    call normal_crit_y(i-1, j, q_sk_multi, max(slope,0.0001), abs(qp(i-1,repeatInterval+1,j)), y_norm, y_crit, &
                !        area_norm, area_crit)
                !    print*, 'Warning: Depth calculated is too high trial at j=', j, 'i=', i-1 , 'and',i, &
                !    'depths are',newY(i-1,j)-z(i-1,j),newY(i,j)-z(i,j), 'slope=', slope, 'dx=',&
                !    dx(i-1,j), 'Q=', qp(i-1,repeatInterval+1,j),'norm depth=', y_norm-z(i-1,j),&
                !    'crit depth=',y_crit-z(i-1,j)
                !end if

                if (newY(i-1,j)-z(i-1,j) .le. 0.) then
                    print*, 'depth is negative at j= ', j,'i=',i-1,'newY=',(newY(jj,j),jj=1,ncomp)
                    print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
                    print*, 'newQ',(qp(jj,repeatInterval+1,j),jj=1,ncomp)
                    print*, 'Bed',(z(jj,j),jj=1,ncomp)
                    print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
                    pause 777
                end if

                !print*, 'All', j, i

            end if      ! end of if (i .gt. 1) || end of WL calculation at j reach

            ! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
            routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1

        end do

        !print*, j, newY(:,j); pause

        !celerity(1:ncomp,j) = minval(celerity2(1:ncomp)) ! change 20210423
        celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp  ! change 20210524

        !if (celerity(1,j) .lt. 0.5) celerity(1:ncomp,j) = 0.5

        do i = 1, ncomp
            celerity(i,j) = max(celerity(i,j),0.5)   !!! Applying Celerity lower limit
		end do

        diffusivity(1:ncomp,j)=sum(diffusivity2(1:ncomp)) / ncomp

		do i = 1, ncomp
			if (diffusivity(i,j) .gt. maxDiffuLm) diffusivity(i,j) = maxDiffuLm !!! Applying diffusivity upper limit
			if (diffusivity(i,j) .lt. minDiffuLm) diffusivity(i,j) = minDiffuLm !!! Applying diffusivity lower limit
		end do



		!!! test of min_X_plus and min_T_plus
		do i=1,ncomp-1
            if (celerity(i,j)/diffusivity(i,j)*dx(i,j) .le. 0.1) diffusivity(1:ncomp,j) = celerity(i,j)/0.1*dx(i,j)
            if (celerity(i,j)*celerity(i,j)/diffusivity(i,j)*dtini .le. 0.3) &
                diffusivity(1:ncomp,j) = celerity(i,j)*celerity(i,j)/0.3*dtini
            min_X_plus = min(min_X_plus,celerity(i,j)/diffusivity(i,j)*dx(i,j))
            min_T_plus = min(min_T_plus,celerity(i,j)*celerity(i,j)/diffusivity(i,j)*dtini)
        end do

!!++++++++++++++++++++ Diffusive wave Backward sweep ends +++++++++++++++++++!!




end subroutine mesh_diffusive_backward
