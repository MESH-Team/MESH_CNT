!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
program mesh

    use constants_module
    use arrays_module
    use arrays_section_module
    use var_module
    use matrix_module
    use sgate_module
    use xsec_attribute_module
    use subtools

    implicit none

    ! Local storage
    integer :: i, j, k, ppn, qqn, n, ntim,totalTimeSteps, igate, pp, boundaryFileMaxEntry, saveFrequency, repeatInterval
    integer :: linknb_ds, linknb_us
    integer :: lateralFLowAdditional ! accommodating this type for some scenarios where one node has two lateral flow connections
    real :: qnp1_ds, qnp1_us, qsum, y_ds, ini_time

    real :: cour, da, dq, x, saveInterval, width
    real :: qn, xt, maxCourant, dtini_given, nodenb, linknb
    real :: frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ, currentQ, area_ds
    real :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue, latFlowValue2

    real :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
    !doubleprecision ::  t, r_interpol_time, tfin, t1, t2, t0 !t0 start time

    integer :: tableLength, timestep, kkk
    real :: area_0, width_0, errorY, hydR_0, q_sk_multi, sumCelerity

    real :: r_interpo_nn

    character(len=128) :: output_path, other_input, ndep_path
    character(len=128) :: path

    !real, allocatable :: aa(:),bb(:),cc(:)

    !real :: pere(500)


    ! open file for input data



    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lower_Mississippi/input/input_BR_2_BC_2009.txt",status='unknown')
    !open(unit=1,file="../lateralFlow_test/input/input_dynamic_lateralFlow.txt",status='unknown')
    !open(unit=1,file="../Mississippi_River_11_years_20200511/input/input_Mississippi_BR2SWP_dynamic_new.txt",status='unknown')
    !open(unit=1,file="../Vermelion_River/input/input_Vermelion_dynamic_20200526.txt",status='unknown')
    !open(unit=1,file="../../../MESH_code/4-A-1_US_2bound/input_naturalChannel_tidal.txt",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/test.txt",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_tidal.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_8channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_11channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/input/input_multiJunction_NHD_mixedRouting_16channel.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/Large_NHD_Geometry/input_NHD_mixedRouting_29channel.txt",status='unknown')
    !open(unit=1,file="../Synthetic_Network_Test/input",status='unknown')
    !open(unit=1,file="../Rectangular_Y_Channel/input/input_naturalChannel_exact.txt",status='unknown')
    !open(unit=1,file="../Multijunction_Network/Large_NHD_Geometry_temp/input_NHD_mixedRouting_29channel_temp.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_exact_20201012.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_test.txt",status='unknown')
    !open(unit=1,file="../lateralFlow_test/input/input_crank_nicolson_test_lateralFlow.txt",status='unknown')
    !open(unit=1,file="../NHD_Y_Channel/input/input_naturalChannel_test_1chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Geometry\input_file_737",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Geometry\input_file_737_650_changed_temp.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Geometry\input_file_737_650_changed",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_LatFlow_from_structure.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_LatFlow_from_structure_one_channel.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_LatFlow_from_structure_&
    !    one_channel_intrpl.txt",status='unknown')
    !open(unit=1,file="D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\Rectangular_Y_Channel\input\input_naturalChannel_tidal_2.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_&
    !LatFlow_from_structure_interpol_shrt.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_chn5_shrt.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular3.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular1.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular1_5.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular1_3.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular1_4.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\input_rectangular1_2.txt",status='unknown')
    !open(unit=1,file="D:\One_Drive_Tulane\OneDrive - Tulane University\Kyle_edit_notWorking\Mesh_F_Kyle\&
    !Mississippi_River_11_years_20200511\input\input_Mississippi_BR2SWP_dynamic_new_run_20210123.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\CS_rectangular\input_rectangular.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_&
    !LatFlow_from_structure_interpol_diffu.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_dynamic_5_50m.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_dynamic_4_50m.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\&
    !CS_rectangular3\input_rectangular.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_&
    !LatFlow_from_structure_interpol_diffu_latQ_as_bound.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_LatFlow_from_structure_&
    !interpol_dyna_allRectang_test3.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\J_3_5\input_chn_3_5.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\&
    !RectangularCS\J_3_5\input_chn_3_5_netwrk.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\J_3_5\input_chn_3_5_analy.txt",status='unknown')
    !Open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\X-sec\RectangularCS\Test_20210310\input.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_added_LatFlow_&
    !from_structure_Right_channel_intrpl.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\&
    !input_ARBNM_added_LatFlow_from_structure_interpol_diffu.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_remapDx_dyna_allRectang_varWidth.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\routeLink_model\input_file_1",status='unknown')


    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\input_ARBNM_latQ_as_bound_3zones.txt",status='unknown')    ! surveyed section run
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\routeLink_model\&
    !    Geometry_RouteLink_1_2_3_4_5\input_orig.txt",status='unknown')    ! original NHD section run
    !open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\routeLink_model\Geometry_RouteLink_1_2_3_4_5\input.txt",status='unknown')    ! modified NHD section run

    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_737_final_interpolatedSections",status='unknown')     ! Florence final, interpolated cross sections
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_737_final_nonInterpolatedSections",status='unknown')     ! Florence two node chn, non-interpolated cross sections
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_737_final_interpolatedSectionsTest",status='unknown')      ! Florence Test case for diffu, celerity, roughness
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_2n\input_file_737",status='unknown')      ! Florence Test case for checking the number of interpolated sections
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_4n\input_file_737",status='unknown')      ! Florence Test case for checking the number of interpolated sections
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_5n\input_file_737",status='unknown')      ! Florence Test case for checking the number of interpolated sections
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_738_dummy",status='unknown')

    !open(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\input.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\input_Y_chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\input_multi_chn_flrnc.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\Florence_singleLine_ntwrk_test2\&
    !input_multi_chn_flrnc2.txt",status='unknown')


    !(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\input.txt",status='unknown')

    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    !input_file_added_reaches1982.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    !input_file_added_reaches1987.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    !input_file_added_reaches1997.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    !input_file_added_reaches1998.txt",status='unknown')

    !open(unit=1,file="D:\Project_Works\JTTI\codeTest_rectangle\compoundTrapoCS\input.txt",status='unknown')

    !open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_Data\input_file_1",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_DividedReach\input_file_4",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_DividedReach_2\input_file_7",status='unknown')  ! NHD Mississippi Final

    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    !moussa1982.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\compoundTrapoCS\input_Y_chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\Single_Trapo\input_single_chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_2n\&
    !input_file_737_Moussa_etal.txt",status='unknown')

    ! Test after applying CNT lat flow
    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Original_&
    !reachDistribution\input_file_9_CNT_original_ql_distribution.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_2n\&
    !input_file_737_Moussa_etal.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_Data_CNT_method_with_ql\input_file_1_CNT",status='unknown')  ! NHD Mississippi Final

    !! Tests after new ghost node !!
    open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\Single_Trapo\input_single_chn.txt",status='unknown')
    open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\compoundTrapoCS\input_Y_chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_Data_CNT_method_with_ql\input_file_1_CNT",status='unknown')  ! NHD Mississippi Final
    !open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Original_&
    !reachDistribution\input_file_9_CNT_original_ql_distribution.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\Test_interpolatedSection_2n\&
    !input_file_737_Moussa_etal.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\Single_Trapo\input_single_chn_nonUniDx.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\CNT_DongHa\input_Y_chn.txt",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\&
    !input_file_737_final_nonInterpolatedSections_Moussa",status='unknown')
    !open(unit=1,file="D:\Project_Works\JTTI\Vermilion_CNT\input_Vermelion_CNT.txt",status='unknown')

    !!! run for paper !!!
    open(unit=1,file="D:\Project_Works\JTTI\ARBNM\Model\routeLink_model\Geometry_RouteLink_1_2_3_4_5\&
    input_CNT.txt",status='unknown')    ! modified NHD section run
    open(unit=1,file="D:\Project_Works\JTTI\Goodwin_Creek_Experimental_Watershed\Devided_reaches\&
    input_file_added_reaches1982.txt",status='unknown')
    open(unit=1,file="D:\Project_Works\JTTI\Mississippi_NHD\NHD_Data_CNT_method_with_ql\input_file_1_CNT",status='unknown')  ! NHD Mississippi Final
    open(unit=1,file="D:\Project_Works\JTTI\Vermilion_CNT\input_Vermelion_CNT.txt",status='unknown')
    open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_737_final_nonInterpolatedSections_Moussa",status='unknown')

    !!! testing of new diffusivity limit as shown in Moussa et.at (1996)
    open(unit=1,file="D:\Project_Works\JTTI\Moussa_etal\Test_Case\compoundTrapoCS\input_Y_chn.txt",status='unknown')
    open(unit=1,file="D:\Project_Works\JTTI\Florence_NC\Model\input_file_737_final_nonInterpolatedSections_Moussa",status='unknown')


    print*, 'Reading input file'

    ! read data
    read(1,*) dtini_given     ! in seconds
    dtini = dtini_given; lastKnownDiffuDT = dtini_given         !; print*, dtini; pause 500
    read(1,*) dxini
    read(1,*) t0        ! in hours
    read(1,*) tfin      ! in hours
	totalTimeSteps = floor( (tfin - t0) / dtini * 3600)+1                 ! the min value of totalTimeSteps = 1

	!print*, totalTimeSteps; pause

	ntim = int(70.*60./dtini_given)+1
	!print*, ntim; pause

	read(1,*) nlinks
	allocate(nx1(nlinks))
	read(1,*) (nx1(i), i=1, nlinks)
    read(1,*) phi
    read(1,*) theta
    read(1,*) thetas
    read(1,*) thesinv
    read(1,*) alfa2
    read(1,*) alfa4
    read(1,*) f
    read(1,*) skk
	allocate(ini_y(nlinks))
	allocate(ini_q(nlinks))
	read(1,*) (ini_y(i), i=1, nlinks)!; print*, yy
	read(1,*) (ini_q(i), i=1, nlinks)
    read(1,*) cfl
    read(1,*) ots
    read(1,*) yw
    read(1,*) bw
    read(1,*) w
    read(1,*) option
    read(1,*) yn
    read(1,*) qn
    read(1,*) igate

    allocate(notSwitchRouting(nlinks))
    allocate(currentROutingDiffusive(nlinks))
    allocate(xSection_path(nlinks))
    allocate(bankLocation_path(nlinks))

    do i=1,nlinks
        read(1,*) xSection_path(i)
    end do

    do i=1,nlinks
        read(1,*) bankLocation_path(i)
    end do

    allocate(manning_strickler_path(nlinks))
    do i=1,nlinks
        read(1,*) manning_strickler_path(i)
    end do

    read(1,*) nupbnds                       ! No of u/s boundary data files
    allocate(upBoundTableEntry(nupbnds))
    allocate(upstream_path(nupbnds))
    do i=1,nupbnds
        read(1,*) upstream_path(i)
    end do

    read(1,*) ndnbnds                      ! No of d/s boundary data files
    allocate(downBoundTableEntry(ndnbnds))
    allocate(downstream_path(ndnbnds))
    do i=1,ndnbnds
        read(1,*) downstream_path(i)
    end do

    allocate(QSKtablePath(nlinks))
    do i=1,nlinks
        read(1,*) QSKtablePath(i)
    end do

    allocate(dx_path(nlinks))
    do i=1,nlinks
        read(1,*) dx_path(i)
    end do

    allocate(lateralFlow_path(nlinks))
    do i=1,nlinks
        read(1,*) lateralFlow_path(i)
    end do

    read(1,*) lateralFLowAdditional

    read(1,*) output_path
    read(1,*) option_dsbc
    read(1,*) maxTableLength
    read(1,*) nel
    read(1,*) timesDepth
    read(1,*) other_input
    read(1,*) boundaryFileMaxEntry
    read(1,*) saveInterval ; saveFrequency = saveInterval / dtini_given !; print*, saveFrequency; pause

    ! Reading lateral flow data starts
    allocate(noLatFlow(nlinks))
	read(1,*) (noLatFlow(i), i=1, nlinks)

    allocate(latFlowLocations(maxval(noLatFlow),nlinks))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(maxval(noLatFlow),nlinks))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(maxval(noLatFlow),nlinks))       ! no of x-secs at the downstream that the lateral flow is applied

    do j = 1,nlinks
        ncomp=nx1(j)
        if (noLatFlow(j) .gt. 0) then
            read(1,*) (latFlowLocations(i,j), i=1, noLatFlow(j))
        else
            read(1,*)
        end if
    end do
    do j = 1,nlinks
        ncomp=nx1(j)
        if (noLatFlow(j) .gt. 0) then
            read(1,*) (latFlowType(i,j), i=1, noLatFlow(j))
            do i=1,noLatFlow(j)
                if (latFlowType(i,j) .eq. 1) then
                    print*, 'Lateral flow at node = ', latFlowLocations(i,j), ', is a time series at reach ', j
                elseif (latFlowType(i,j) .eq. 2) then
                    print*, 'Lateral flow at node = ', latFlowLocations(i,j), ', is a function of upstream flow at reach ', j
                else
                    print*, 'Wrong lateral flow type is provided. Type ', latFlowType(i,j), 'is not a valid type at reach ', j
                    stop
                end if
            end do
        else
            read(1,*)
        end if
    end do

    do j = 1,nlinks
        ncomp = nx1(j)
        if (noLatFlow(j) .gt. 0) then
            read(1,*) (latFlowXsecs(i,j), i=1, noLatFlow(j))
        else
            read(1,*)
        end if
    end do
    ! Reading lateral flow data ends


    ! Reading Q-SK table data data starts
	allocate(noQSKtable(nlinks))                                ! how many tables are there in each river reach
	read(1,*) (noQSKtable(i), i=1, nlinks)

    allocate(eachQSKtableNodeRange(2,maxval(noQSKtable),nlinks))
    do j = 1,nlinks
        if (noQSKtable(j) .gt. 0) then
            read(1,*) (eachQSKtableNodeRange(1,i,j), i=1, noQSKtable(j))    ! upper limit of the river node number assigned to current table
            read(1,*) (eachQSKtableNodeRange(2,i,j), i=1, noQSKtable(j))    ! lower limit of the river node number assigned to current table
        !!! Need to test so that one section does not corresponds to more than one table
            do i = 2, noQSKtable(j)
                if ( eachQSKtableNodeRange(2,i-1,j) .ge. eachQSKtableNodeRange(1,i,j) ) then
                    print*, 'Wrong range of nodes applied for Q-Sk table.'
                    print*, 'Lower limit of Table ', i-1,'must be smaller than the upper limit of Table ', i, ' of reach', j
                    stop
                end if
            end do
        else
            read(1,*)
            read(1,*)
        end if
    end do
    ! Reading Q-SK table data data ends

    read(1,*) ndep_path ! Reading the location of network file

    read(1,*) applyNaturalSection  ! if 1, then attribute table will be activated, if 0, then rectangular channel will be applied
!print*, ndep_path; pause
    read(1,*) minDiffuLm
    read(1,*) maxDiffuLm
    close (1)       ! all input data read is finished

    ! Allocate arrays
    call setup_arrays(ntim, maxval(nx1), maxTableLength, boundaryFileMaxEntry, maxval(noLatFlow), maxval(noQSKtable), nlinks)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, maxval(nx1),nlinks)

    dt = dtini
    minDx = 1e6

    do j = 1,nlinks
        ncomp = nx1(j)
        open(unit=90, file=trim(dx_path(j)), status='unknown')
            do i=1,ncomp-1
                read(90, *) x, dx(i,j)
            end do
            if (minval(dx(1:ncomp-1,j)) .le. minDx) minDx=minval(dx(1:ncomp-1,j))
        close(90)
        print*, j, 'dx', dx(1:ncomp-1,j)
    end do

  ! reading Strickler's coefficient at each section
  !  do j = 1,nlinks
  !      ncomp = nx1(j)
  !      open(unit=85,file=trim(manning_strickler_path(j)), status='unknown') !! //'Mannings_Stricklers_coeff.txt', status='unknown')
  !      do i=1,ncomp
  !          read(85, *) x, sk(i,j)
  !          call readXsection(i,(1.0/sk(i,j)),timesDepth,j)
  !          ! This subroutine creates attribute table for each cross sections and saves in the hdd
  !          ! setting initial condition
  !          oldY(i,j) = ini_y(j) + z(i,j)
  !          !if (oldY(i,j) .lt. 0.) oldY(i,j) = 0.
  !          oldQ(i,j) = ini_q(j)
  !      end do
  !      close(85)
  !  end do

    ! reading Strickler's coefficient at each section
    print*, 'Reading Geometry'



    ! reading bank locations
    if (applyNaturalSection .eq. 1) then
        do j = 1,nlinks
            ncomp = nx1(j)
            open(unit=12,file=trim(bankLocation_path(j)), status='unknown')  !! read bed level
            do i=1,ncomp
                read(12, *) x, leftBank(i,j), rightBank(i,j)
            end do
            close(12)
        end do
    end if


    do j = 1,nlinks
        ncomp = nx1(j)
        open(unit=85,file=trim(manning_strickler_path(j)), status='unknown') !! //'Mannings_Stricklers_coeff.txt', status='unknown')



!       Skip headings
        !read(11,*)
        if (applyNaturalSection .eq. 0) then
            open(unit=11,file=trim(xSection_path(j)), status='unknown')  !! read bed level
            do i=1,ncomp
                read(85, *) x, skMain(i,j)
                read(11, *) x, z(i,j) ,bo(i,j)
                oldY(i,j) = ini_y(j) + z(i,j)
                oldQ(i,j) = ini_q(j)
            end do
        else

            do i=1,ncomp
                read(85, *) x, skLeft(i,j), skMain(i,j), skRight(i,j)

                !print*, x, skLeft(i,j), skMain(i,j), skRight(i,j)
                !print*, i,(1.0/skLeft(i,j)),(1.0/skMain(i,j)),(1.0/skRight(i,j)),leftBank(i,j), rightBank(i,j),timesDepth,j
                call readXsection(i,(1.0/skLeft(i,j)),(1.0/skMain(i,j)),(1.0/skRight(i,j)),&
                    leftBank(i,j), rightBank(i,j),timesDepth,j)
                oldY(i,j) = ini_y(j) + z(i,j)
                oldQ(i,j) = ini_q(j)
            end do
        end if
        close(85)

        print*, j, 'bed', z(1:ncomp,j)
        !print*, j, 'width', bo(1:ncomp,j)
        print*, j, 'initial_wl', oldY(1:ncomp,j)
    end do

    NAnum = -100
    ityp = 1

    ! setting initial condition
    ! setting initial condition from previous work
    !open(unit=91,file=trim(output_path)//'initialCondition.txt', status='unknown')
    ! read(91, *)
    !do i=1,ncomp
    !    read(91, *) oldQ(i), oldY(i), oldArea(i)
    !end do
    !close(91)

    !q(1, :) = qq
    !oldQ = qq

    ! reading Q-Strickler's coefficient multiplier table
    do j = 1,nlinks
        ncomp = nx1(j)
        do i=1,noQSKtable(j)
            write(file_num,'(i4.4)')i
            open(86,file=trim(QSKtablePath(j))//'Q_Mannings_table_'//file_num//'.txt')
            do n=1,maxTableLength
                read(86,*,end=300) Q_sk_Table(1, n, i, j), Q_sk_Table(2, n, i, j)
            end do
300         close(86)
            Q_sk_tableEntry(i,j) = n-1
        end do
    end do



    !!! This part is the code from DongHa for network and then modified by Nazmul
    open(unit=80, file=trim(ndep_path), status='unknown')
    ! this file contains the information for network connectivity.
    ! left column shows river reach number (j) and
    ! the right column shows  the number of links that are immediately upstream of link j
    ! skip one header line
    read(80, *)
    do j = 1,nlinks
        read(80, *) x, ndep(j)
    end do

    ! ++++ Y channel connectivity ++++++!
    !* the number of links that are immediately upstream of link j
    !ndep(1)=0; ndep(2)=0; ndep(3)=0; ndep(4)=0; ndep(5)=2; ndep(6)=3; ndep(7)=1; ndep(8)=1


    allocate(uslinks(maxval(ndep),nlinks))
    uslinks = NAnum
    ! Next part of the ndep file reads the link number of k_th link that is immediately upstream of link j
    ! The first column is the link number j, and the rest of the columns are the k_th upstream link of j
    ! skip one header line or blank line
    read(80, *)
    do j = 1,nlinks
        read(80,*) x, (uslinks(i,j), i=1, ndep(j)) !; print*, x, (uslinks(i,j), i=1, ndep(j))
    end do
    !pause

    ! Next part of the ndep file reads the link number of k_th link that is immediately downstream of link j
    ! skip one header line or blank line
    dslink = 0
    read(80, *)
    do j = 1,nlinks
        read(80,*) x, dslink(j) !; print*, x, dslink(j)
    end do
    !pause

    ! Next part of the ndep file reads the boundary condition of j_th link
    ! the first coulmn is the link number j, the 2nd column indicates the u/s condition of j,
    ! and the 3rd column indicates the d/s condition of j
    !* when data at either upper or lower end of link j is available,
    !* instrdflag(j,1)=1 when water level data is known at the upper end of link j
    !* instrdflag(j,1)=2 when discharge data is known
    !* instrdflag(j,1)=3 when rating curve is known

    !* instrdflag(j,2)=1 when water level data is known at the lower end of link j
    !* instrdflag(j,2)=2 when discharge data is known
    !* instrdflag(j,2)=3 when rating curve is known
    !* Otherwise, instrdflag(j,1/2)=0

    ! skip one header line or blank line
    read(80, *)
    do j = 1,nlinks
        read(80, *) x, instrdflag(j,1), instrdflag(j,2)  !; print*, x, instrdflag(j,1), instrdflag(j,2)
    end do
    !pause

    close(80)



    ! Read hydrograph input Upstream
    do j = 1, nupbnds
        open(unit=87, file=upstream_path(j))
        do n=1,boundaryFileMaxEntry
            read(87,*,end=301) USBoundary(1, n,j), USBoundary(2, n,j)       !! time column is in minutes
        end do
301     close(87)
        upBoundTableEntry(j) = n-1   ! stores the number of entries in the boundary file
    end do
    !pause

    ! Read hydrograph input Downstream
    do j = 1, ndnbnds
        open(unit=88, file=downstream_path(j))
        do n=1,boundaryFileMaxEntry
          read(88,*,end=302)  DSBoundary(1, n, j), DSBoundary(2, n, j)       !! time column is in minutes
        end do
302     close(88)
        downBoundTableEntry(j) = n-1   ! stores the number of entries in the boundary file
    end do

    t=t0*60.0     !! from now on, t is in minute

    ! applying boundary
    ! interpolation of boundaries at the initial time step
    !! Need to define which channel has terminal boundary
    ppn = 1; qqn = 1        ! ppn and qqn indicates the sequence no of the boundary data
    do j = 1, nlinks
        ncomp = nx1(j)
        if (instrdflag(j,1) .eq. 2) then
            ! interpolation of boundaries at the desired time step
            newQ(1,1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t)
            ppn = ppn +1
        end if
        if (instrdflag(j,2) .eq. 1) then
            ! interpolation of boundaries at the desired time step
            oldY(ncomp,j)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
                DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),t)
            qqn = qqn +1
        end if
    end do

    !oldY = newY


 !   if (instrdflag(j,1) .eq. 2) then        !! I.e. for river 1 and 2
 !       ! interpolation of boundaries at the desired time step at upstream Q boundaries from given time series
 !       newQ(1,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
 !           USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),t+dtini/60.)
 !       ppn = ppn +1
 !   end if

! correcting the WL initial condition based on the WL boundary
! so that the initial WL is higher than or equal to the WL boundary, at j = nlinks, i=ncomp
    print*, 'Correcting initial WL based on d/s boundary', oldY(nx1(nlinks),nlinks)
    do j = 1,nlinks
        ncomp = nx1(j)
        do i=1,ncomp
            oldY(i,j) = max(oldY(i,j),oldY(nx1(nlinks),nlinks))     ! corrected 20210524
        end do
        print*, j, (oldY(i,j),i=1,ncomp)
    end do


    ! read lateral flow conditions
    do j=1, nlinks
        ncomp = nx1(j)
        do i=1,noLatFlow(j)
            write(file_num,'(i4.4)')latFlowLocations(i,j)
            open(89,file=trim(lateralFlow_path(j))//'lateral_'//file_num//'.txt')
            read(89, *)     ! skipping the header row
            do n=1,boundaryFileMaxEntry
                read(89,*,end=303) lateralFlowTable(1, n, i,j), lateralFlowTable(2, n, i,j)
            end do
303         close(89)
            dataInEachLatFlow(i,j) = n-1
        end do
    end do
    !pause






    ! Open files for output
    path = trim(output_path) // 'output_wl.txt'
    open(unit=8, file=trim(path), status='unknown')

    path = trim(output_path) // 'q.txt'
    open(unit=9, file=trim(path), status='unknown')

    path = trim(output_path) // 'celerity.txt'
    open(unit=95, file=trim(path), status='unknown')

    path = trim(output_path) // 'diffusivity.txt'
    open(unit=96, file=trim(path), status='unknown')

    path = trim(output_path) // 'depth.txt'
    open(unit=991, file=trim(path), status='unknown')

    path = trim(output_path) // 'velocity.txt'
    open(unit=993, file=trim(path), status='unknown')

    path = trim(output_path) // 'added_Q.txt'
    open(unit=994, file=trim(path), status='unknown')

    path = trim(output_path) // 'courant.txt'
    open(unit=9921, file=trim(path), status='unknown')


	! Some essential initial parameters for Diffusive Wave
	theta = 1.0
	qpx = 0.

    width = 10. !!! average width of MS from BR to BC
    celerity = -999.
    maxCelerity = 1.0
    diffusivity = -999.
    maxCelDx = maxCelerity / minDx      ! Change 20210408

    added_Q = 0.


    !!! setting initial values of dimensionless parameters
    !dimensionless_Cr, dimensionless_Fo, dimensionless_Fi, dimensionless_Fc, dimensionless_Di, dimensionless_D
    dimensionless_Fi = 10.1
    dimensionless_Fc = 10.1
    dimensionless_D  = 0.1
    currentROutingDiffusive = 1

    ! parameters for diffusive vs partial diffusive
    currentRoutingNormal = 0
    routingNotChanged = 0

    ! initialization of Q, celerity and diffusivity
    do j=1,nlinks
        ncomp = nx1(j)
        newQ(1:ncomp,1,j) = ini_q(j)
        celerity(1:ncomp,j) = 1.0   !0.93 !real(j)
        diffusivity(1:ncomp,j) = 100.    !50. ! 500*real(j)
    end do

    min_Q = 0.028316    ! 1 cfs = 0.028316 m3/s
    !min_Q = 0.0028316

    call cpu_time( t1 )
        ! Output initial conditions
    do j=1, nlinks
        ncomp = nx1(j)

        write(8, 10)  t, j, (oldY(i,j), i=1,maxval(nx1))
        !write(9, 10)  t, j, (newQ(i,1,j), i=1, maxval(nx1))
        !write(51, 10) t, j, (oldArea(i,j), i=1, maxval(nx1))
        !write(882, 10) t, j, (oldQ(i,1,j)/oldArea(i,j), i=1, maxval(nx1))
        !write(941, 10) t, j, (dimensionless_Cr(i,j), i=1, maxval(nx1)-1)
        !write(942, 10) t, j, (dimensionless_Fo(i,j), i=1, maxval(nx1)-1)
        !write(943, 10) t, j, (dimensionless_Fi(i,j), i=1, maxval(nx1)-1)
        !write(944, 10) t, j, (dimensionless_Di(i,j), i=1, maxval(nx1)-1)
        !write(945, 10) t, j, (dimensionless_Fc(i,j), i=1, maxval(nx1)-1)
        !write(946, 10) t, j, (dimensionless_D(i,j), i=1, maxval(nx1)-1)
        write(95, 10) t, j, (celerity(i,j), i=1, maxval(nx1))
        write(96, 10) t, j, (diffusivity(i,j), i=1, maxval(nx1))
        !write(97, *) t, j, currentROutingDiffusive(j)
        !write(98, *) t, j, (currentRoutingNormal(i,j), i=1, maxval(nx1)-1)
        !write(99, *) t, j, (routingNotChanged(i,j), i=1, maxval(nx1)-1)
        write(991, 10) t, j, (oldY(i,j)-z(i,j), i=1,maxval(nx1))
        !write(992, *) t, j, (courant(i), i=1, maxval(nx1)-1)
        !write(9921, *) t, j, (courant(i), i=1, maxval(nx1))

        !write(881, 10) t,j, (normalDepthAtNodes(i,j), i=1, maxval(nx1))

        write(993, 10) t, j, (celerity(i,j), i=1, maxval(nx1))

        write(9921, 10) t, j, (celerity(i,j), i=1, maxval(nx1))
    end do

    !do j=1,nlinks
    !    print*, (z(i,j),i=1,maxval(nx1))
    !end do
    !pause

!    frus2 = 9999.
!    notSwitchRouting=0
!    minNotSwitchRouting = 10000         ! works between Dynamic and Diffusive switching
!    minNotSwitchRouting2 = 000        ! works between full Diffusive and partial Diffusive switching

    !

    ini_E = 1.0
    ini_F = 0.0

    do j = 1, nlinks
        ncomp = nx1(j)
        ini_q_repeat(1:ncomp,j) = ini_q(j)
    end do



    ! repeat time = 1hr
    !repeatInterval = 60
	repeatInterval = int(60.*60./dtini_given)


    !!! test of min_X_plus and min_T_plus
    min_X_plus = 999.
    min_T_plus = 999.


    !! calculation loop starts !!
    do kkk = 1,totalTimeSteps-1, repeatInterval

    ini_time = real(kkk-1)*dtini/60.+t0*60.

    ! Loop in space
    do j = 1, nlinks

        ppn = 1
        qqn = 1

        ncomp = nx1(j)

        lateralFlow(:,:,j) = 0

        ! applying the boundary conditions in the matrix newQ
        do timestep=1, ntim
            if (instrdflag(j,1) .eq. 2) then        !! I.e. for river 1 and 2
                ! interpolation of boundaries at the desired time step at upstream Q boundaries from given time series
                newQ(1,timestep,j)=r_interpol_time(USBoundary(1, 1:upBoundTableEntry(ppn), ppn), &
                    USBoundary(2, 1:upBoundTableEntry(ppn), ppn),upBoundTableEntry(ppn),ini_time+dtini/60.*(timestep-1))

            end if

            !!START+++++++ If the channel has boundary originated from a junction+++++++
            !+++----------------------------------------------------------------
            !+ Hand over water from upstream to downstream properly according
            !+ to the nature of link connections, i.e., serial or branching.
            !+ Refer to p.52,RM1_MESH
            !+++----------------------------------------------------------------
            if (ndep(j).gt.0) then  !     !* the number of links that are immediately upstream of link j. Example: j = 3
                !*total water discharge at n+1 at the end nodes of upstream links that join link j
                newQ(1,timestep,j) = 0.
                do k=1, ndep(j)
                    linknb=uslinks(k,j); nodenb=nx1(linknb)
                    newQ(1,timestep,j) = newQ(1,timestep,j) + newQ(nodenb,timestep,linknb)

                end do
            end if
          !!END+++++++ If the channel has boundary originated from a junction+++++++


            do i=1,noLatFlow(j)
                latFlowValue = r_interpol_time(lateralFlowTable(1, 1:dataInEachLatFlow(i,j), i,j), &
                    lateralFlowTable(2, 1:dataInEachLatFlow(i,j), i,j),dataInEachLatFlow(i,j),ini_time+dtini/60.*(timestep-1))



                ! added condition for lateral flow at the upstream boundary node
                if (latFlowLocations(i,j) .eq. 1) then
                    latFlowValue = latFlowValue / dx(1,j)
                    lateralFlow(1,timestep,j)=latFlowValue
                    newQ(1,timestep,j) = newQ(1,timestep,j)+lateralFlow(1,timestep,j)*dx(1,j)
                else
                    latFlowValue = latFlowValue / dx(latFlowLocations(i,j),j)
                    lateralFlow(latFlowLocations(i,j),timestep,j)=latFlowValue
                end if
                !print*, j, i,latFlowLocations(i,j), timestep, lateralFlow(latFlowLocations(i,j),timestep,j); pause 1003
            end do

        end do

        lateralFlow(1,:,j) = 0.       ! lateral flow at i = 1 is already added to the u/s

        qp(:,:,j) = newQ(:,:,j)

        call diffusive_CNT(j,ntim,repeatInterval)

        newQ(:,:,j) = qp(:,:,j)

    end do

    do timestep = 1,repeatInterval,saveFrequency
        do j=1,nlinks

            write(9,  10) ini_time+real(timestep-1)*dtini/60.,j, (newQ(i,timestep,j), i=1, maxval(nx1))
            write(994,10) ini_time+real(timestep-1)*dtini/60.,j, (added_Q(i,timestep,j), i=1, maxval(nx1))
        end do
    end do

    added_Q = 0.


    ! running the backward diffusive to calculate the depth, celerity and diffusivity
    do j=nlinks,1,-1

        ncomp = nx1(j)

        ! applying downstream water level at the end time of the current loop
        if (instrdflag(j,2) .eq. 1) then    ! water level is known at the downstream
            ! interpolation of boundaries at the desired time step at downstream WL boundaries from given time series
            newY(nx1(nlinks),nlinks)=r_interpol_time(DSBoundary(1, 1:downBoundTableEntry(qqn), qqn), &
            DSBoundary(2, 1:downBoundTableEntry(qqn), qqn),downBoundTableEntry(qqn),ini_time+dtini/60.*(repeatInterval))
        elseif (instrdflag(j,2) .eq. 0) then    ! water level is calculated from the downstream river reach
            linknb=dslink(j)    ! Which river reach is immediately downstream of reach j
            newY(ncomp,j)= newY(1,linknb)   ! taking the WL from the d/s reach
        end if

        call mesh_diffusive_backward(j,ntim,repeatInterval)
    end do

    do j=1,nlinks
        write(8, 10) ini_time+real(repeatInterval)*dtini/60.,j, (newY(i,j), i=1, maxval(nx1))
        write(95, 10) ini_time+real(repeatInterval)*dtini/60.,j, (celerity(i,j), i=1, maxval(nx1))
        write(96, 10) ini_time+real(repeatInterval)*dtini/60.,j, (diffusivity(i,j), i=1, maxval(nx1))
        write(991, 10) ini_time+real(repeatInterval)*dtini/60.,j, (newY(i,j)-z(i,j), i=1, maxval(nx1))
        write(993, 10) ini_time+real(repeatInterval)*dtini/60.,j, (velocity(i,j), i=1, maxval(nx1))

        ncomp = nx1(j)
        do i=1,ncomp-1
            courant(i) = celerity(i,j) * dtini / dx(i,j)
        end do
        courant(ncomp:maxval(nx1))=-999
        write(9921, 10) ini_time+real(repeatInterval)*dtini/60.,j, (courant(i), i=1, maxval(nx1))
    end do

    oldY = newY




    !pause 7


    print*, 'Finished timestep', kkk+repeatInterval-1, 'out of ', totalTimeSteps-1, 'min_X_plus',min_X_plus,'min_T_plus',min_T_plus

    end do ! end kkk loop

    close(8)
    close(9)
    close(95)
    close(96)
    close(991)
    close(993)
    close(994)
    close(9921)

10  format(f12.5 ,i6, 2200f14.4)
11  format(f12.5 ,i6, 2200i6)

    call cpu_time( t2 )

    print*, t2-t1, 'sec'

    print*, 'min_X_plus',min_X_plus,'min_T_plus',min_T_plus

end program mesh























function r_interpol_time(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: x(jj), y(jj)
    real, intent(in) :: xt
    real :: yt
    !real(kind=8), intent(out) :: r_interpol_time


    if (xt.le.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        print*, 'jj', jj
        print*, 'x', (x(i), i=1, jj)
        print*, 'y', (y(i), i=1, jj)
        stop
        !if (xt.le. minval(x)) yt=minval(y)
        !if (xt.ge. maxval(x)) yt=maxval(y)
    end if
    r_interpol_time = yt
    ! print*,xt
    return
end function

function r_interpo_nn(x,y,jj,xt)

    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)

    real :: yt

    ! nn means nearest neighbour

    if (xt.le. x(1)) then
        yt=y(1)
    elseif (xt.ge. x(jj)) then
        yt=y(jj)
    else
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    end if
    r_interpo_nn = yt
    return
end function

