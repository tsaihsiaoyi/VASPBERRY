! PROGRAM VASPBERRY Version 1.0 (f77) for VASP
! Written by Hyun-Jung Kim (angpangmokjang@hanmail.net, Infant@kias.re.kr)
!  Korea Institute for Advanced Study (KIAS)
!  Dep. of Phys., Hanyang Univ.
! Copyright 2015. Hyun-Jung Kim All rights reserved.
! Evaluate berry curvature omega(k) for a closed loop on a small patches in k-space
! version 0.1 rough version. not working at all           : 2015. Mar. 17. H.-J. Kim
! version 0.2 error fix for k-loop finding                : 2015. Mar. 18. H.-J. Kim
! version 0.3 ISPIN=1, ISPIN=2 available                  : 2015. Mar. 23. H.-J. Kim
! version 0.4 sign error fix for BZ boundary              : 2015. Mar. 25. H.-J. Kim
! version 0.5 bug fix for defining overlap matrix S(k,k') : 2015. Mar. 31. H.-J. Kim
! version 0.6 including the option for circular dichroism : 2015. Apr. 10. H.-J. Kim
! version 0.7 set coefficients type to be complex*16 and  : 2015. Apr. 12. H.-J. Kim
!             implemented simple routine for grid extending
! version 0.8 routine for velocity expectation value      : 2015. Apr. 30. H.-J. Kim
! version 0.9 Z2 invariant via Fukui's method:not finished: 2016. Feb. 18. H.-J. Kim
! version 0.91 routine for wavefunction plot: -wf n -k k  : 2016. May. 10. H.-J. Kim
! version 0.92 MPI parallelization implemented (MPI_USE)  : 2016. Jun. 27. H.-J. Kim & Y.-K. Kang &
!              (during the CAC workshop & KIAS)                            S.-B. Cho & S.-W. Kim  &
!              -for Berrycurvature & Chern number evaluation               Y.-J. Choi & S.-H. Lee
! version 0.93 MPI parallelization implemented (MPI_USE)  : 2016. Jun. 29. H.-J. Kim
!              -for Z2 invariant evaluation (routines are modified..)
! routine change: routine for getting a determiant -> get_det : 2018, Jun. 28. H.-J. Kim

! last update and bug fixes : 2020. Jun. 10. by H.-J. Kim

!#define MPI_USE
!#undef  MPI_USE

PROGRAM VASPBERRY

    IMPLICIT REAL * 8(a - h, o - z)
    COMPLEX*8, ALLOCATABLE :: coeff(:)
    COMPLEX*16, ALLOCATABLE :: Siju(:, :), Sijd(:, :), Sijt(:, :)
    COMPLEX*16, ALLOCATABLE :: coeff1u(:), coeff1d(:)
    COMPLEX*16, ALLOCATABLE :: coeff2u(:), coeff2d(:)
    COMPLEX*16, ALLOCATABLE :: cener(:)
    REAL*8, ALLOCATABLE :: berrycurv(:), recivec(:, :), wklp(:, :, :)
    REAL*8, ALLOCATABLE :: berrycurv_tot(:)
    REAL*8, ALLOCATABLE :: rnfield(:), rnfield_tot(:)
    REAL*8, ALLOCATABLE :: recivec_tot(:, :)
    REAL*8, ALLOCATABLE :: recilat_tot(:, :)
    REAL*8, ALLOCATABLE :: xrecivec(:, :), xberrycurv(:), wklist(:, :)
    REAL*8, ALLOCATABLE :: wnklist(:, :), xrnfield(:)
    REAL*8, ALLOCATABLE :: recilat(:, :), xrecilat(:, :), occ(:)
    REAL*8, ALLOCATABLE :: selectivity(:), xselectivity(:)
    REAL*16, ALLOCATABLE :: ener(:)
    INTEGER, ALLOCATABLE :: ig(:, :), nplist(:)
    DIMENSION selectivitymax(4), selectivitymin(4)
    DIMENSION a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), a2xa3(3), sumkg(3)
    DIMENSION wk(3), wkk(3, 5), ikk(5), isgg(2, 5), npl(5), itr(5), itrim(5)
    DIMENSION nbmax(3), xb(3), berrymax(4), berrymin(4), ng(3), rs(3)
    COMPLEX * 16 csum1, csum2
    COMPLEX * 16 detS(4), detA, detLOOP
    INTEGER k, n, nkx, nky, nini, nmax, ns, ne, icd, ivel
    CHARACTER * 75 filename, foname, fonameo, fbz, ver_tag, vdirec
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    REAL * 8 rfield, rnnfield, rnnfield_bottom
    REAL * 8 rnnfield_tot, rnnfield_bottom_tot
    REAL*8, ALLOCATABLE :: w_half_klist(:, :)
    INTEGER, ALLOCATABLE :: i_half_klist(:), i_trim_klist(:)
    INTEGER :: myrank, nprocs, ierr
#ifdef MPI_USE
    INCLUDE 'mpif.h'
    INTEGER status(MPI_STATUS_SIZE)
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    IF (myrank == 0) THEN
        WRITE (6, *) "THIS IS ROOT:", myrank
    END IF
    IF (myrank == 0) THEN
        time_1 = MPI_WTIME()
    END IF
#else
    nprocs = 1
    myrank = 0
#endif

    ver_tag = "# VASPBERRY (Ver 1.0), by Hyun-Jung Kim. 2018. Aug. 23."
    pi = 4.*ATAN(1.)

    !default settings
    kperiod = 2
    nkx = 12; nky = 12
    ispinor = 2  ! 2 for soc, 1 for non-soc
    itr = 0; itrim = 0

!!$*  reading general informations
    CALL parse(filename, foname, nkx, nky, ispinor, icd, ixt, fbz, &
               ivel, iz, ihf, nini, nmax, kperiod, it, iskp, ine, ver_tag, &
               iwf, ikwf, ng, rs, imag)
    IF (myrank == 0) CALL creditinfo(ver_tag)

    IF (ixt .NE. 0) CALL extendingBZ(fbz, ixt)
    IF (it .EQ. 1) CALL test
    IF (myrank == 0) WRITE (6, '(A,A)') "# File reading... : ", filename
    CALL inforead(irecl, ispin, nk, nband, ecut, a1, a2, a3, filename)
    IF (ispin .EQ. 2) ispinor = 1
    ALLOCATE (cener(nband), ener(nband), occ(nband))
    IF (myrank == 0) WRITE (6, '(A,I9)') "# TOTAL RECORD LENGTH = ", irecl
    ALLOCATE (wklist(3, nk), nplist(nk))
    IF (iz == 0) THEN
        iz2 = 1 ! enhancing factor for variable size define in subroutines
        ALLOCATE (berrycurv(nk))
        ALLOCATE (berrycurv_tot(nk))
        ALLOCATE (xberrycurv(kperiod * 2 * kperiod * 2 * nk))
    ELSEIF (iz == 1) THEN
        iz2 = 4 ! enhancing factor for variable size define in subroutines
        ALLOCATE (rnfield(nk * iz2), rnfield_tot(nk * iz2))
        ALLOCATE (xrnfield(kperiod * 2 * kperiod * 2 * nk * iz2))
        ALLOCATE (wnklist(3, nk * iz2))
        ALLOCATE (w_half_klist(3, nk), i_half_klist(nk), i_trim_klist(nk))
    END IF
    DO ik = 1, nk
        ne_temp0 = ne
        ne = 0
        irec = 3 + (ik - 1) * (nband + 1)
        READ (10, rec=irec) xnplane, (wk(i), i=1, 3), &
            (cener(nn), occ(nn), nn=1, nband)
        wklist(:, ik) = wk(:)
        nplist(ik) = NINT(xnplane)
        DO n = 1, nband; ne = ne + NINT(occ(n)); END DO
        IF (ik .GT. 1 .AND. ne_temp0 .NE. ne .AND. ine .EQ. 0) THEN
            WRITE (6, *) "error. !!! ne(K) /= ne(K') !!!", ne, ne_temp0, ik; STOP
        END IF
    END DO
    IF (ine .NE. 0) ne = ine ! manually specified ne ; useful for the semimetal
    ! check whether multi or single band calculation is performed
    IF ((nini .EQ. nmax)) THEN
        nini = nmax
    ELSE IF (nmax .EQ. 999999) THEN
        IF (ispin .EQ. 2 .AND. icd .EQ. 0) nmax = ne
        IF (ispin .EQ. 1 .AND. ispinor .EQ. 2 .AND. icd .EQ. 0) nmax = ne
        IF (ispin .EQ. 1 .AND. ispinor .EQ. 1 .AND. icd .EQ. 0) nmax = ne / 2
        IF (ispin .EQ. 2 .AND. icd .EQ. 1) THEN; nini = ne; nmax = nini + 1; END IF
        IF (ispin .EQ. 1 .AND. ispinor .EQ. 1 .AND. icd .EQ. 1) THEN
            nini = ne / 2
            nmax = nini + 1
        END IF
        IF (ispinor .EQ. 2 .AND. icd .EQ. 1) THEN
            nini = ne
            nmax = nini + 1
        END IF
    END IF ! check multi or single ?
    ns = nmax - nini + 1
    IF (myrank == 0) THEN
        WRITE (6, '(A,I6)') "# NELECT     : ", ne * ispin
        IF (ispinor .EQ. 2) THEN
            WRITE (6, '(A,I6,A)') "# ISPIN      : ", ispin, " (LSORBIT =.TRUE.)"
        ELSE
            WRITE (6, '(A,I6,A)') "# ISPIN      : ", ispin, " (LSORBIT =.FALSE.)"
        END IF
        WRITE (6, '(A,F11.4)') "# ENCUT (eV) : ", ecut
        WRITE (6, '(A,I6)') "# NKPOINT    : ", nk
        WRITE (6, '(A,I6,A,I4)') "#  K-GRID    : ", nkx, "   X", nky
        WRITE (6, '(A,I6)') "# NBANDS     : ", nband
        WRITE (6, '(A,3F13.6)') "# LATTVEC A1 : ", (a1(i), i=1, 3)
        WRITE (6, '(A,3F13.6)') "# LATTVEC A2 : ", (a2(i), i=1, 3)
        WRITE (6, '(A,3F13.6)') "# LATTVEC A3 : ", (a3(i), i=1, 3)
    END IF
    CALL recilatt(b1, b2, b3, dSkxky, a1, a2, a3, nkx, nky)
    IF (myrank == 0) THEN
        WRITE (6, '(A,3F13.6)') "# RECIVEC B1 : ", (b1(i), i=1, 3)
        WRITE (6, '(A,3F13.6)') "# RECIVEC B2 : ", (b2(i), i=1, 3)
        WRITE (6, '(A,3F13.6)') "# RECIVEC B3 : ", (b3(i), i=1, 3)
        IF (icd .EQ. 0) WRITE (6, '(A,F13.6)') "#  dk^2 = |dk1xk2| = ", dSkxky
    END IF
    CALL reciproperty(nbmax, npmax, b1, b2, b3, ecut, ispinor)
    IF (myrank == 0) THEN
        WRITE (6, '(A,I6)') "# NPMAX      : ", npmax
        WRITE (6, '(A,I6)') "#  NB1MAX    : ", nbmax(1)
        WRITE (6, '(A,I6)') "#  NB2MAX    : ", nbmax(2)
        WRITE (6, '(A,I6)') "#  NB3MAX    : ", nbmax(3)
    END IF
#ifdef MPI_USE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    ALLOCATE (ig(3, npmax))
    ALLOCATE (coeff(npmax))
    IF (iz == 1) THEN
        ALLOCATE (wklp(3, 5, nk * 4))
        ALLOCATE (recivec(3, 4 * nk), recilat(3, 4 * nk))
        ALLOCATE (recivec_tot(3, 4 * nk), recilat_tot(3, 4 * nk))
        ALLOCATE (xrecivec(3, kperiod * 2 * kperiod * 2 * nk * 4))
        ALLOCATE (xrecilat(3, kperiod * 2 * kperiod * 2 * nk * 4))
    ELSEIF (iz == 0) THEN
        ALLOCATE (wklp(3, 5, nk))
        ALLOCATE (recivec(3, nk), recilat(3, nk))
        ALLOCATE (recivec_tot(3, nk), recilat_tot(3, nk))
        ALLOCATE (xrecivec(3, kperiod * 2 * kperiod * 2 * nk))
        ALLOCATE (xrecilat(3, kperiod * 2 * kperiod * 2 * nk))
    END IF
    nbtot = (2 * nbmax(1) + 2) * (2 * nbmax(2) + 2) * (2 * nbmax(3) + 1)
    ALLOCATE (coeff1u(nbtot), coeff1d(nbtot))
    ALLOCATE (coeff2u(nbtot), coeff2d(nbtot))
    IF (icd .EQ. 1) ALLOCATE (selectivity(nk))
    IF (icd .EQ. 1) ALLOCATE (xselectivity(kperiod * 2 * kperiod * 2 * nk))
    IF (ivel .EQ. 1 .AND. myrank == 0) THEN
        ni = nini
        nj = nini
        CALL vel_expectation(a1, a2, a3, b1, b2, b3, &
                             kperiod, nbmax, npmax, ecut, ispinor, ispin, nband, ne, nk, wklist, &
                             nplist, ni, nj, irecl, filename, foname)
    END IF

    IF (iwf .GE. 1 .AND. myrank == 0) THEN ! plotting wavefunction
        ilp = 1
        ikk(ilp) = ikwf
        iband = iwf
        npl(ilp) = nplist(ikk(ilp))
        wk(:) = wklist(:, ikk(ilp))
        nx = ng(1)
        ny = ng(2)
        nz = ng(3)
        IF (nx * ny * nz .EQ. 0) THEN
            nx = nbmax(1)
            ny = nbmax(2)
            nz = nbmax(3)
        END IF
        WRITE (6, '(A)') " "
        WRITE (6, '(A)') "# Evaluating wavefunction..."
        CALL wf_plot(a1, a2, a3, b1, b2, b3, wk, &
                     nbmax, npmax, ecut, ispinor, ispin, nband, nk, &
                     npl, nmax, irecl, filename, foname, nx, ny, nz, iband, ikk, 1, rs, imag)

    END IF !iwf end

    IF (iz .EQ. 1) THEN  !get Z2 invariant using Fukui's method
        del = 1E-6 ! criterion for k-point find
        CALL set_BZ_fukui(nnk, wnklist, w_half_klist, i_half_klist, &
                          i_trim_klist, nhf, ispin, wklist, del, nk, iz2)

        DO isp = 1, ispin !ispin start
            IF (myrank == 0) THEN
                WRITE (6, *) " "
                IF (isp .EQ. 1 .AND. ispin .EQ. 2) WRITE (6, '(A,I2)') "# SPIN : UP"
                IF (isp .EQ. 2 .AND. ispin .EQ. 2) WRITE (6, '(A,I2)') "# SPIN : DN"
                WRITE (6, *) " "
            END IF
#ifdef MPI_USE
            nk_rank = nnk / nprocs
            nk_rank_mod = MOD(nnk, nprocs)

            IF (myrank < nk_rank_mod) THEN
                ink = myrank * (nk_rank + 1) + 1
                ifk = (myrank + 1) * (nk_rank + 1)
            ELSE
                ink = myrank * nk_rank + nk_rank_mod + 1
                ifk = (myrank + 1) * nk_rank + nk_rank_mod
            END IF

            recivec = 0.
            recilat = 0.
            recivec_tot = 0.
            recilat_tot = 0.
            rnnfield = 0.
            rnnfield_bottom = 0.
            rnnfield_tot = 0.
            rnnfield_bottom_tot = 0.
            rnfield_tot = 0.
            rnfield = 0.
#else
            ink = 1; ifk = nnk
#endif
            DO ik = ink, ifk !ik - loop
                CALL get_nfield(rnfield(ik), wklp, isp, ik, nk, nband, nkx, nky, &
                                wnklist, ihf, nplist, nnk, iz, ns, kperiod, &
                                w_half_klist, i_half_klist, i_trim_klist, nhf, &
                                nbmax, npmax, ecut, ispinor, ispin, ne, irecl, &
                                nmax, nini, b1, b2, b3, iz2)
                IF ((wklp(2, 1, ik) + wklp(2, 3, ik)) / 2. .GT. -del .AND. &
                    (wklp(2, 1, ik) + wklp(2, 3, ik)) / 2. .LT. 0.5 + del) THEN
                    rnnfield = rnnfield + rnfield(ik)
                ELSE
                    rnnfield_bottom = rnnfield_bottom + rnfield(ik)
                END IF

!       write(6,'(A)')"#__________________________
!    &________________________________ "
!       write(6,'(A)')"# Berry Curvature F (A^-2) :
!    & -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
!       write(6,'(A)')"# N-field strength       :
!    & Sum_s{Im[log(det(S(K_s,k_s+1)))]}/2pi - F/2pi "

                DO j = 1, 3
                    recivec(j, ik) = (wklp(1, 1, ik) + wklp(1, 3, ik)) * b1(j) * 0.5 + &
                                     (wklp(2, 1, ik) + wklp(2, 3, ik)) * b2(j) * 0.5 + &
                                     (wklp(3, 1, ik) + wklp(3, 3, ik)) * b3(j) * 0.5
                    recilat(j, ik) = (wklp(j, 1, ik) + wklp(j, 3, ik)) * 0.5
                END DO

!       write(6,'(A)')"#      kx        ky        kz(A^-1)
!    &    n-field strength      ,IK"
!       write(6,'(A,3F11.6,A,F8.3,I4,A)')'#',(recivec(i,ik),i=1,3),
!    &                             "      ",rnfield(ik),ik,"th-K"

            END DO   ! ik    end

#ifdef MPI_USE
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(rnnfield, rnnfield_tot, 1, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(rnnfield_bottom, rnnfield_bottom_tot, 1, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(rnfield, rnfield_tot, iz2 * nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(recivec, recivec_tot, 3 * iz2 * nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(recilat, recilat_tot, 3 * iz2 * nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)

            rnnfield = rnnfield_tot
            rnnfield_bottom = rnnfield_bottom_tot
            rnfield = rnfield_tot
            recivec = recivec_tot
            recilat = recilat_tot
#endif

            IF (myrank == 0) THEN
                WRITE (6, *) " "
                WRITE (6, '(A)') "# DONE!"

                IF (nini .EQ. nmax) THEN
                    WRITE (6, '(A,I4)') "# Z2 invariant for the BAND : ", nini
                ELSE
                    WRITE (6, '(A,I4,A,I4)') "# Z2 invariant for the BANDS : ", nini &
                        , " - ", nmax
                END IF
                WRITE (6, '(A,I2)') "# Z2 Invariant =    ", &
                    MOD((NINT(rnnfield)), 2)
                WRITE (6, '(A,I2)') "# Z2 Invariant(bottom) =    ", &
                    MOD((NINT(rnnfield_bottom)), 2)

                CALL get_ext_variable(xrecivec, xrecilat, xrnfield, kext, &
                                      recivec, recilat, rnfield, nnk, kperiod, nk, iz2, &
                                      b1, b2, b3) ! extending over exteded-BZ

                CALL get_sorted_xrvari(xrecivec, xrecilat, xrnfield, &
                                       kext, kperiod, nk, iz2) ! sorting
                CALL write_result(isp, ispin, ispinor, fonameo, foname, filename, &
                                  irecl, ecut, nk, nkx, nky, nband, b1, b2, b3, kperiod, &
                                  dSkxky, nini, nmax, xrecivec, xrecilat, kext, &
                                  xrnfield, rnnfield, rnnfield_bottom, 0, 0, &
                                  iz2, icd, iz, ivel, nprocs)
            END IF !myrank ==0
        END DO !ispin end

        GOTO 9999
    END IF !if Z2 end

    DO isp = 1, ispin   ! ispin start
        IF (icd + ivel + iz + iwf .EQ. 0) THEN
            chernnumber = 0.
            IF (myrank == 0) THEN
                WRITE (6, *) " "
                IF (isp .EQ. 1 .AND. ispin .EQ. 2) WRITE (6, '(A,I2)') "# SPIN : UP"
                IF (isp .EQ. 2 .AND. ispin .EQ. 2) WRITE (6, '(A,I2)') "# SPIN : DN"
                WRITE (6, *) " "
            END IF

#ifdef MPI_USE
            nk_rank = nk / nprocs
            nk_rank_mod = MOD(nk, nprocs)
            IF (myrank < nk_rank_mod) THEN
                ink = myrank * (nk_rank + 1) + 1
                ifk = (myrank + 1) * (nk_rank + 1)
            ELSE
                ink = myrank * nk_rank + nk_rank_mod + 1
                ifk = (myrank + 1) * nk_rank + nk_rank_mod
            END IF
            IF (myrank == 0) THEN
                time_2 = MPI_Wtime()
            END IF
            berrycurv = 0.
            berrycurv_tot = 0.
            recivec = 0.
            recilat = 0.
            recivec_tot = 0.
            recilat_tot = 0.
#else
            ink = 1; ifk = nk
#endif
            DO ik = ink, ifk     !ik loop.
                CALL klpfind(wkk, isp, ik, nk, nband, nkx, nky, 0, 0, iz, iz2) !k-loop(4pts) of ik-kpoint (K)
                WRITE (6, *) " "
                WRITE (6, 490) ik, (wkk(i, 1), i=1, 3)
490             FORMAT("#* Closed loop for KPOINT", I5, " : (", 3F9.5, "  )")
                CALL kindxfind(ikk, isgg, npl, itr, itrim, wklist, nplist, wkk, nk, &
                               nband, &
                               iz, w_half_klist, i_half_klist, 0, 0, iz2) !find K-loop for ik
                WRITE (6, 500) (ikk(i), i=1, 5)
500             FORMAT('# K1>K2>K3>K4>K5=K1 (k-index):', &
                       I5, '  >', I5, '  >', I5, '  >', I5, '  >', I5)
                DO j = 1, 5
                    WRITE (6, 510) j, (wkk(i, j), i=1, 3), j, (isgg(i, j), i=1, 2)
510                 FORMAT('# K', I1, ' = (', 3F9.5, '  )', &
                           ', G', I1, ' = (n1,n2) = (', I1, ',', I1, ')')
                END DO
                wklp(:, :, ik) = wkk(:, :)  ! closed loop (C) for each K : RECIVEC

!!$* get overlap matrix S_ij(k,k+1) over the C and PI_s [det S(k_s,k_s+1)], s=1,4
                ALLOCATE (Siju(ns, ns), Sijd(ns, ns), Sijt(ns, ns))
                DO ilp = 1, 4  ! berry curvature loop for ik
                    Siju = (0., 0.)
                    Sijd = (0., 0.)
                    Sijt = (0., 0.)
                    ! construct overlap matrix S_ij(ikk1,ikk2)
!         write(6,570)ilp,ilp+1,ilp,ilp+1
! 570     format('#  Constructing overlap matrix, S_ij(K',I1,',K',I1,
!    &          ') = <u_i,K',I1,'(r)|u_j,K',I1,'(r)>')
                    DO ni = nini, nmax   ! calculate upto valence band maximum
                        ncnt = 0; coeff1u = (0., 0.); coeff1d = (0., 0.); coeff = (0., 0.)
                        READ (10, rec=(3 + (ikk(ilp) - 1) * (nband + 1) + &
                                       nk * (nband + 1) * (isp - 1) + ni)) (coeff(i), i=1, npl(ilp))
                        DO ig3 = 0, 2 * nbmax(3); ig3p = ig3
                            IF (ig3 .GT. nbmax(3)) ig3p = ig3 - 2 * nbmax(3) - 1
                            DO ig2 = 0, 2 * nbmax(2); ig2p = ig2
                                IF (ig2 .GT. nbmax(2)) ig2p = ig2 - 2 * nbmax(2) - 1
                                DO ig1 = 0, 2 * nbmax(1); ig1p = ig1
                                    IF (ig1 .GT. nbmax(1)) ig1p = ig1 - 2 * nbmax(1) - 1
                                    CALL get_etot(etot, ilp, ni, b1, b2, b3, isgg, ig1p, ig2p, ig3p, &
                                                  ig1, ig2, ig3, wkk, itrim, itr)

                                    IF (etot .LT. ecut) THEN; ncnt = ncnt + 1
                                        CALL get_incnt(incnt, ilp, ni, ig1p, ig2p, ig3p, nbmax, &
                                                       itr, itrim, isgg, wkk)
                                        IF (ispinor .EQ. 2) THEN
                                            coeff1d(incnt) = coeff(ncnt + npl(ilp) / 2)
                                        END IF
                                        coeff1u(incnt) = coeff(ncnt)
                                    END IF
                                END DO  !loop for ig1
                            END DO   !loop for ig2
                        END DO    !loop for ig3
                        IF (ispinor * ncnt .NE. npl(ilp)) THEN
                            WRITE (0, *) '*ni error - computed NPL=', 2 * ncnt, &
                                ' != input nplane=', npl(ilp); STOP
                        END IF

                        DO nj = nini, nmax
                            ncnt = 0; coeff2u = (0., 0.); coeff2d = (0., 0.); coeff = (0., 0.)
                            READ (10, rec=(3 + (ikk(ilp + 1) - 1) * (nband + 1) + nk * (nband + 1) * (isp - 1) &
                                           + nj)) (coeff(i), i=1, npl(ilp + 1))
                            DO ig3 = 0, 2 * nbmax(3); ig3p = ig3
                                IF (ig3 .GT. nbmax(3)) ig3p = ig3 - 2 * nbmax(3) - 1
                                DO ig2 = 0, 2 * nbmax(2); ig2p = ig2
                                    IF (ig2 .GT. nbmax(2)) ig2p = ig2 - 2 * nbmax(2) - 1
                                    DO ig1 = 0, 2 * nbmax(1); ig1p = ig1
                                        IF (ig1 .GT. nbmax(1)) ig1p = ig1 - 2 * nbmax(1) - 1
                                        CALL get_etot(etot, ilp + 1, nj, b1, b2, b3, isgg, ig1p, ig2p, ig3p, &
                                                      ig1, ig2, ig3, wkk, itrim, itr)
                                        IF (etot .LT. ecut) THEN; ncnt = ncnt + 1
                                            CALL get_incnt(incnt, ilp + 1, nj, ig1p, ig2p, ig3p, nbmax, &
                                                           itr, itrim, isgg, wkk)
                                            IF (ispinor .EQ. 2) THEN
                                                coeff2d(incnt) = coeff(ncnt + npl(ilp + 1) / 2)
                                            END IF
                                            coeff2u(incnt) = coeff(ncnt)
                                        END IF
                                    END DO   !loop for ig1
                                END DO    !loop for ig2
                            END DO     !loop for ig3
                            IF (ispinor * ncnt .NE. npl(ilp + 1)) THEN
                                WRITE (0, *) '*nj error - computed NPL=', 2 * ncnt, &
                                    ' != input nplane=', npl(ilp + 1); STOP
                            END IF
                            IF (ispinor .EQ. 2) THEN
                                Siju(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1u, coeff2u)
                                Sijd(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1d, coeff2d)
                                Sijt(ni - nmax + ns, nj - nmax + ns) = Siju(ni - nmax + ns, nj - nmax + ns) + &
                                                                       Sijd(ni - nmax + ns, nj - nmax + ns)
                            ELSE IF (ispinor .EQ. 1) THEN
                                Sijt(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1u, coeff2u)
                            END IF
                        END DO  ! loop for nj
                    END DO   ! loop for ni

                    IF (nini .EQ. nmax) THEN   ! get determinant : det(S)
                        detS(ilp) = Sijt(1, 1)
                    ELSE
!           call getdetA(detA, Sijt,ns)
                        CALL get_det(detA, Sijt, ns)
                        detS(ilp) = detA; detA = (0., 0.)
                    END IF
                END DO    ! loop for ilp

                detLOOP = detS(1) * detS(2) * detS(3) * detS(4)
                WRITE (6, '(A)') "# "
                WRITE (6, 600) detLOOP
600             FORMAT('#===>PI_S[det(S(K_s,K_s+1))] =', F16.8, '  +', F16.8, ' i')

                berrycurv(ik) = -1.*AIMAG(LOG(detLOOP)) / dSkxky

                WRITE (6, '(A)') "#__________________________&
       &     ________________________________ "
                WRITE (6, '(A)') "# Berry Curvature (A^-2) :&
       &      -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
                DO j = 1, 3
                    recivec(j, ik) = (wklp(1, 1, ik) + wklp(1, 3, ik)) * b1(j) * 0.5 + &
                                     (wklp(2, 1, ik) + wklp(2, 3, ik)) * b2(j) * 0.5 + &
                                     (wklp(3, 1, ik) + wklp(3, 3, ik)) * b3(j) * 0.5
                    recilat(j, ik) = (wklp(j, 1, ik) + wklp(j, 3, ik)) * 0.5
                END DO
                WRITE (6, '(A)') "#      kx        ky        kz(A^-1)&
       &        Berry Curvature (A^-2)"
                WRITE (6, '(A,3F11.6,A,F16.6)') '#', (recivec(i, ik), i=1, 3), "     ", &
                    berrycurv(ik)

#ifdef MPI_USE
                chernnumber = chernnumber + berrycurv(ik) * dSkxky / (2.*pi)
#else
                chernnumber = chernnumber + berrycurv(ik) * dSkxky / (2.*pi)
                IF (ik .EQ. 1) THEN
                    berrymax(4) = berrycurv(ik)
                    berrymax(1) = recilat(1, ik)
                    berrymax(2) = recilat(2, ik)
                    berrymax(3) = recilat(3, ik)
                    berrymin(4) = berrycurv(ik)
                    berrymin(1) = recilat(1, ik)
                    berrymin(2) = recilat(2, ik)
                    berrymin(3) = recilat(3, ik)
                ELSE IF (ik .GE. 2 .AND. berrycurv(ik) .GE. berrymax(4)) THEN
                    berrymax(4) = berrycurv(ik)
                    berrymax(1) = recilat(1, ik)
                    berrymax(2) = recilat(2, ik)
                    berrymax(3) = recilat(3, ik)
                ELSE IF (ik .GE. 2 .AND. berrycurv(ik) .LE. berrymin(4)) THEN
                    berrymin(4) = berrycurv(ik)
                    berrymin(1) = recilat(1, ik)
                    berrymin(2) = recilat(2, ik)
                    berrymin(3) = recilat(3, ik)
                END IF
#endif

                DEALLOCATE (Sijt)
                DEALLOCATE (Siju)
                DEALLOCATE (Sijd)
            END DO         ! ik loop over

#ifdef MPI_USE
            IF (myrank == 0) THEN
                time_3 = MPI_Wtime()
                WRITE (6, *) " "
                WRITE (6, '(A)') "# DONE!"
            END IF
#else
            WRITE (6, *) " "
            WRITE (6, '(A)') "# DONE!"
#endif

!!$*  SUMMARIZTION and output the results..
            IF (myrank == 0) THEN
                WRITE (6, '(A)') "# Chern Number is sum of Berry Curvature on 1BZ"
                IF (nini .EQ. nmax) THEN
                    WRITE (6, '(A,I4)') "# Chern Number for the BAND : ", nini
                ELSE
                    WRITE (6, '(A,I4,A,I4)') "# Chern Number for the BANDS : ", nini &
                        , " - ", nmax
                END IF
            END IF

#ifdef MPI_USE
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(chernnumber, chernnumber_tot, 1, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(berrycurv, berrycurv_tot, nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(recivec, recivec_tot, 3 * nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(recilat, recilat_tot, 3 * nk, MPI_REAL8, &
                            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            chernnumber = chernnumber_tot
            berrycurv(:) = berrycurv_tot(:)
            recivec(:, :) = recivec_tot(:, :)
            recilat(:, :) = recilat_tot(:, :)
#endif
            IF (myrank == 0) THEN
                WRITE (6, '(A,F16.6)') "# Chern Number =    ", chernnumber
            END IF
!KKKKK NOTE : START-sorting

            IF (myrank == 0) THEN
                CALL get_ext_variable(xrecivec, xrecilat, xberrycurv, kext, &
                                      recivec, recilat, berrycurv, nk, kperiod, nk, iz2, &
                                      b1, b2, b3) ! extending over exteded-BZ

                CALL get_sorted_xrvari(xrecivec, xrecilat, xberrycurv, &
                                       kext, kperiod, nk, iz2) ! sorting
#ifdef MPI_USE
                CALL write_result(isp, ispin, ispinor, fonameo, foname, filename, &
                                  irecl, ecut, nk, nkx, nky, nband, b1, b2, b3, kperiod, &
                                  dSkxky, nini, nmax, xrecivec, xrecilat, kext, &
                                  xberrycurv, chernnumber, 0., 0, 0, &
                                  iz2, icd, iz, ivel, nprocs)
#else
                CALL write_result(isp, ispin, ispinor, fonameo, foname, filename, &
                                  irecl, ecut, nk, nkx, nky, nband, b1, b2, b3, kperiod, &
                                  dSkxky, nini, nmax, xrecivec, xrecilat, kext, &
                                  xberrycurv, chernnumber, 0., berrymax, berrymin, &
                                  iz2, icd, iz, ivel, nprocs)
#endif
            END IF

!!! ######### LOOP for OPTICAL SELECTIVITY ###########################################
        ELSEIF (icd .EQ. 1 .AND. myrank == 0) THEN   ! calculate optical selectivity
            CALL optical_selectivity(selectivity, &
                                     selectivitymax, selectivitymin, &
                                     b1, b2, b3, wklist, isp, &
                                     nband, ecut, ispinor, nplist, nbmax, npmax, nk, nini, nmax)
            DO ik = 1, nk
                DO j = 1, 3
                    recivec(j, ik) = wklist(1, ik) * b1(j) + &
                                     wklist(2, ik) * b2(j) + &
                                     wklist(3, ik) * b3(j)
                    recilat(j, ik) = wklist(j, ik)
                END DO
            END DO

            CALL get_ext_variable(xrecivec, xrecilat, xselectivity, kext, &
                                  recivec, recilat, selectivity, nk, kperiod, nk, iz2, &
                                  b1, b2, b3) ! extending over exteded-BZ

            CALL get_sorted_xrvari(xrecivec, xrecilat, xselectivity, &
                                   kext, kperiod, nk, iz2) ! sorting
            CALL write_result(isp, ispin, ispinor, fonameo, foname, filename, &
                              irecl, ecut, nk, nkx, nky, nband, b1, b2, b3, kperiod, &
                              dSkxky, nini, nmax, xrecivec, xrecilat, kext, &
                              xselectivity, 0, 0., selectivitymax, selectivitymin, &
                              iz2, icd, iz, ivel, nprocs)

        END IF !icd over
!!! ######### LOOP END for OPTICAL SELECTIVITY ###########################################

    END DO !ispin loop over

    IF (iskp .EQ. 1 .AND. myrank == 0) CALL write_special_kpoint(b1, b2, b3)
#ifdef MPI_USE
    IF (myrank == 0) THEN
#endif
        WRITE (6, '(A)') "# DONE! "
        DO isp = 1, ispin
            IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
                WRITE (6, '(A,A,A)') "#  Results are summarized in ", TRIM(foname), &
                    ".UP.dat for spin-1"
            ELSE IF (isp .EQ. 2 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
                WRITE (6, '(A,A,A)') "#  Results are summarized in ", &
                    TRIM(foname), ".DN.dat for spin-2"
            ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 2) THEN
                WRITE (6, '(A,A,A)') "#  Results are summarized in ", &
                    TRIM(foname), ".dat"
            ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 1) THEN
                WRITE (6, '(A,A,A)') "#  Results are summarized in ", &
                    TRIM(foname), ".dat"
            END IF
        END DO
#ifdef MPI_USE
    END IF
#endif
9999 IF (myrank == 0) WRITE (6, *) "end of program"

#ifdef MPI_USE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    DEALLOCATE (ig)
    DEALLOCATE (coeff)
    DEALLOCATE (coeff1u)
    DEALLOCATE (coeff1d)
    DEALLOCATE (coeff2u)
    DEALLOCATE (coeff2d)
    DEALLOCATE (recilat)
    DEALLOCATE (recivec)
    DEALLOCATE (xrecilat)
    DEALLOCATE (xrecivec)
    DEALLOCATE (wklp)
    DEALLOCATE (wklist)
    DEALLOCATE (cener)
    DEALLOCATE (ener)
    IF (iz == 0) THEN
        DEALLOCATE (berrycurv)
        DEALLOCATE (berrycurv_tot)
        DEALLOCATE (xberrycurv)
    ELSEIF (iz == 1) THEN
        DEALLOCATE (rnfield)
        DEALLOCATE (rnfield_tot)
        DEALLOCATE (wnklist)
        DEALLOCATE (xrnfield)
        DEALLOCATE (i_trim_klist)
        DEALLOCATE (i_half_klist)
        DEALLOCATE (w_half_klist)
    END IF
    DEALLOCATE (occ)
    DEALLOCATE (nplist)
    IF (icd .EQ. 1) DEALLOCATE (selectivity)
    IF (icd .EQ. 1) DEALLOCATE (xselectivity)
#ifdef MPI_USE
    IF (myrank == 0 .AND. iz + ivel + icd + iwf .EQ. 0) THEN
        time_4 = MPI_Wtime()
        WRITE (6, *) "Data reading      : ", time_2 - time_1
        WRITE (6, *) "Parallel sequence : ", time_3 - time_1
        WRITE (6, *) "End sequence      : ", time_4 - time_3
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    CALL MPI_FINALIZE(ierr)
#endif
END PROGRAM

!!$*  subroutine for writing results
SUBROUTINE write_result(isp, ispin, ispinor, fonameo, foname, filename, &
                        irecl, ecut, nk, nkx, nky, nband, b1, b2, b3, kperiod, &
                        dSkxky, nini, nmax, xrecivec, xrecilat, kext, &
                        xrvari, rvari, rvari2, rvari3, rvari4, &
                        iz2, icd, iz, ivel, nprocs)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION b1(3), b2(3), b3(3)
    REAL * 8 xrecivec(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xrecilat(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xrvari(kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 rvari, rvari2
    REAL * 8 xb(3), xtemp, rvari3(4), rvari4(4)
    CHARACTER * 75 filename, foname, fonameo

    IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
        WRITE (fonameo, '(A,A)') TRIM(foname), '.UP.dat'
    ELSE IF (isp .EQ. 2 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
        WRITE (fonameo, '(A,A)') TRIM(foname), '.DN.dat'
    ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 2) THEN
        WRITE (fonameo, '(A,A)') TRIM(foname), '.dat'
    ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 1) THEN
        WRITE (fonameo, '(A,A)') TRIM(foname), '.dat'
    END IF
    OPEN (32, file=fonameo, status='unknown')
    WRITE (32, '(A,I4,A)') "# Job running on ", nprocs, " total cores"
    WRITE (32, '(A,A)') "# File reading...  : ", filename
    WRITE (32, '(A,I9)') "# TOTAL RECORD LENGTH = ", irecl
    IF (ispinor .EQ. 2) THEN
        WRITE (32, '(A,I6,A)') "# ISPIN            : ", ispin, &
            " (LSORBIT = .TRUE.)"
    ELSE
        WRITE (32, '(A,I6,A)') "# ISPIN            : ", ispin, &
            " (LSORBIT = .FALSE.)"
    END IF
    WRITE (32, '(A,F11.4)') "# ENCUT (eV)       : ", ecut
    WRITE (32, '(A,I6)') "# NKPOINT          : ", nk
    WRITE (32, '(A,I6,A,I4)') "#  K-GRID          : ", nkx, "   X", nky
    WRITE (32, '(A,I6)') "# NBANDS           : ", nband
    WRITE (32, '(A,3F13.6)') "# RECIVEC B1 (A^-1): ", (b1(i), i=1, 3)
    WRITE (32, '(A,3F13.6)') "# RECIVEC B2       : ", (b2(i), i=1, 3)
    WRITE (32, '(A,3F13.6)') "# RECIVEC B3       : ", (b3(i), i=1, 3)
    WRITE (32, '(A,F13.6)') "#  dk^2 = |dk1xk2| = ", dSkxky
    WRITE (32, *) " "

    IF (nini .EQ. nmax) THEN
        IF (iz == 1) THEN
            WRITE (32, '(A,I4)') "# Z2 invariant for the BAND : ", nini
        ELSEIF (iz + ivel + icd .EQ. 0) THEN
            WRITE (32, '(A,I4)') "# Chern Number for the BAND : ", nmax
        END IF
    ELSE
        IF (iz == 1) THEN
            WRITE (32, '(A,I4,A,I4)') "# Z2 invariant for the BANDS : ", nini, &
                " - ", nmax
        ELSEIF (iz + ivel + icd .EQ. 0) THEN
            WRITE (32, '(A,I4,A,I4)') "# Chern Number for the BANDS : ", nini, &
                "    -  ", nmax
        END IF
    END IF

    IF (iz == 1) THEN !Z2 INVARIANT
        WRITE (32, '(A,I2)') "# Z2 Invariant (top) =    ", &
            MOD((NINT(rvari)), 2)
        WRITE (32, '(A,I2)') "# Z2 Invariant (bottom) =    ", &
            MOD((NINT(rvari2)), 2)
        WRITE (32, '(A)') "# Berry Curvature F (A^-2) :&
&       -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
        WRITE (32, '(A)') "# N-field strength       :&
&       Sum_s{Im[log(det(S(K_s,k_s+1)))]}/2pi - F/2pi "
        WRITE (32, '(A)') "# (cart) kx        ky        kz(A^-1)&
&         n-field strength      ,   (recip)kx        ky        kz"

    ELSEIF (icd == 1) THEN !OPTICAL SELECTIVITY
        WRITE (32, '(A,I4,A,I4)') "# OPTICAL SELECTIVITY BETWEEN BANDS: ", &
            nini, "    -  ", nmax
        WRITE (32, '(A)') "# n(k,w_cv)= |P(k,s,cv,+)|^2 - |P(k,s,cv,-)|^2"
        WRITE (32, '(A)') "#            ---------------------------------"
        WRITE (32, '(A)') "#            |P(k,s,cv,+)|^2 + |P(k,s,cv,-)|^2"
        WRITE (32, '(A)') "#  The TRANSITION MATRIX ELEMENT P ="
        WRITE (32, '(A)') "#   P(k,s,cv,+ or -) = 1/sqrt(2)[p_x(k,cv,s) +&
&     (or -) i*p_y(k,cv,s)]"
        WRITE (32, '(A)') "#  THE INTERBAND TRANSITION MATRIX p_x,y ="
        WRITE (32, '(A)') "#   p_x,y(k,cv,s)=<psi(k,c,s)|-i*hbar*1/dx(y)|&
&     psi(k,v,s>"
        WRITE (32, '(A,4F16.6)') "# MAXVAL of SELECTIVITY at kx,ky,kz &
&     (in reci)= ", (rvari3(i), i=1, 4)
        WRITE (32, '(A,4F16.6)') "# MINVAL of SELECTIVITY at kx,ky,kz &
&     (in reci)= ", (rvari4(i), i=1, 4)
        WRITE (32, '(A)') "# (cart) kx        ky        kz(A^-1)&
&        selectivity(n(k)),        (recip)kx        ky        kz"

    ELSE !BERRYCURVATURE
        WRITE (32, '(A)') "# Chern Number is sum of &
&     Berry Curvature over 1BZ"
        WRITE (32, '(A,F16.6)') "# Chern Number =   ", rvari
        WRITE (32, '(A,4F16.6)') "# MAXVAL of BERRYCURV at kx,ky,kz &
&     (in reci)= ", (rvari3(i), i=1, 4)
        WRITE (32, '(A,4F16.6)') "# MINVAL of BERRYCURV at kx,ky,kz &
&     (in reci)= ", (rvari4(i), i=1, 4)
        WRITE (32, '(A)') "# Berry Curvature (A^-2) :&
&     -Im[logPI_S(det(S(K_s,K_s+1)))]/dk^2"
        WRITE (32, '(A)') "# (cart) kx        ky        kz(A^-1)&
&        Berry Curvature (A^-2),   (recip)kx        ky        kz"
    END IF

    DO ik = 1, kext
        WRITE (32, '(3F11.6,A,F11.6,A,3F11.6)') (xrecivec(i, ik), i=1, 3), &
            "     ", xrvari(ik), "            ", &
            (xrecilat(i, ik), i=1, 3)
    END DO

    CLOSE (32)
    RETURN
END SUBROUTINE write_result

!!$*  subroutine for sorting
SUBROUTINE get_sorted_xrvari(xrecivec, xrecilat, xrvari, &
                             kext, kperiod, nk, iz2)
    IMPLICIT REAL * 8(a - h, o - z)
    REAL * 8 xrecivec(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xrecilat(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xrvari(kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xb(3), xtemp

    DO k = kext - 1, 1, -1  ! sorting kx
        DO j = 1, k
            IF (xrecivec(1, j + 1) .GT. xrecivec(1, j)) THEN
                xb(:) = xrecivec(:, j)
                xrecivec(:, j) = xrecivec(:, j + 1)
                xrecivec(:, j + 1) = xb(:)
                xb(:) = xrecilat(:, j)
                xrecilat(:, j) = xrecilat(:, j + 1)
                xrecilat(:, j + 1) = xb(:)
                xtemp = xrvari(j)
                xrvari(j) = xrvari(j + 1)
                xrvari(j + 1) = xtemp
            END IF
        END DO
    END DO
    DO k = kext - 1, 1, -1  ! sorting ky
        DO j = 1, k
            IF (xrecivec(1, j + 1) .EQ. xrecivec(1, j)) THEN
                IF (xrecivec(2, j + 1) .GT. xrecivec(2, j)) THEN
                    xb(:) = xrecivec(:, j)
                    xrecivec(:, j) = xrecivec(:, j + 1)
                    xrecivec(:, j + 1) = xb(:)
                    xb(:) = xrecilat(:, j)
                    xrecilat(:, j) = xrecilat(:, j + 1)
                    xrecilat(:, j + 1) = xb(:)
                    xtemp = xrvari(j)
                    xrvari(j) = xrvari(j + 1)
                    xrvari(j + 1) = xtemp
                END IF
            END IF
        END DO
    END DO
    RETURN
END SUBROUTINE get_sorted_xrvari

!!$*  subroutine for extending data_set distribution over extended BZ: 2D array
SUBROUTINE get_ext_variable(xrecivec, xrecilat, xrvari, kext, &
                            recivec, recilat, rvari, nnk, kperiod, nk, iz2, &
                            b1, b2, b3) !iz2=4 for z2, other 1
    IMPLICIT REAL * 8(a - h, o - z)
    REAL * 8 recivec(3, nk * iz2), recilat(3, nk * iz2)
    REAL * 8 xrecivec(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 xrecilat(3, kperiod * 2 * kperiod * 2 * nk * iz2)
    REAL * 8 rvari(nk * iz2), xrvari(kperiod * 2 * kperiod * 2 * nk * iz2)
    DIMENSION b1(3), b2(3), b3(3)

    kk = 0   ! extend variable distribution over extended BZ
    DO ib2 = -1 * (kperiod - 1) + 1, kperiod
        DO ib1 = -1 * (kperiod - 1) + 1, kperiod
            DO ik = 1, nnk
                kk = kk + 1
                xrecivec(1, kk) = recivec(1, ik) + (ib1 - 1) * b1(1) + (ib2 - 1) * b2(1)
                xrecivec(2, kk) = recivec(2, ik) + (ib1 - 1) * b1(2) + (ib2 - 1) * b2(2)
                xrecivec(3, kk) = recivec(3, ik) + (ib1 - 1) * b1(3) + (ib2 - 1) * b2(3)
                xrecilat(1, kk) = recilat(1, ik) + (ib1 - 1)
                xrecilat(2, kk) = recilat(2, ik) + (ib2 - 1)
                xrecilat(3, kk) = recilat(3, ik)
                xrvari(kk) = rvari(ik)
                kext = kk
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE

!!$*  subroutine for getting planewave coefficient for given k & n
SUBROUTINE wf_plot(a1, a2, a3, b1, b2, b3, wk, &
                   nbmax, npmax, ecut, ispinor, ispin, nband, nk, &
                   npl, nmax, irecl, filename, foname, nx, ny, nz, iband, ikk, ilp, rs, imag)

    IMPLICIT REAL * 8(a - h, o - z)
    COMPLEX * 8 coeff(npmax), coeffi(npmax)
    COMPLEX * 16 csum(ispinor, nx * ny * nz)
    REAL * 8 wklist(3, nk), w_half_klist(3, nk), pi, pi2
    REAL * 8 recivec(3), recilat(3, nk * 4)
    DIMENSION a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), a2xa3(3), sumkg(3)
    DIMENSION wk(3), wkg(3), rr(3), rs(3)
    DIMENSION wkgr(npmax)
    INTEGER nplist(nk)
    INTEGER ikk(5), isgg(2, 5), npl(5), itrim(5), ig(3, npmax)
    INTEGER nbmax(3)
    INTEGER nk, nband, npmax, kperiod, ispin, irecl
    CHARACTER * 75 filename, foname, fonameo, fonameoi
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    pi = 4.*ATAN(1.)
    pi2 = pi * 2.

    itrim = 0
    CALL vcross(a2xa3, a2, a3)
    vol = DOT_PRODUCT(a1, a2xa3)
    ngrid = nx * ny * nz
    nline = INT(ngrid / 5.)
    nresi = MOD(ngrid, 5)

    DO isp = 1, ispin
        WRITE (6, '(A,I1)') " ##Processing for spin-", isp
        csum = (0., 0.); wkgr = 0.
        WRITE (fonameo, '(A,A,I1)') TRIM(foname), "-SPIN", isp
        CALL get_coeff(coeff, isgg, 1, iband, ikk, npl, nband, nk, isp, &
                       npmax, itrim)
        CALL plindx(ig, ncnt, ispinor, wk, b1, b2, b3, nbmax, npl(1), ecut, npmax)

        OPEN (unit=18, file=fonameo, status='unknown')
        CALL write_CHGCAR_head(18, fonameo, ikk, isp, iband, a1, a2, a3, &
                               nx, ny, nz, rs)

        IF (imag .EQ. 1) THEN
            WRITE (fonameoi, '(A,A,I1)') TRIM(foname), "-IM-SPIN", isp
            OPEN (unit=19, file=fonameoi, status='unknown')
            CALL write_CHGCAR_head(19, fonameoi, ikk, isp, iband, a1, a2, a3, &
                                   nx, ny, nz, rs)
        END IF

        ii = 0
        DO i3 = 1, nz
            DO i2 = 1, ny
                DO i1 = 1, nx
                    ii = ii + 1
                    DO iispinor = 1, ispinor
                        wkgr(:) = (wk(1) + ig(1, 1:ncnt)) * ((i1 - 1) / DBLE(nx) + rs(1)) + &
                                  (wk(2) + ig(2, 1:ncnt)) * ((i2 - 1) / DBLE(ny) + rs(2)) + &
                                  (wk(3) + ig(3, 1:ncnt)) * ((i3 - 1) / DBLE(nz) + rs(3))

                        csum(iispinor, ii) = SUM( &
                                             coeff(1 + npl(1) / 2 * (iispinor - 1):ncnt + npl(1) / 2 * (iispinor - 1)) * &
                                             cdexp(pi2 * CMPLX(0., 1.) * wkgr(1:ncnt)) * dsqrt(vol))
                    END DO !ispinor
                END DO !i1
            END DO !i2
            IF (MOD(i3, NINT(nz / 10.)) == 0) THEN
                WRITE (6, '(F5.1,A)', advance='yes') i3 / DBLE(nz) * 100, "%"
            END IF
        END DO !i3

        DO iline = 1, nline
            WRITE (18, '(5E15.7)') (REAL(csum(1, (iline - 1) * 5 + i)), i=1, 5)
        END DO
        IF (nresi .GE. 1) THEN
            WRITE (18, '(5E15.7)') (REAL(csum(1, (nline) * 5 + i)), i=1, nresi)
        END IF

        IF (ispinor .GE. 2) THEN
            WRITE (18, *) nx, ny, nz
            DO iline = 1, nline
                WRITE (18, '(5E15.7)') (REAL(csum(2, (iline - 1) * 5 + i)), i=1, 5)
            END DO
            IF (nresi .GE. 1) THEN
                WRITE (18, '(5E15.7)') (REAL(csum(2, (nline) * 5 + i)), i=1, nresi)
            END IF
        END IF

        IF (imag .EQ. 1) THEN
            DO iline = 1, nline
                WRITE (19, '(5E15.7)') (dimag(csum(1, (iline - 1) * 5 + i)), i=1, 5)
            END DO
            IF (nresi .GE. 1) THEN
                WRITE (19, '(5E15.7)') (dimag(csum(1, (nline) * 5 + i)), i=1, nresi)
            END IF

            IF (ispinor .GE. 2) THEN
                WRITE (19, *) nx, ny, nz
                DO iline = 1, nline
                    WRITE (19, '(5E15.7)') (dimag(csum(2, (iline - 1) * 5 + i)), i=1, 5)
                END DO
                IF (nresi .GE. 1) THEN
                    WRITE (19, '(5E15.7)') (dimag(csum(2, (nline) * 5 + i)), i=1, nresi)
                END IF
            END IF
        END IF

    END DO !ispin

    CLOSE (18)
    IF (imag .EQ. 1) CLOSE (19)
    STOP
END SUBROUTINE wf_plot

!!$*  subroutine for CHGCAR header (lattice info + atomic coordinates) writing
SUBROUTINE write_CHGCAR_head(ID, fonameo, ikk, isp, iband, a1, a2, a3, &
                             nx, ny, nz, rs)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION a1(3), a2(3), a3(3), ikk(1), coord(3)
    DIMENSION n_atom(10), rs(3), rss(3)
    CHARACTER * 4 at_name(10), const(3)
    CHARACTER * 75 fonameo
    CHARACTER dummy

    !get total number of atoms : read EIGENVAL header
    OPEN (unit=ID + 10, file='EIGENVAL', status='old', iostat=IERR)
    IF (IERR .NE. 0) WRITE (6, *) 'open error EIGENVAL - iostat =', IERR
    READ (ID + 10, *) natom, idummy, idummy, idummy
    CLOSE (ID + 10)

    !get atomic info: read POSCAR
    OPEN (unit=ID + 10, file='POSCAR', status='old', iostat=IERR)
    IF (IERR .NE. 0) WRITE (6, *) 'open error POSCAR - iostat =', IERR
    DO i = 1, 6; READ (ID + 10, *) dummy; END DO

    itot = 0
    itype = 1
    iselect = 0
    icont = 0
    n_atom = 0
    idirect = 0

    DO WHILE (ifinish .NE. 1)
        READ (ID + 10, *) n_atom(1:itype)
        itot = SUM(n_atom)
        IF (itot .NE. natom) THEN
            ifinish = 0
            itype = itype + 1
        ELSEIF (itot .EQ. natom) THEN
            ifinish = 1
        ELSEIF (itot .GT. natom) THEN
            WRITE (6, '(A,I4)') "ERROR: total number of atom exeeds ", natom
        END IF
        IF (ifinish .NE. 1) BACKSPACE (unit=ID + 10)
    END DO
    BACKSPACE (unit=ID + 10)
    BACKSPACE (unit=ID + 10)
    READ (ID + 10, *) at_name(1:itype)
    READ (ID + 10, *) dummy

    WRITE (ID, '(A,I3,A,I3,A,I1)') "WAVEFUNCTION: BAND= ", iband, &
        " ,KP= ", ikk(1), " SPIN= ", isp
    WRITE (ID, *) 1.0000
    WRITE (ID, '(3F20.16)') a1(:)
    WRITE (ID, '(3F20.16)') a2(:)
    WRITE (ID, '(3F20.16)') a3(:)
    WRITE (ID, *) at_name(1:itype)
    WRITE (ID, *) n_atom(1:itype)

    DO WHILE (icont .EQ. 0) ! check selective or dynamics
        READ (ID + 10, *) dummy
        IF (dummy == "S" .OR. dummy == "s") THEN
            WRITE (ID, '(A)') "Selective dynamics"
            icont = 0; iselect = 1
        ELSEIF (dummy == "D" .OR. dummy == "d") THEN
            WRITE (ID, '(A)') "Direct"
            icont = 1; idirect = 1
            rss(:) = rs(:)
        ELSEIF (dummy == "C" .OR. dummy == "c" .OR. dummy == "k" .OR. &
                dummy == "K") THEN
            WRITE (ID, '(A)') "Cartesian"
            icont = 1
            DO j = 1, 3
                rss(j) = rs(1) * a1(j) + rs(2) * a2(j) + rs(3) * a3(j)
            END DO
        END IF
    END DO

    DO i = 1, natom
        IF (iselect .EQ. 1) THEN
            READ (ID + 10, *) coord(1:3), const(1:3)
            IF (idirect .EQ. 1) THEN
                IF (coord(1) + rss(1) .GE. 1.) rss(1) = rss(1) - 1.
                IF (coord(2) + rss(2) .GE. 1.) rss(2) = rss(2) - 1.
                IF (coord(3) + rss(3) .GE. 1.) rss(3) = rss(3) - 1.
            END IF

            WRITE (ID, '(3F20.16,3X, A, A, A)') coord(:) + rss(:), const(1:3)
        ELSEIF (iselect .EQ. 0) THEN
            READ (ID + 10, *) coord(1:3)
            WRITE (ID, '(3F20.16)') coord(:) + rss(:)
        END IF
    END DO

    WRITE (ID, *) " "
    WRITE (ID, *) nx, ny, nz

    RETURN
END SUBROUTINE write_CHGCAR_head

!!$*  subroutine for setting k-points in BZ
SUBROUTINE set_BZ_fukui(nnk, wnklist, w_half_klist, i_half_klist, &
                        i_trim_klist, nhf, ispin, wklist, del, nk, iz2)

    IMPLICIT REAL * 8(a - h, o - z)
    REAL * 8 wklist(3, nk), wnklist(3, nk * iz2), w_half_klist(3, nk)
    INTEGER i_half_klist(nk), i_trim_klist(nk)
    DIMENSION wk(3)

! set new kpoints : TR partners
    i_trim_klist = 0
    DO isp = 1, ispin ! ispin start
        nnk = 0; nhf = 0
        DO ik = 1, nk !ik loop.
            wk(:) = wklist(:, ik)
            IF (ABS(wk(1)) .LE. del .AND. ABS(wk(2)) .LE. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !gamma
                w_half_klist(:, nhf) = wk(:) !save only B+ k-point's index
                i_half_klist(nhf) = ik
                i_trim_klist(nhf) = 4
            ELSEIF (ABS(wk(1) - .5) .LE. del .AND. ABS(wk(2) - .0) .LE. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !M1
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik
                i_trim_klist(nhf) = 1
            ELSEIF (ABS(wk(1) - .0) .LE. del .AND. ABS(wk(2) - .5) .LE. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !M2
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik
                i_trim_klist(nhf) = 2
            ELSEIF (ABS(wk(1) - .5) .LE. del .AND. ABS(wk(2) - .5) .LE. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !M3
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik
                i_trim_klist(nhf) = 3
            ELSEIF (wk(2) .GT. .0 .AND. wk(2) .LT. .5 &
                    .AND. (.5 - wk(2)) .GT. del &
                    .AND. (.5 - wk(1)) .GT. del &
                    .AND. (wk(2) - .0) .GT. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:) ! B+ interior
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik

                !find TR partner
                nnk = nnk + 1
                wnklist(:, nnk) = -wk(:)

!       elseif(wk(1).gt. del .and. (.5-wk(1)) .gt. del
            ELSEIF (wk(1) .LE. -del .AND. (ABS(wk(2) - .5) .LE. del)) THEN
!    &                          .and. (abs(wk(2) -.5) .le. del))then
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !B+ edge (upper)
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik

                !find TR + T(G2) partner
                nnk = nnk + 1
                wnklist(:, nnk) = -wk(:)
                wnklist(2, nnk) = wnklist(2, nnk) + 1.

            ELSEIF (wk(1) .GT. del .AND. (.5 - wk(1)) .GT. del &
                    .AND. (ABS(wk(2) - .0) .LE. del)) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !B+ edge (lower)
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik

                !find TR partner
                nnk = nnk + 1
                wnklist(:, nnk) = -wk(:)

            ELSEIF (ABS(.5 - wk(1)) .LE. del .AND. (.5 - wk(2)) .GT. del &
                    .AND. (wk(2) - .0) .GT. del) THEN
                nnk = nnk + 1
                nhf = nhf + 1
                wnklist(:, nnk) = wk(:)  !B+ edge (right)
                w_half_klist(:, nhf) = wk(:)
                i_half_klist(nhf) = ik

                !find TR + T(G1) partner
                nnk = nnk + 1
                wnklist(:, nnk) = -wk(:)
                wnklist(1, nnk) = wnklist(1, nnk) + 1.

            END IF
        END DO !ik
    END DO !ispin end

!     do i=1,nnk
!     write(6,*)(wnklist(jj,i),jj=1,3)
!     enddo
!     stop
    RETURN
END SUBROUTINE set_BZ_fukui

!!$*  subroutine for get nfield strength
SUBROUTINE get_nfield(rfield, wklp, isp, ik, nk, nband, nkx, nky, &
                      wnklist, ihf, nplist, nnk, iz, ns, kperiod, &
                      w_half_klist, i_half_klist, i_trim_klist, nhf, &
                      nbmax, npmax, ecut, ispinor, ispin, ne, irecl, &
                      nmax, nini, b1, b2, b3, iz2)
    IMPLICIT REAL * 8(a - h, o - z)
    COMPLEX * 8 coeff(npmax)
    COMPLEX*16, ALLOCATABLE :: Siju(:, :), Sijd(:, :), Sijt(:, :)
    DIMENSION nbmax(3), wkk(3, 5), wklp(3, 5, nk * iz2)
    COMPLEX * 16 coeff1u((2 * nbmax(1) + 2) * (2 * nbmax(2) + 2) * (2 * nbmax(3) + 1))
    COMPLEX * 16 coeff1d((2 * nbmax(1) + 2) * (2 * nbmax(2) + 2) * (2 * nbmax(3) + 1))
    COMPLEX * 16 coeff2u((2 * nbmax(1) + 2) * (2 * nbmax(2) + 2) * (2 * nbmax(3) + 1))
    COMPLEX * 16 coeff2d((2 * nbmax(1) + 2) * (2 * nbmax(2) + 2) * (2 * nbmax(3) + 1))
    COMPLEX * 16 detS(4), detA, detLOOP
    REAL * 8 wnklist(3, nk * iz2), w_half_klist(3, nk), rfield
    DIMENSION b1(3), b2(3), b3(3)
    INTEGER nplist(nk), i_half_klist(nk), i_trim_klist(nk)
    INTEGER itr(5), itrim(5), ikk(5), isgg(2, 5), npl(5)
    INTEGER ni, nj, ne, nk, nband, np, npmax, kperiod, ispin, irecl
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    pi = 4.*ATAN(1.)

    CALL klpfind(wkk, isp, ik, nk, nband, nkx, nky, wnklist, ihf, iz, iz2) !k-loop(4pts) of ik-kpoint (K),counter-clock
    CALL kindxfind(ikk, isgg, npl, itr, itrim, wnklist, nplist, wkk, nnk, &
                   nband, iz, w_half_klist, i_half_klist, i_trim_klist, &
                   nhf, iz2) !find index

    wklp(:, :, ik) = wkk(:, :)  ! closed loop (C) for each K : RECIVEC

!!$* get overlap matrix S_ij(k,k+1) over the C and PI_s [detS(k_s,k_s+1)], s=1,4
    ALLOCATE (Siju(ns, ns), Sijd(ns, ns), Sijt(ns, ns))
    DO ilp = 1, 4  ! loop for ilp
        Siju = (0., 0.)  !initialize
        Sijd = (0., 0.)
        Sijt = (0., 0.)
        ! construct overlap matrix S_ij(ikk1,ikk2)
        DO ni = nini, nmax   ! calculate upto valence band maximum
            ncnt = 0; coeff1u = (0., 0.); coeff1d = (0., 0.); coeff = (0., 0.)
            CALL get_coeff(coeff, isgg, ilp, ni, ikk, npl, nband, nk, isp, &
                           npmax, itrim)
            DO ig3 = 0, 2 * nbmax(3); ig3p = ig3
                IF (ig3 .GT. nbmax(3)) ig3p = ig3 - 2 * nbmax(3) - 1
                DO ig2 = 0, 2 * nbmax(2); ig2p = ig2
                    IF (ig2 .GT. nbmax(2)) ig2p = ig2 - 2 * nbmax(2) - 1
                    DO ig1 = 0, 2 * nbmax(1); ig1p = ig1
                        IF (ig1 .GT. nbmax(1)) ig1p = ig1 - 2 * nbmax(1) - 1
                        CALL get_etot(etot, ilp, ni, b1, b2, b3, isgg, ig1p, ig2p, ig3p, &
                                      ig1, ig2, ig3, wkk, itrim, itr)
                        IF (etot .LT. ecut) THEN; ncnt = ncnt + 1
                            CALL get_incnt(incnt, ilp, ni, ig1p, ig2p, ig3p, nbmax, &
                                           itr, itrim, isgg, wkk)
                            IF (ispinor .EQ. 2) THEN !spinor-dn
                                IF (itr(ilp) .EQ. 0) THEN
                                    IF (itrim(ilp) .EQ. 0) THEN
                                        coeff1d(incnt) = coeff(ncnt + npl(ilp) / 2) !for |u(dn)>
                                    ELSEIF (itrim(ilp) .GE. 1 .AND. MOD(ni, 2) .EQ. 1) THEN
                                        coeff1d(incnt) = coeff(ncnt + npl(ilp) / 2) !for |u(dn)>
                                    ELSEIF (itrim(ilp) .GE. 1 .AND. MOD(ni, 2) .EQ. 0) THEN
                                        coeff1d(incnt) = CONJG(coeff(ncnt)) !TRIM for |u*(2n-1,up)>
                                    END IF
                                ELSEIF (itr(ilp) .EQ. 1) THEN
                                    coeff1d(incnt) = CONJG(coeff(ncnt)) !TRS for |u*(up)>
                                END IF
                            END IF

                            IF (itr(ilp) .EQ. 0) THEN  !spinor-up
                                IF (itrim(ilp) .EQ. 0) THEN
                                    coeff1u(incnt) = coeff(ncnt) !for |u(up)>
                                ELSEIF (itrim(ilp) .GE. 1 .AND. MOD(ni, 2) .EQ. 1) THEN
                                    coeff1u(incnt) = coeff(ncnt) !for |u(up)>
                                ELSEIF (itrim(ilp) .GE. 1 .AND. MOD(ni, 2) .EQ. 0) THEN
                                    IF (ispinor .EQ. 1) THEN
                                        coeff1u(incnt) = CONJG(coeff(ncnt)) !TRIM for |u*(up)>
                                    ELSEIF (ispinor .EQ. 2) THEN
                                        coeff1u(incnt) = -CONJG(coeff(ncnt + npl(ilp) / 2)) !TRIM for -|u*(dn,2n-1)>
                                    END IF
                                END IF
                            ELSEIF (itr(ilp) .EQ. 1) THEN
                                IF (ispinor .EQ. 2) THEN
                                    coeff1u(incnt) = -CONJG(coeff(ncnt + npl(ilp) / 2)) !TRS for -|u*(dn)>
                                ELSEIF (ispinor .EQ. 1) THEN
                                    coeff1u(incnt) = CONJG(coeff(ncnt)) !TRS for |u*>
                                END IF
                            END IF

                        END IF
                    END DO  !loop for ig1
                END DO   !loop for ig2
            END DO    !loop for ig3
            IF (ispinor * ncnt .NE. npl(ilp)) THEN
                WRITE (0, *) '*ni error - computed NPL=', 2 * ncnt, &
                    ' != input nplane=', npl(ilp); STOP
            END IF
            DO nj = nini, nmax
                ncnt = 0; coeff2u = (0., 0.); coeff2d = (0., 0.); coeff = (0., 0.)
                CALL get_coeff(coeff, isgg, ilp + 1, nj, ikk, npl, nband, nk, isp, &
                               npmax, itrim)
                DO ig3 = 0, 2 * nbmax(3); ig3p = ig3
                    IF (ig3 .GT. nbmax(3)) ig3p = ig3 - 2 * nbmax(3) - 1
                    DO ig2 = 0, 2 * nbmax(2); ig2p = ig2
                        IF (ig2 .GT. nbmax(2)) ig2p = ig2 - 2 * nbmax(2) - 1
                        DO ig1 = 0, 2 * nbmax(1); ig1p = ig1
                            IF (ig1 .GT. nbmax(1)) ig1p = ig1 - 2 * nbmax(1) - 1
                            CALL get_etot(etot, ilp + 1, nj, b1, b2, b3, isgg, ig1p, ig2p, ig3p, &
                                          ig1, ig2, ig3, wkk, itrim, itr)
                            IF (etot .LT. ecut) THEN; ncnt = ncnt + 1
                                CALL get_incnt(incnt, ilp + 1, nj, ig1p, ig2p, ig3p, nbmax, &
                                               itr, itrim, isgg, wkk)

                                IF (ispinor .EQ. 2) THEN  !spinor-dn
                                    IF (itr(ilp + 1) .EQ. 0) THEN
                                        IF (itrim(ilp + 1) .GE. 0) THEN
                                            coeff2d(incnt) = coeff(ncnt + npl(ilp + 1) / 2) !c(dn,n)
                                        ELSEIF (itrim(ilp + 1) .GE. 1 .AND. MOD(nj, 2) .EQ. 1) THEN
                                            coeff2d(incnt) = coeff(ncnt + npl(ilp + 1) / 2) !c(dn,n)
                                        ELSEIF (itrim(ilp + 1) .GE. 1 .AND. MOD(nj, 2) .EQ. 0) THEN
                                            coeff2d(incnt) = CONJG(coeff(ncnt)) !c*(up,n-1)
                                        END IF
                                    ELSEIF (itr(ilp + 1) .EQ. 1) THEN
                                        coeff2d(incnt) = CONJG(coeff(ncnt)) !c*(up)
                                    END IF
                                END IF

                                IF (itr(ilp + 1) .EQ. 0) THEN  !spinor-up
                                    IF (itrim(ilp + 1) .EQ. 0) THEN
                                        coeff2u(incnt) = coeff(ncnt)
                                    ELSEIF (itrim(ilp + 1) .GE. 1 .AND. MOD(nj, 2) .EQ. 1) THEN
                                        coeff2u(incnt) = coeff(ncnt)
                                    ELSEIF (itrim(ilp + 1) .GE. 1 .AND. MOD(nj, 2) .EQ. 0) THEN
                                        IF (ispinor .EQ. 1) THEN
                                            coeff2u(incnt) = CONJG(coeff(ncnt))  !c(n-1)*
                                        ELSEIF (ispinor .EQ. 0) THEN
                                            coeff2u(incnt) = -CONJG(coeff(ncnt + npl(ilp + 1) / 2)) !-c*(dn,n-1)
                                        END IF
                                    END IF
                                ELSEIF (itr(ilp + 1) .EQ. 1) THEN
                                    IF (ispinor .EQ. 2) THEN
                                        coeff2u(incnt) = -CONJG(coeff(ncnt + npl(ilp + 1) / 2))
                                    ELSEIF (ispinor .EQ. 1) THEN
                                        coeff2u(incnt) = CONJG(coeff(ncnt))
                                    END IF
                                END IF

                            END IF
                        END DO   !loop for ig1
                    END DO    !loop for ig2
                END DO     !loop for ig3
                IF (ispinor * ncnt .NE. npl(ilp + 1)) THEN
                    WRITE (0, *) '*nj error - computed NPL=', 2 * ncnt, &
                        ' != input nplane=', npl(ilp + 1); STOP
                END IF
                IF (ispinor .EQ. 2) THEN
                    Siju(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1u, coeff2u)
                    Sijd(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1d, coeff2d)
                    Sijt(ni - nmax + ns, nj - nmax + ns) = Siju(ni - nmax + ns, nj - nmax + ns) + &
                                                           Sijd(ni - nmax + ns, nj - nmax + ns)
                ELSE IF (ispinor .EQ. 1) THEN
                    Sijt(ni - nmax + ns, nj - nmax + ns) = DOT_PRODUCT(coeff1u, coeff2u)
                END IF
            END DO !nj
        END DO !ni

        IF (nini .EQ. nmax) THEN   ! get determinant : det(S)
            detS(ilp) = Sijt(1, 1)
        ELSE
!         call getdetA(detA, Sijt,ns)
            CALL get_det(detA, Sijt, ns)
            detS(ilp) = detA; detA = (0., 0.)
        END IF
    END DO !ilp

    detLOOP = detS(1) * detS(2) * detS(3) * detS(4)
    rfield = +SUM(AIMAG(LOG(detS(1:4)))) / (2.*pi) &
             - (AIMAG(LOG(detLOOP))) / (2.*pi)
!       write(6,'(A)')"# "
! 704   format('#===>PI_S[det(S(K_s,K_s+1))] =',F16.8,'  +',F16.8,' i'
!    &         ,I4,'th-K')
    IF (SUM(itrim(1:4)) .GE. 1) THEN
        IF (itrim(1) + itrim(3) .GE. 1) nparity = 1
        IF (itrim(2) + itrim(4) .GE. 1) nparity = -1
        itr_sum = SUM(itrim(1:4)) * nparity
    ELSE
        itr_sum = 0
    END IF
    WRITE (6, 704) ik, rfield, LOG(detLOOP), ikk(:), wkk(:, ilp), &
        (wkk(:, 1) + wkk(:, 3)) / 2., itr_sum, SUM(LOG(detS(:)))
704 FORMAT('#IK=', I4, ' NK=', F5.2, ' P(detS)=', F16.8, '+', F16.8, 'i', &
           ' KLP=', 5I4, ' WK=', 3F8.4, ' WKC=', 3F8.4, ' STRM=', I2, &
           " S(detS)=", F16.8, '+', F16.8, 'i')

    DEALLOCATE (Sijt)
    DEALLOCATE (Siju)
    DEALLOCATE (Sijd)

    RETURN
END SUBROUTINE get_nfield

!!$*  subroutine for getting planewave coefficient for given k & n
SUBROUTINE get_coeff(coeff, isgg, iilp, nni, ikk, npl, nband, nk, isp, &
                     npmax, itrim)

    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION itrim(5), isgg(2, 5), ikk(5)
    COMPLEX * 8 coeff(npmax)
    INTEGER npl(5)
    irecl = 3 + (ikk(iilp) - 1) * (nband + 1) + nk * (nband + 1) * (isp - 1) + nni
    IF (itrim(iilp) .EQ. 0) THEN
        READ (10, rec=irecl) (coeff(i), i=1, npl(iilp))
    ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 1) THEN
        READ (10, rec=irecl) (coeff(i), i=1, npl(iilp))
    ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 0) THEN
        READ (10, rec=irecl - 1) (coeff(i), i=1, npl(iilp))
    END IF

    RETURN
END SUBROUTINE

!!$*  subroutine for getting incnt with given G vector
SUBROUTINE get_incnt(incnt, iilp, nni, ig1p, ig2p, ig3p, nbmax, &
                     itr, itrim, isgg, wkk)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION nbmax(3), itrim(5), itr(5), isgg(2, 5)
    REAL * 8 wkk(3, 5)
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    IF (itr(iilp) .EQ. 0) THEN
        IF (itrim(iilp) .EQ. 0) THEN
            incnt = (ig3p + nbmax(3)) * (2 * nbmax(2) + 1) * (2 * nbmax(1) + 1) + &
                    (ig2p - isgg(2, iilp) + nbmax(2)) * (2 * nbmax(1) + 1) + &
                    (ig1p - isgg(1, iilp) + nbmax(1)) + 1
        ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 1) THEN  ! TRIM & 2N-1 state
            incnt = (+ig3p + nbmax(3)) * (2 * nbmax(2) + 1) * (2 * nbmax(1) + 1) + &
                    (+ig2p - isgg(2, iilp) * 0 + nbmax(2)) * (2 * nbmax(1) + 1) + &
                    (+ig1p - isgg(1, iilp) * 0 + nbmax(1)) + 1
        ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 0) THEN  ! TRIM & 2N state
            incnt = (-ig3p + nbmax(3)) * (2 * nbmax(2) + 1) * (2 * nbmax(1) + 1) + &
                    (-ig2p - isgg(2, iilp) + nbmax(2)) * (2 * nbmax(1) + 1) + &
                    (-ig1p - isgg(1, iilp) + nbmax(1)) + 1
        END IF
    ELSEIF (itr(iilp) .EQ. 1) THEN !TIME REVERSAL
        incnt = (-ig3p + nbmax(3)) * (2 * nbmax(2) + 1) * (2 * nbmax(1) + 1) + &
                (-ig2p - isgg(2, iilp) + nbmax(2)) * (2 * nbmax(1) + 1) + &
                (-ig1p - isgg(1, iilp) + nbmax(1)) + 1
    END IF

    RETURN
END SUBROUTINE

!!$*  subroutine for getting etot within given K & G vector
SUBROUTINE get_etot(etot, iilp, nni, b1, b2, b3, isgg, ig1p, ig2p, ig3p, &
                    ig1, ig2, ig3, wkk, itrim, itr)

    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION itrim(5), itr(5), isgg(2, 5), sumkg(3)
    REAL * 8 b1(3), b2(3), b3(3), wkk(3, 5)
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    DO j = 1, 3
        IF (itr(iilp) .EQ. 0) THEN
            IF (itrim(iilp) .EQ. 0) THEN
                sumkg(j) = (wkk(1, iilp) + ig1p - isgg(1, iilp)) * b1(j) + &
                           (wkk(2, iilp) + ig2p - isgg(2, iilp)) * b2(j) + &
                           (wkk(3, iilp) + ig3p) * b3(j)
            ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 1) THEN
                sumkg(j) = (wkk(1, iilp) + ig1p - isgg(1, iilp) * 0) * b1(j) + &
                           (wkk(2, iilp) + ig2p - isgg(2, iilp) * 0) * b2(j) + &
                           (wkk(3, iilp) + ig3p) * b3(j)
            ELSEIF (itrim(iilp) .GE. 1 .AND. MOD(nni, 2) .EQ. 0) THEN
                sumkg(j) = (wkk(1, iilp) - ig1p - isgg(1, iilp) * 0) * b1(j) + &
                           (wkk(2, iilp) - ig2p - isgg(2, iilp) * 0) * b2(j) + &
                           (wkk(3, iilp) - ig3p) * b3(j)
            END IF
        ELSEIF (itr(iilp) .EQ. 1) THEN !TIME REVERSAL
            sumkg(j) = (-wkk(1, iilp) - ig1p - isgg(1, iilp)) * b1(j) + &
                       (-wkk(2, iilp) - ig2p - isgg(2, iilp)) * b2(j) + &
                       (-wkk(3, iilp) - ig3p) * b3(j)
        END IF
    END DO
    gtot = SQRT(DOT_PRODUCT(sumkg, sumkg))
    etot = gtot**2 / c

    RETURN
END SUBROUTINE

!!$*  subroutine for computing optical selectivity in the given k-point
SUBROUTINE optical_selectivity(selectivity, &
                               selectivitymax, selectivitymin, &
                               b1, b2, b3, wklist, isp, &
                               nband, ecut, ispinor, nplist, nbmax, npmax, nk, nini, nmax)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION nbmax(3), nplist(nk), wk(3)
    REAL * 8 selectivity(nk), b1(3), b2(3), b3(3), wklist(3, nk)
    REAL * 8 selectivitymax(4), selectivitymin(4)
    COMPLEX * 16 ctrans_mtrx_left, ctrans_mtrx_right
    COMPLEX * 16 cinter_mtrx_x, cinter_mtrx_y
    COMPLEX * 16 coeffv(npmax), coeffc(npmax)
    COMPLEX * 8 coeff(npmax)
    COMPLEX * 16 coeffvu(npmax), coeffvd(npmax)
    COMPLEX * 16 coeffcu(npmax), coeffcd(npmax)
    INTEGER :: ig(3, npmax)
!     data hbar/6.58211928e-16/
    DATA hbar/1./ !Here, I will set hbar = 1. for the simplicity
    DO ik = 1, nk
        cinter_mtrx_x = (0., 0.)
        cinter_mtrx_y = (0., 0.)
        coeff = (0., 0.)
        coeffv = (0., 0.); coeffc = (0., 0.)
        coeffvu = (0., 0.); coeffcu = (0., 0.)
        coeffvd = (0., 0.); coeffcd = (0., 0.)
        wk(:) = wklist(:, ik)
        np = nplist(ik)
        CALL plindx(ig, ncnt, ispinor, wk, b1, b2, b3, nbmax, np, ecut, npmax)
        READ (10, rec=(3 + (ik - 1) * (nband + 1) + &
                       nk * (nband + 1) * (isp - 1) + nini)) (coeff(i), i=1, np)
        coeffv = coeff; coeff = (0., 0.)
        READ (10, rec=(3 + (ik - 1) * (nband + 1) + &
                       nk * (nband + 1) * (isp - 1) + nmax)) (coeff(i), i=1, np)
        coeffc = coeff; coeff = (0., 0.)
        DO iplane = 1, ncnt
            xkgx = (wk(1) + ig(1, iplane)) * b1(1) + &
                   (wk(2) + ig(2, iplane)) * b2(1) + &
                   (wk(3) + ig(3, iplane)) * b3(1)
            xkgy = (wk(1) + ig(1, iplane)) * b1(2) + &
                   (wk(2) + ig(2, iplane)) * b2(2) + &
                   (wk(3) + ig(3, iplane)) * b3(2)
            IF (ispinor .EQ. 2) THEN
                coeffvu(iplane) = coeffv(iplane)
                coeffvd(iplane) = coeffv(iplane + ncnt)
                coeffcu(iplane) = coeffc(iplane)
                coeffcd(iplane) = coeffc(iplane + ncnt)
                cinter_mtrx_x = cinter_mtrx_x + &
                                hbar * CONJG(coeffcu(iplane)) * xkgx * coeffvu(iplane) + &
                                hbar * CONJG(coeffcd(iplane)) * xkgx * coeffvd(iplane)
                cinter_mtrx_y = cinter_mtrx_y + &
                                hbar * CONJG(coeffcu(iplane)) * xkgy * coeffvu(iplane) + &
                                hbar * CONJG(coeffcd(iplane)) * xkgy * coeffvd(iplane)
            ELSE IF (ispinor .EQ. 1) THEN
                cinter_mtrx_x = cinter_mtrx_x + &
                                hbar * CONJG(coeffc(iplane)) * xkgx * coeffv(iplane)
                cinter_mtrx_y = cinter_mtrx_y + &
                                hbar * CONJG(coeffc(iplane)) * xkgy * coeffv(iplane)
            END IF
        END DO ! iplane loop end

        ctrans_mtrx_left = cinter_mtrx_x + (0., 1.) * cinter_mtrx_y
        ctrans_mtrx_right = cinter_mtrx_x - (0., 1.) * cinter_mtrx_y
        selectivity(ik) = &
            ((ABS(ctrans_mtrx_left))**2 - (ABS(ctrans_mtrx_right))**2) / &
            ((ABS(ctrans_mtrx_left))**2 + (ABS(ctrans_mtrx_right))**2)
        WRITE (6, '(A,I4,4F11.6)') "# IK, K(reci), SELECTIVITY(n) : ", &
            ik, wk, selectivity(ik)
        IF (ik .eq. 1) THEN
            selectivitymax(4) = selectivity(ik)
            selectivitymax(1) = wklist(1, ik)
            selectivitymax(2) = wklist(2, ik)
            selectivitymax(3) = wklist(3, ik)
            selectivitymin(4) = selectivity(ik)
            selectivitymin(1) = wklist(1, ik)
            selectivitymin(2) = wklist(2, ik)
            selectivitymin(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. selectivity(ik) .GE. selectivitymax(4)) THEN
            selectivitymax(4) = selectivity(ik)
            selectivitymax(1) = wklist(1, ik)
            selectivitymax(2) = wklist(2, ik)
            selectivitymax(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. selectivity(ik) .LE. selectivitymin(4)) THEN
            selectivitymin(4) = selectivity(ik)
            selectivitymin(1) = wklist(1, ik)
            selectivitymin(2) = wklist(2, ik)
            selectivitymin(3) = wklist(3, ik)
        END IF
    END DO !ik loop end

    RETURN
END SUBROUTINE optical_selectivity

!!$*  subroutine for computing velocity expectation value for state psi(n,k)
SUBROUTINE vel_expectation(a1, a2, a3, b1, b2, b3, &
                           kperiod, nbmax, npmax, ecut, ispinor, ispin, nband, ne, nk, wklist, &
                           nplist, ni, nj, irecl, filename, foname)
    IMPLICIT REAL * 8(a - h, o - z)
    COMPLEX * 8 coeff(npmax)
    COMPLEX * 16 coeffi(npmax), coeffj(npmax)
    COMPLEX * 16 coeffiu(npmax), coeffid(npmax)
    COMPLEX * 16 coeffju(npmax), coeffjd(npmax)
    COMPLEX * 16 vel_x, vel_y
    REAL * 8 vel_expt(2, nk, ispin), wklist(3, nk)
    REAL * 8 recivec(3, nk), recilat(3, nk)
    REAL * 8 xvel_expt(2, kperiod * 2 * kperiod * 2 * nk, ispin)
    REAL * 8 xrecivec(3, kperiod * 2 * kperiod * 2 * nk)
    REAL * 8 xrecilat(3, kperiod * 2 * kperiod * 2 * nk)
    REAL * 8 vel_x_expt_max(4), vel_y_expt_max(4)
    REAL * 8 vel_x_expt_min(4), vel_y_expt_min(4)
    DIMENSION a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), a2xa3(3), sumkg(3)
    DIMENSION wk(3), nbmax(3), xb(3)
    INTEGER ig(3, npmax), nplist(nk)
    INTEGER ni, nj, ne, nk, nband, np, npmax, kperiod, ispin, irecl
    CHARACTER * 75 filename, foname, fonameo
    DATA c/0.262465831D0/ ! constant c = 2m/hbar**2 [1/eV Ang^2]
    DATA hbar/6.58211928E-16/ !h/2pi [eV * s]
    DATA xm/0.510998910E+6/ ! electron mass (eV/c^2)
    DATA anginv/1.0E+10/     ! inverse angstrom (1/Ang)
    pi = 4.*ATAN(1.)
    DO isp = 1, ispin
        DO ik = 1, nk
            coeff = (0., 0.)
            coeffi = (0., 0.); coeffj = (0., 0.)
            coeffiu = (0., 0.); coeffid = (0., 0.)
            coeffju = (0., 0.); coeffjd = (0., 0.)
            wk(:) = wklist(:, ik)
            np = nplist(ik)
            IF (ni .NE. nj) THEN
                CALL plindx(ig, ncnt, ispinor, wk, b1, b2, b3, nbmax, np, ecut, npmax)
                READ (10, rec=(3 + (ik - 1) * (nband + 1) + &
                               nk * (nband + 1) * (isp - 1) + ni)) (coeff(i), i=1, np)
                coeffi = coeff; coeff = (0., 0.)
                READ (10, rec=(3 + (ik - 1) * (nband + 1) + &
                               nk * (nband + 1) * (isp - 1) + nj)) (coeff(i), i=1, np)
                coeffj = coeff; coeff = (0., 0.)
            ELSE IF (ni .EQ. nj) THEN
                CALL plindx(ig, ncnt, ispinor, wk, b1, b2, b3, nbmax, np, ecut, npmax)
                READ (10, rec=(3 + (ik - 1) * (nband + 1) + &
                               nk * (nband + 1) * (isp - 1) + ni)) (coeff(i), i=1, np)
                coeffi = coeff; coeff = (0., 0.)
            END IF
            vel_x = (0., 0.); vel_y = (0., 0.)
            DO iplane = 1, ncnt
                xkgx = (wk(1) + ig(1, iplane)) * b1(1) + &
                       (wk(2) + ig(2, iplane)) * b2(1) + &
                       (wk(3) + ig(3, iplane)) * b3(1)
                xkgy = (wk(1) + ig(1, iplane)) * b1(2) + &
                       (wk(2) + ig(2, iplane)) * b2(2) + &
                       (wk(3) + ig(3, iplane)) * b3(2)
                IF (ispinor .EQ. 2) THEN
                    coeffiu(iplane) = coeffi(iplane)
                    coeffid(iplane) = coeffi(iplane + ncnt)
                    IF (ni .EQ. nj) THEN
                        coeffju(iplane) = coeffi(iplane)
                        coeffjd(iplane) = coeffi(iplane + ncnt)
                    ELSE
                        coeffju(iplane) = coeffj(iplane)
                        coeffjd(iplane) = coeffj(iplane + ncnt)
                    END IF
                    ! -i*hbar/m <psi|d/dx|psi>
                    vel_x = vel_x + &
                            anginv * hbar / xm * (CONJG(coeffiu(iplane)) * xkgx * coeffju(iplane) + &
                                                  CONJG(coeffid(iplane)) * xkgx * coeffjd(iplane))
                    ! -i*hbar/m <psi|d/dy|psi>
                    vel_y = vel_y + &
                            anginv * hbar / xm * (CONJG(coeffiu(iplane)) * xkgy * coeffju(iplane) + &
                                                  CONJG(coeffid(iplane)) * xkgy * coeffjd(iplane))
                ELSE IF (ispinor .EQ. 1) THEN
                    coeffiu(iplane) = coeffi(iplane)
                    IF (ni .EQ. nj) THEN
                        coeffju(iplane) = coeffi(iplane)
                    ELSE
                        coeffju(iplane) = coeffj(iplane)
                    END IF
                    vel_x = vel_x + &
                            anginv * hbar / xm * CONJG(coeffiu(iplane)) * xkgx * coeffju(iplane)
                    vel_y = vel_y + &
                            anginv * hbar / xm * CONJG(coeffiu(iplane)) * xkgy * coeffju(iplane)
                END IF ! ispinor
            END DO  ! iplane
            vel_expt(1, ik, isp) = REAL(vel_x, 8)
            vel_expt(2, ik, isp) = REAL(vel_y, 8)
            WRITE (6, '(A,I4,5F11.6)') "# IK, K(reci), VEL_EXPT(x,y) : ", &
                ik, wk, (vel_expt(i, ik, isp), i=1, 2)
        END DO  !ik

        IF (ik .EQ. 1) THEN
            vel_x_expt_max(4) = vel_expt(1, ik, isp)
            vel_x_expt_max(1) = wklist(1, ik)
            vel_x_expt_max(2) = wklist(2, ik)
            vel_x_expt_max(3) = wklist(3, ik)
            vel_x_expt_min(4) = vel_expt(1, ik, isp)
            vel_x_expt_min(1) = wklist(1, ik)
            vel_x_expt_min(2) = wklist(2, ik)
            vel_x_expt_min(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. vel_expt(1, ik, isp) .GE. vel_x_expt_max(4)) THEN
            vel_x_expt_max(4) = vel_expt(1, ik, isp)
            vel_x_expt_max(1) = wklist(1, ik)
            vel_x_expt_max(2) = wklist(2, ik)
            vel_x_expt_max(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. vel_expt(1, ik, isp) .LE. vel_x_expt_min(4)) THEN
            vel_x_expt_min(4) = vel_expt(1, ik, isp)
            vel_x_expt_min(1) = wklist(1, ik)
            vel_x_expt_min(2) = wklist(2, ik)
            vel_x_expt_min(3) = wklist(3, ik)
        END IF

        IF (ik .EQ. 1) THEN
            vel_y_expt_max(4) = vel_expt(2, ik, isp)
            vel_y_expt_max(1) = wklist(1, ik)
            vel_y_expt_max(2) = wklist(2, ik)
            vel_y_expt_max(3) = wklist(3, ik)
            vel_y_expt_min(4) = vel_expt(2, ik, isp)
            vel_y_expt_min(1) = wklist(1, ik)
            vel_y_expt_min(2) = wklist(2, ik)
            vel_y_expt_min(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. vel_expt(2, ik, isp) .GE. vel_y_expt_max(4)) THEN
            vel_y_expt_max(4) = vel_expt(2, ik, isp)
            vel_y_expt_max(1) = wklist(1, ik)
            vel_y_expt_max(2) = wklist(2, ik)
            vel_y_expt_max(3) = wklist(3, ik)
        ELSE IF (ik .GE. 2 .AND. vel_expt(2, ik, isp) .LE. vel_y_expt_min(4)) THEN
            vel_y_expt_min(4) = vel_expt(2, ik, isp)
            vel_y_expt_min(1) = wklist(1, ik)
            vel_y_expt_min(2) = wklist(2, ik)
            vel_y_expt_min(3) = wklist(3, ik)
        END IF
        DO ik = 1, nk
            DO j = 1, 3
                recivec(j, ik) = wklist(1, ik) * b1(j) + &
                                 wklist(2, ik) * b2(j) + &
                                 wklist(3, ik) * b3(j)
                recilat(j, ik) = wklist(j, ik)
            END DO
        END DO

        kk = 0   ! extend berry curvature distribution over extended BZ
        DO ib2 = -1 * (kperiod - 1) + 1, kperiod
            DO ib1 = -1 * (kperiod - 1) + 1, kperiod  ! you may adjust these values as you wish..
                DO ik = 1, nk
                    kk = kk + 1
                    xrecivec(1, kk) = recivec(1, ik) + (ib1 - 1) * b1(1) + (ib2 - 1) * b2(1)
                    xrecivec(2, kk) = recivec(2, ik) + (ib1 - 1) * b1(2) + (ib2 - 1) * b2(2)
                    xrecivec(3, kk) = recivec(3, ik) + (ib1 - 1) * b1(3) + (ib2 - 1) * b2(3)
                    xrecilat(1, kk) = recilat(1, ik) + (ib1 - 1)
                    xrecilat(2, kk) = recilat(2, ik) + (ib2 - 1)
                    xrecilat(3, kk) = recilat(3, ik)
                    xvel_expt(1, kk, isp) = vel_expt(1, ik, isp)
                    xvel_expt(2, kk, isp) = vel_expt(2, ik, isp)
                    kext = kk
                END DO
            END DO
        END DO

!$$*  sorting k-points and the corresponding optical selectivity on the
!     periodically repeated data grid
        WRITE (6, *) " "
        WRITE (6, '(A)') "# SORTING K-grids..."
        DO k = kext - 1, 1, -1  ! sorting kx
            DO j = 1, k
                IF (xrecivec(1, j + 1) .GT. xrecivec(1, j)) THEN
                    xb(:) = xrecivec(:, j)
                    xrecivec(:, j) = xrecivec(:, j + 1)
                    xrecivec(:, j + 1) = xb(:)
                    xb(:) = xrecilat(:, j)
                    xrecilat(:, j) = xrecilat(:, j + 1)
                    xrecilat(:, j + 1) = xb(:)
                    xtemp = xvel_expt(1, j, isp)
                    xvel_expt(1, j, isp) = xvel_expt(1, j + 1, isp)
                    xvel_expt(1, j + 1, isp) = xtemp
                    xtemp = xvel_expt(2, j, isp)
                    xvel_expt(2, j, isp) = xvel_expt(2, j + 1, isp)
                    xvel_expt(2, j + 1, isp) = xtemp
                END IF
            END DO
        END DO
        DO k = kext - 1, 1, -1  ! sorting ky
            DO j = 1, k
                IF (xrecivec(1, j + 1) .EQ. xrecivec(1, j)) THEN
                    IF (xrecivec(2, j + 1) .GT. xrecivec(2, j)) THEN
                        xb(:) = xrecivec(:, j)
                        xrecivec(:, j) = xrecivec(:, j + 1)
                        xrecivec(:, j + 1) = xb(:)
                        xb(:) = xrecilat(:, j)
                        xrecilat(:, j) = xrecilat(:, j + 1)
                        xrecilat(:, j + 1) = xb(:)
                        xtemp = xvel_expt(1, j, isp)
                        xvel_expt(1, j, isp) = xvel_expt(1, j + 1, isp)
                        xvel_expt(1, j + 1, isp) = xtemp
                        xtemp = xvel_expt(2, j, isp)
                        xvel_expt(2, j, isp) = xvel_expt(2, j + 1, isp)
                        xvel_expt(2, j + 1, isp) = xtemp
                    END IF
                END IF
            END DO
        END DO

        IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
            WRITE (fonameo, '(A,A)') TRIM(foname), '.UP.dat'
        ELSE IF (isp .EQ. 2 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 2) THEN
            WRITE (fonameo, '(A,A)') TRIM(foname), '.DN.dat'
        ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 2) THEN
            WRITE (fonameo, '(A,A)') TRIM(foname), '.dat'
        ELSE IF (isp .EQ. 1 .AND. ispinor .EQ. 1 .AND. ispin .EQ. 1) THEN
            WRITE (fonameo, '(A,A)') TRIM(foname), '.dat'
        END IF
        OPEN (61, file=fonameo, status='unknown')
        WRITE (61, '(A,A)') "# File reading... : ", filename
        WRITE (61, '(A,I9)') "# TOTAL RECORD LENGTH = ", irecl
        IF (ispinor .EQ. 2) THEN
            WRITE (61, '(A,I6,A)') "# ISPIN            : ", ispin, &
                " (LSORBIT = .TRUE.)"
        ELSE
            WRITE (61, '(A,I6,A)') "# ISPIN            : ", ispin, &
                " (LSORBIT = .FALSE.)"
        END IF
        WRITE (61, '(A,F11.4)') "# ENCUT (eV)       : ", ecut
        WRITE (61, '(A,I6)') "# NKPOINT          : ", nk
        WRITE (61, '(A,I6)') "# NBANDS           : ", nband
        WRITE (61, '(A,3F13.6)') "# RECIVEC B1 (A^-1): ", (b1(i), i=1, 3)
        WRITE (61, '(A,3F13.6)') "# RECIVEC B2       : ", (b2(i), i=1, 3)
        WRITE (61, '(A,3F13.6)') "# RECIVEC B3       : ", (b3(i), i=1, 3)
        WRITE (61, *) " "

        WRITE (61, '(A,I4)') "# VEOLOCITY EXPECTATION VALUE of BAND:", ni
        WRITE (61, '(A)') "# <v(n,k)>= 1/hbar dE(k)/dk =&
&      1/m_e<psi(n,k)|p|psi(n,k)>, p=-i*hbar*d/dx,m_e=elect_rest_mass"
        WRITE (61, '(A,4F16.6)') "# MAXVAL of VEL_EXPT <v_x> &
&     (in reci)= ", (vel_x_expt_max(i), i=1, 4)
        WRITE (61, '(A,4F16.6)') "# MINVAL of VEL_EXPT <v_x> &
&     (in reci)= ", (vel_x_expt_min(i), i=1, 4)
        WRITE (61, '(A,4F16.6)') "# MAXVAL of VEL_EXPT <v_y> &
&     (in reci)= ", (vel_y_expt_max(i), i=1, 4)
        WRITE (61, '(A,4F16.6)') "# MINVAL of VEL_EXPT <v_y> &
&     (in reci)= ", (vel_y_expt_min(i), i=1, 4)
        WRITE (61, '(A)') "# (cart) kx     ky     kz(A^-1)&
&         vel_expt(vx(n,k), vy(n,k))(m/s)  (recip)kx      ky      kz"
        DO ik = 1, kext
            WRITE (61, '(3F11.6,A,2F10.6,A,3F11.6)') (xrecivec(i, ik), i=1, 3), &
                "   ", (xvel_expt(i, ik, isp), i=1, 2), "       ", &
                (xrecilat(i, ik), i=1, 3)
        END DO
        CLOSE (61)

    END DO   !isp

!     write(6,'(A)')"# DONE! "
!     do isp=1, ispin
!      if(isp .eq. 1 .and. ispinor .eq. 1 .and. ispin .eq. 2) then
!       write(6,'(A,A,A)')"#  Results are summarized in ",TRIM(foname),
!    &                    ".UP.dat for spin-1"
!       else if (isp.eq.2 .and. ispinor.eq.1 .and. ispin.eq.2)then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".DN.dat for spin-2"
!       else if (isp .eq. 1 .and. ispinor .eq. 2) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!       else if (isp.eq.1 .and. ispinor.eq.1 .and. ispin.eq.1) then
!        write(6,'(A,A,A)')"#  Results are summarized in ",
!    &                     TRIM(foname),".dat"
!      endif
!     enddo
    CLOSE (10)
    RETURN
END SUBROUTINE vel_expectation

!!$*  subroutine for computing reciprocal lattice vector
SUBROUTINE recilatt(b1, b2, b3, dSkxky, a1, a2, a3, nkx, nky)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION a1(3), a2(3), a3(3), b1(3), b2(3), b3(3), a2xa3(3)
    DIMENSION b1xb2(3)
    pi = 4.*ATAN(1.)
    CALL vcross(a2xa3, a2, a3)
    Vcell = DOT_PRODUCT(a1, a2xa3)
    a3mag = dsqrt(DOT_PRODUCT(a3, a3))
    CALL vcross(b1, a2, a3); CALL vcross(b2, a3, a1); CALL vcross(b3, a1, a2)
    b1 = 2.*pi * b1 / Vcell; b2 = 2.*pi * b2 / Vcell; b3 = 2.*pi * b3 / Vcell

    CALL vcross(b1xb2, b1, b2)
    dSkxky = dsqrt(DOT_PRODUCT(b1xb2, b1xb2)) / REAL(nkx) / REAL(nky)
    RETURN
END SUBROUTINE recilatt

!!$*  subroutine for computing vector cross-product
SUBROUTINE vcross(a, b, c)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION a(3), b(3), c(3)
    a(1) = b(2) * c(3) - b(3) * c(2)
    a(2) = b(3) * c(1) - b(1) * c(3)
    a(3) = b(1) * c(2) - b(2) * c(1)
    RETURN
END SUBROUTINE vcross

!!$*  subroutine for computing reciprocal properties
SUBROUTINE reciproperty(nbmax, npmax, b1, b2, b3, ecut, ispinor)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION b1(3), b2(3), b3(3), vtmp(3), nbmax(3)
    DATA c/0.262465831D0/
    pi = 4.*ATAN(1.)

    b1mag = dsqrt(b1(1)**2 + b1(2)**2 + b1(3)**2)
    b2mag = dsqrt(b2(1)**2 + b2(2)**2 + b2(3)**2)
    b3mag = dsqrt(b3(1)**2 + b3(2)**2 + b3(3)**2)

    phi12 = ACOS((b1(1) * b2(1) + b1(2) * b2(2) + b1(3) * b2(3)) / (b1mag * b2mag))
    CALL vcross(vtmp, b1, b2)
    vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (b3(1) * vtmp(1) + b3(2) * vtmp(2) + b3(3) * vtmp(3)) / (vmag * b3mag)
    nb1maxA = (dsqrt(ecut * c) / (b1mag * ABS(SIN(phi12)))) + 1
    nb2maxA = (dsqrt(ecut * c) / (b2mag * ABS(SIN(phi12)))) + 1
    nb3maxA = (dsqrt(ecut * c) / (b3mag * ABS(sinphi123))) + 1
    npmaxA = NINT(4.*pi * nb1maxA * nb2maxA * nb3maxA / 3.)

    phi13 = ACOS((b1(1) * b3(1) + b1(2) * b3(2) + b1(3) * b3(3)) / (b1mag * b3mag))
    CALL vcross(vtmp, b1, b3)
    vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (b2(1) * vtmp(1) + b2(2) * vtmp(2) + b2(3) * vtmp(3)) / (vmag * b2mag)
    phi123 = ABS(ASIN(sinphi123))
    nb1maxB = (dsqrt(ecut * c) / (b1mag * ABS(SIN(phi13)))) + 1
    nb2maxB = (dsqrt(ecut * c) / (b2mag * ABS(sinphi123))) + 1
    nb3maxB = (dsqrt(ecut * c) / (b3mag * ABS(SIN(phi13)))) + 1
    npmaxB = NINT(4.*pi * nb1maxB * nb2maxB * nb3maxB / 3.)

    phi23 = ACOS((b2(1) * b3(1) + b2(2) * b3(2) + b2(3) * b3(3)) / (b2mag * b3mag))
    CALL vcross(vtmp, b2, b3)
    vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (b1(1) * vtmp(1) + b1(2) * vtmp(2) + b1(3) * vtmp(3)) / (vmag * b1mag)
    phi123 = ABS(ASIN(sinphi123))
    nb1maxC = (dsqrt(ecut * c) / (b1mag * ABS(sinphi123))) + 1
    nb2maxC = (dsqrt(ecut * c) / (b2mag * ABS(SIN(phi23)))) + 1
    nb3maxC = (dsqrt(ecut * c) / (b3mag * ABS(SIN(phi23)))) + 1
    npmaxC = NINT(4.*pi * nb1maxC * nb2maxC * nb3maxC / 3.)

    nbmax(1) = max0(nb1maxA, nb1maxB, nb1maxC)       ! maximum
    nbmax(2) = max0(nb2maxA, nb2maxB, nb2maxC)
    nbmax(3) = max0(nb3maxA, nb3maxB, nb3maxC)

      !! multiply 'ispinor' to handle two component spinors
    npmax = ispinor * min0(npmaxA, npmaxB, npmaxC)
    RETURN
END SUBROUTINE reciproperty

!!$*  subroutine for computing planewave G index
SUBROUTINE plindx(ig, ncnt, &
                  ispinor, wk, b1, b2, b3, nbmax, np, ecut, npmax)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION wk(3), sumkg(3), b1(3), b2(3), b3(3), nbmax(3)
    INTEGER :: ig(3, npmax)
    DATA c/0.262465831D0/
    ncnt = 0
    DO ig3 = 0, 2 * nbmax(3)
        ig3p = ig3
        IF (ig3 .GT. nbmax(3)) ig3p = ig3 - 2 * nbmax(3) - 1
        DO ig2 = 0, 2 * nbmax(2)
            ig2p = ig2
            IF (ig2 .GT. nbmax(2)) ig2p = ig2 - 2 * nbmax(2) - 1
            DO ig1 = 0, 2 * nbmax(1)
                ig1p = ig1
                IF (ig1 .GT. nbmax(1)) ig1p = ig1 - 2 * nbmax(1) - 1
                DO j = 1, 3
                    sumkg(j) = (wk(1) + ig1p) * b1(j) + &
                               (wk(2) + ig2p) * b2(j) + &
                               (wk(3) + ig3p) * b3(j)
                END DO
                gtot = SQRT(DOT_PRODUCT(sumkg, sumkg))
                etot = gtot**2 / c
                IF (etot .LT. ecut) THEN
                    ncnt = ncnt + 1
                    ig(1, ncnt) = ig1p
                    ig(2, ncnt) = ig2p
                    ig(3, ncnt) = ig3p
                END IF
            END DO
        END DO
    END DO
    IF (ispinor * ncnt .NE. np) THEN
        WRITE (0, *) '*** error - computed ispinor*ncnt=', ispinor * ncnt, &
            ' != input nplane=', np; STOP
    END IF
    RETURN
END SUBROUTINE plindx

!!$*  subroutine for finding k-loop set for certain "k"-point
SUBROUTINE klpfind(wkk, isp, ik, nk, nband, nkx, nky, wnklist, ihf, iz, &
                   iz2)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION wkk(3, 5), wk(3), dk1(3), dk2(3)
    REAL * 8 occ(nband), wnklist(3, nk * iz2)
    REAL * 16 ener(nband)
    dk1 = 0.; dk2 = 0. !;ne=0
    dk1(1) = 1./REAL(nkx, 8); dk2(2) = 1./REAL(nky, 8)  !vector for k-grid
! asign k-loop set wklp for each k-point
    IF (iz .EQ. 0) THEN
        IF (ihf .EQ. 0) THEN
            CALL kread(wk, nplane, ener, occ, isp, ik, nk, nband)
        ELSEIF (ihf .EQ. 1) THEN
            wk(:) = wnklist(:, ik)
        END IF
        wkk(:, 1) = wk(:)
        wkk(:, 2) = wk(:) + dk1(:)
        wkk(:, 3) = wk(:) + dk1(:) + dk2(:)
        wkk(:, 4) = wk(:) + dk2(:)
        wkk(:, 5) = wk(:)
    ELSEIF (iz .EQ. 1) THEN
        IF (ihf .EQ. 0) THEN
            CALL kread(wk, nplane, ener, occ, isp, ik, nk, nband)
        ELSEIF (ihf .EQ. 1) THEN
            wk(:) = wnklist(:, ik)
        END IF
        wkk(:, 1) = wk(:)
        wkk(:, 2) = wk(:) - dk1(:)
        wkk(:, 3) = wk(:) - dk1(:) - dk2(:)
        wkk(:, 4) = wk(:) - dk2(:)
        wkk(:, 5) = wk(:)
    END IF

!     do i=1,5
!     write(6,'(i4,3F10.4)')ik, (wkk(jj,i),jj=1,3)
!     enddo
    RETURN
END SUBROUTINE klpfind

!!$*  subroutine for finding k-point index for certain wk
SUBROUTINE kindxfind(ikk, isgg, npl, itr, itrim, wklist, nplist, wkk, &
                     nk, nband, iz, w_half_klist, &
                     i_half_klist, i_trim_klist, nhf, iz2)

    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION wkk(3, 5), ikk(5), isgg(2, 5), wklist(3, nk * iz2), itr(5)
    DIMENSION npl(5), nplist(nk), itrim(5)
    REAL * 8 w_half_klist(3, nk)
    INTEGER i_half_klist(nk), i_trim_klist(nk)
    del = 1E-6 ! criterion for k-point find

    ikk = 0; isgg = 0; itr = 0; itrim = 0
!     ! find k-point index ikk for wkk set
    IF (iz .EQ. 0) THEN
        DO ilp = 1, 5   !ilp
            DO ik = 1, nk  !ik
                d1 = wklist(1, ik) - wkk(1, ilp)
                d2 = wklist(2, ik) - wkk(2, ilp)
                d3 = wklist(3, ik) - wkk(3, ilp)
                distk = dsqrt(d1**2 + d2**2 + d3**2)
                IF (distk .LT. 1E-5) THEN
                    ikk(ilp) = ik
                    npl(ilp) = nplist(ik)
                ELSE IF (ABS(1.-distk) .LT. 1E-5) THEN
                    IF (ABS(1.-dsqrt(d1**2)) .LT. 1E-5) THEN
                        ikk(ilp) = ik
                        npl(ilp) = nplist(ik)
                        isgg(1, ilp) = 1
                    ELSE IF (ABS(1.-dsqrt(d2**2)) .LT. 1E-5) THEN
                        ikk(ilp) = ik
                        npl(ilp) = nplist(ik)
                        isgg(2, ilp) = 1
                    END IF
                ELSE IF (ABS(SQRT(2.) - distk) .LT. 1E-5) THEN
                    ikk(ilp) = ik
                    isgg(1, ilp) = 1
                    isgg(2, ilp) = 1
                    npl(ilp) = nplist(ik)
                END IF
            END DO !ik
        END DO  !ilp

    ELSEIF (iz .EQ. 1) THEN
        DO ilp = 1, 5 !ilp
            DO ik = 1, nhf !ik
                d1 = w_half_klist(1, ik) - wkk(1, ilp)
                d2 = w_half_klist(2, ik) - wkk(2, ilp)
                d3 = w_half_klist(3, ik) - wkk(3, ilp)
                distk = dsqrt(d1**2 + d2**2 + d3**2)      !B+ region including TRIM and edge (all nhf)

                d1t = -w_half_klist(1, ik) - wkk(1, ilp)
                d2t = -w_half_klist(2, ik) - wkk(2, ilp)
                d3t = -w_half_klist(3, ik) - wkk(3, ilp)
                distkt = dsqrt(d1t**2 + d2t**2 + d3t**2)  !time-reversal

                d1tx = -w_half_klist(1, ik) + 1.-wkk(1, ilp)
                d2tx = -w_half_klist(2, ik) - wkk(2, ilp)
                d3tx = -w_half_klist(3, ik) - wkk(3, ilp)
                distktx = dsqrt(d1tx**2 + d2tx**2 + d3tx**2) !time-reversal+G1(B- right edge or B- right outside)

                d1ty = -w_half_klist(1, ik) - wkk(1, ilp)
                d2ty = -w_half_klist(2, ik) + 1.-wkk(2, ilp)
                d3ty = -w_half_klist(3, ik) - wkk(3, ilp)
                distkty = dsqrt(d1ty**2 + d2ty**2 + d3ty**2) !time-reversal+G2(B+ upper edge or B+ upper ouside)

                d1txy = -w_half_klist(1, ik) + 1.-wkk(1, ilp)
                d2txy = -w_half_klist(2, ik) + 1.-wkk(2, ilp)
                d3txy = -w_half_klist(3, ik) - wkk(3, ilp)
                distktxy = dsqrt(d1txy**2 + d2txy**2 + d3txy**2) !time-reversal+G1+G2 (B+ upper left outside corner)

                d1x = w_half_klist(1, ik) + 1.-wkk(1, ilp)
                d2x = w_half_klist(2, ik) - wkk(2, ilp)
                d3x = w_half_klist(3, ik) - wkk(3, ilp)
                distkx = dsqrt(d1x**2 + d2x**2 + d3x**2) !G1 (B+ right edge)

                d1x_ = w_half_klist(1, ik) - 1.-wkk(1, ilp)
                d2x_ = w_half_klist(2, ik) - wkk(2, ilp)
                d3x_ = w_half_klist(3, ik) - wkk(3, ilp)
                distkx_ = dsqrt(d1x_**2 + d2x_**2 + d3x_**2) !G1 (B+ left edge)

                d1y = w_half_klist(1, ik) - wkk(1, ilp)
                d2y = w_half_klist(2, ik) - 1.-wkk(2, ilp)
                d3y = w_half_klist(3, ik) - wkk(3, ilp)
                distky_ = dsqrt(d1y**2 + d2y**2 + d3y**2) !-G2 (B- bottom edge)

                d1m1 = -0.5 - wkk(1, ilp)
                d2m1 = 0.0 - wkk(2, ilp)
                d3m1 = 0.0 - wkk(3, ilp)
                distkm1 = dsqrt(d1m1**2 + d2m1**2 + d3m1**2) !-M1

                d1m2 = 0.0 - wkk(1, ilp)
                d2m2 = -0.5 - wkk(2, ilp)
                d3m2 = 0.0 - wkk(3, ilp)
                distkm2 = dsqrt(d1m2**2 + d2m2**2 + d3m2**2) !-M2

                d1m3 = -0.5 - wkk(1, ilp)
                d2m3 = -0.5 - wkk(2, ilp)
                d3m3 = 0.0 - wkk(3, ilp)
                distkm3 = dsqrt(d1m3**2 + d2m3**2 + d3m3**2) !-M3

                d1m31 = 0.5 - wkk(1, ilp)
                d2m31 = -0.5 - wkk(2, ilp)
                d3m31 = 0.0 - wkk(3, ilp)
                distkm31 = dsqrt(d1m31**2 + d2m31**2 + d3m31**2) !-M3+G1

                d1m32 = -0.5 - wkk(1, ilp)
                d2m32 = 0.5 - wkk(2, ilp)
                d3m32 = 0.0 - wkk(3, ilp)
                distkm32 = dsqrt(d1m32**2 + d2m32**2 + d3m32**2) !-M3+G2

                IF (distk .LT. 1E-5) THEN
                    ikk(ilp) = i_half_klist(ik) !Original
                    npl(ilp) = nplist(i_half_klist(ik))
                    IF (i_trim_klist(ik) .GE. 1) THEN
                        itrim(ilp) = i_trim_klist(ik)
                    END IF
                ELSEIF (distkt .LT. 1E-5 .AND. distk .GT. 1E-5) THEN
                    IF (distkm1 .LT. 1E-5) THEN
                        ikk(ilp) = i_half_klist(ik) ! -M1
                        npl(ilp) = nplist(i_half_klist(ik))
                        itrim(ilp) = i_trim_klist(ik)
                        isgg(1, ilp) = -1
                    ELSEIF (distkm2 .LT. 1E-5) THEN
                        ikk(ilp) = i_half_klist(ik) ! -M2
                        npl(ilp) = nplist(i_half_klist(ik))
                        itrim(ilp) = i_trim_klist(ik)
                        isgg(2, ilp) = -1
                    ELSEIF (distkm3 .LT. 1E-5) THEN
                        ikk(ilp) = i_half_klist(ik) ! -M3
                        npl(ilp) = nplist(i_half_klist(ik))
                        itrim(ilp) = i_trim_klist(ik)
                        isgg(:, ilp) = -1
                    ELSE
                        ikk(ilp) = i_half_klist(ik) ! TRS
                        npl(ilp) = nplist(i_half_klist(ik))
                        itr(ilp) = 1
                    END IF
                ELSEIF (distktx .LT. 1E-5 .AND. distk .GT. 1E-5) THEN
                    IF (distkm31 .LT. 1E-5) THEN
                        ikk(ilp) = i_half_klist(ik) ! M3-G2
                        npl(ilp) = nplist(i_half_klist(ik))
                        itrim(ilp) = i_trim_klist(ik)
                        isgg(2, ilp) = -1
                    ELSE
                        ikk(ilp) = i_half_klist(ik) !TRS + G1
                        npl(ilp) = nplist(i_half_klist(ik))
                        itr(ilp) = 1
                        isgg(1, ilp) = 1
                    END IF
                ELSEIF (distkty .LT. 1E-5 .AND. distk .GT. 1E-5) THEN
                    IF (distkm32 .LT. 1E-5) THEN
                        ikk(ilp) = i_half_klist(ik) ! M3-G1
                        npl(ilp) = nplist(i_half_klist(ik))
                        itrim(ilp) = i_trim_klist(ik)
                        isgg(1, ilp) = -1
                    ELSE
                        ikk(ilp) = i_half_klist(ik) ! TRS + G2
                        npl(ilp) = nplist(i_half_klist(ik))
                        itr(ilp) = 1
                        isgg(2, ilp) = 1
                    END IF
                ELSEIF (distkx_ .LT. 1E-5 .AND. distk .GT. 1E-5 &
                        .AND. wkk(2, ilp) .GT. del &
                        .AND. ABS(0.5 - wkk(2, ilp)) .GT. del) THEN
                    ikk(ilp) = i_half_klist(ik) ! B+ left edge
                    npl(ilp) = nplist(i_half_klist(ik))
                    isgg(1, ilp) = -1
                ELSEIF (distky_ .LT. 1E-5 .AND. distk .GT. 1E-5 &
                        .AND. wkk(1, ilp) .GT. -0.5 + del &
                        .AND. wkk(1, ilp) .LT. -del) THEN
                    ikk(ilp) = i_half_klist(ik) ! B- bottom edge
                    npl(ilp) = nplist(i_half_klist(ik))
                    isgg(2, ilp) = -1
                END IF

            END DO !ik
        END DO !ilp
    END IF

    RETURN
END SUBROUTINE kindxfind

!!$*  subroutine for finding k-point vector in the given k index
SUBROUTINE kread(wk, nplane, ener, occ, isp, k, nk, nband)
    IMPLICIT REAL * 8(a - h, o - z)
    DIMENSION wk(3)
    COMPLEX * 16 cener(nband)
    REAL * 8 occ(nband)
    REAL * 16 ener(nband)
    irec = 3 + (k - 1) * (nband + 1) + nk * (nband + 1) * (isp - 1)  !record addres for "k"-point
    READ (10, rec=irec) xnplane, (wk(i), i=1, 3), &
        (cener(nn), occ(nn), nn=1, nband)
    nplane = NINT(xnplane); ener = REAL(cener)
    IF (kpoint .GT. nk) THEN
        WRITE (0, *) '*** error - selected k=', k, ' > max k=', nk; STOP
    END IF

    RETURN
END SUBROUTINE kread

!This function is to calculate determinant of the complex matrix
!The source is adoped from : https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
SUBROUTINE get_det(determinant, mat, N)
    IMPLICIT NONE
    INTEGER*4, INTENT(in) :: N
    COMPLEX * 16 mat(N, N)
    INTEGER * 4 i, info
    INTEGER * 4 ipiv(N)
    COMPLEX * 16 determinant
    REAL * 8 sgn

    ipiv = 0
    CALL zgetrf(N, N, mat, N, ipiv, info)

    determinant = (1D0, 0D0)
    DO i = 1, N
        determinant = determinant * mat(i, i)
    END DO

    sgn = 1D0
    DO i = 1, N
        IF (ipiv(i) /= i) THEN
            sgn = -sgn
        END IF
    END DO
    determinant = sgn * determinant

end subroutine

!$$*  subroutine to get det(A) where A is n x n complex matrix
!   This subroutine is adopted from below,
!   "http://computer-programming-forum.com/49-fortran/9e079d718158c944.htm"
!   The author is : Che-Ping Su
!c****************************************************
!    Description:     Calculate the determinant of
!                     general complex matrix
!    Input:           A - matrix to be calculated
!                     N - the order of matrix A
!    Output:          DT - Determinant of matrix A
!    Last Updated:      03/17/1998
!****************************************************
SUBROUTINE getdetA(DT, S, N)
    IMPLICIT REAL * 8(a - h, o - z)
    COMPLEX * 16 S(N, N), DT, TM, TC, W(N, 1)
    COMPLEX * 16 A(N, N)
    A = (0., 0.)
    A = S
    L = 1
    K = 2
    flag = (0., 0.)
    DT = CMPLX(1.)
10  TM = A(L, L)
    IF (A(L, L) .EQ. (0.0, 0.0)) THEN
        DO I = L + 1, N
            IF (A(L, I) .NE. (0.0, 0.0)) THEN
                DO J = 1, N
                    w(J, 1) = A(J, L)
                    A(J, L) = A(J, I)
                    A(J, I) = W(J, 1)
                END DO
                flag = flag + (1., 0.)
                GOTO 10
            END IF
        END DO
    END IF

    DO 20 J = L, N
        A(L, J) = A(L, J) / TM
20      CONTINUE
        DO 30 I = k, N
            TC = A(I, L)
            DO 30 J = L, N
                A(I, J) = A(I, J) - A(L, J) * TC
30              CONTINUE
                L = L + 1
                k = k + 1
                DT = DT * TM
                IF (L - N) 10, 40, 40
40              DT = (-1., 0.)**flag * DT * A(N, N)
                RETURN
                END SUBROUTINE getdetA

!!$   parse command line arguments
                SUBROUTINE parse(filename, foname, nkx, nky, ispinor, icd, ixt, fbz, &
                                 ivel, iz, ihf, nini, nmax, kperiod, it, iskp, ine, ver_tag, &
                                 iwf, ikwf, ng, rs, imag)
                    IMPLICIT REAL * 8(a - h, o - z)
                    CHARACTER * 75 filename, foname, fbz, ver_tag, vdirec
                    REAL * 8 x, y
                    CHARACTER * 20 option, VALUE
                    INTEGER iarg, narg, ia, nkx, nky, ispinor, iskp, ine, ng(3)
                    DIMENSION rs(3)

                    nini = 1; it = 0; iskp = 0; ine = 0; icd = 0; ixt = 0; ivel = 0; iz = 0; ihf = 0
                    iwf = 0; ikwf = 1; ng = 0; imag = 0; rs = 0.
                    nmax = 999999
                    iarg = iargc()
                    nargs = iarg / 2
                    filename = "WAVECAR"
                    foname = "BERRYCURV"
                    fbz = "BERRYCURV.tot.dat"
                    IF (iarg .NE. 2 * nargs) THEN
                        CALL help(ver_tag)
                    END IF
                    DO ia = 1, nargs
                        CALL getarg(2 * ia - 1, option)
                        CALL getarg(2 * ia, VALUE)
                        IF (option == "-f") THEN
                            READ (VALUE, *) filename
                        ELSE IF (option == "-o") THEN
                            READ (VALUE, *) foname
                        ELSE IF (option == "-kx") THEN
                            READ (VALUE, *) nkx
                        ELSE IF (option == "-ky") THEN
                            READ (VALUE, *) nky
                        ELSE IF (option == "-s") THEN
                            READ (VALUE, *) ispinor
                        ELSE IF (option == "-ii") THEN
                            READ (VALUE, *) nini
                        ELSE IF (option == "-if") THEN
                            READ (VALUE, *) nmax
                        ELSE IF (option == "-is") THEN
                            READ (VALUE, *) nini
                            nini = nini; nmax = nini
                        ELSE IF (option == "-kp") THEN
                            READ (VALUE, *) kperiod
                        ELSE IF (option == "-t") THEN
                            READ (VALUE, *) it
                        ELSE IF (option == "-skp") THEN
                            READ (VALUE, *) iskp
                        ELSE IF (option == "-ne") THEN
                            READ (VALUE, *) ine
                        ELSE IF (option == "-ixt") THEN
                            READ (VALUE, *) ixt
                        ELSE IF (option == "-fbz") THEN
                            READ (VALUE, *) fbz
                        ELSE IF (option == "-cd") THEN
                            READ (VALUE, *) icd
                        ELSE IF (option == "-vel") THEN
                            READ (VALUE, *) ivel
                        ELSE IF (option == "-z2") THEN
                            READ (VALUE, *) iz
                        ELSE IF (option == "-hf") THEN
                            READ (VALUE, *) ihf
                        ELSE IF (option == "-wf") THEN
                            READ (VALUE, *) iwf
                        ELSE IF (option == "-k") THEN
                            READ (VALUE, *) ikwf
                        ELSE IF (option == "-ng") THEN
                            READ (VALUE(:), *) ng(1:3)
                        ELSE IF (option == "-ishift") THEN
                            READ (VALUE(:), *) rs(1:3)
                        ELSE IF (option == "-im") THEN
                            READ (VALUE, *) imag
                        ELSE IF (option == "-h") THEN
                            CALL help(ver_tag)
                        ELSE
                            CALL help(ver_tag)
                        END IF
                    END DO
                    IF (icd .EQ. 1 .AND. TRIM(foname) .NE. 'BERRYCURV') THEN
                        WRITE (foname, '(A,A)') "CIRC_DICHROISM.", TRIM(foname)
                    ELSE IF (icd .EQ. 1 .AND. TRIM(foname) .EQ. 'BERRYCURV') THEN
                        foname = "CIRC_DICHROISM"
                    ELSE IF (icd + ivel .EQ. 0 .AND. TRIM(foname) .NE. 'BERRYCURV') THEN
                        WRITE (foname, '(A,A)') "BERRYCURV.", TRIM(foname)
                    ELSE IF (ivel .EQ. 1 .AND. TRIM(foname) .NE. 'BERRYCURV') THEN
                        WRITE (foname, '(A,A)') "VEL_EXPT.", TRIM(foname)
                    ELSE IF (ivel .EQ. 1 .AND. TRIM(foname) .EQ. 'BERRYCURV') THEN
                        foname = "VEL_EXPT"
                    ELSE IF (iz .EQ. 1 .AND. TRIM(foname) .EQ. 'BERRYCURV') THEN
                        foname = "NFIELD"
                    END IF

                    IF (iwf .GE. 1) THEN
                        ine = iwf
                        IF (iwf .LT. 10) THEN
                            IF (ikwf .LT. 10) THEN
                                WRITE (foname, '(A,I1,A,I1)') "PARCHG-W-K00", ikwf, "-E00", iwf
                            ELSEIF (ikwf .GE. 10 .AND. ikwf .LT. 100) THEN
                                WRITE (foname, '(A,I2,A,I1)') "PARCHG-W-K0", ikwf, "-E00", iwf
                            ELSEIF (ikwf .GE. 100) THEN
                                WRITE (foname, '(A,I3,A,I1)') "PARCHG-W-K", ikwf, "-E00", iwf
                            END IF
                        ELSEIF (iwf .GE. 10 .AND. iwf .LT. 100) THEN
                            IF (ikwf .LT. 10) THEN
                                WRITE (foname, '(A,I1,A,I2)') "PARCHG-W-K00", ikwf, "-E0", iwf
                            ELSEIF (ikwf .GE. 10 .AND. ikwf .LT. 100) THEN
                                WRITE (foname, '(A,I2,A,I2)') "PARCHG-W-K0", ikwf, "-E0", iwf
                            ELSEIF (ikwf .GE. 100) THEN
                                WRITE (foname, '(A,I3,A,I2)') "PARCHG-W-K", ikwf, "-E0", iwf
                            END IF
                        ELSEIF (iwf .GE. 100) THEN
                            IF (ikwf .LT. 10) THEN
                                WRITE (foname, '(A,I1,A,I3)') "PARCHG-W-K00", ikwf, "-E", iwf
                            ELSEIF (ikwf .GE. 10 .AND. ikwf .LT. 100) THEN
                                WRITE (foname, '(A,I2,A,I3)') "PARCHG-W-K0", ikwf, "-E", iwf
                            ELSEIF (ikwf .GE. 100) THEN
                                WRITE (foname, '(A,I3,A,I3)') "PARCHG-W-K", ikwf, "-E", iwf
                            END IF
                        END IF
                    END IF

                    RETURN
                END SUBROUTINE parse

!!$*  subroutine for reading basic information
                SUBROUTINE inforead(irecl, ispin, nk, nband, ecut, a1, a2, a3, filename)
                    IMPLICIT REAL * 8(a - h, o - z)
                    CHARACTER * 75 filename
                    DIMENSION a1(3), a2(3), a3(3)

                    irecl = 24
                    OPEN (unit=10, file=filename, access='direct', recl=irecl, &
                          iostat=iost, status='old')
                    IF (iost .NE. 0) WRITE (6, *) '0.open error - iostat =', iost

                    READ (10, rec=1) xirecl, xispin, xiprec !RDUM,RISPIN,RTAG(in real type)
                    CLOSE (10)
                    irecl = NINT(xirecl); ispin = NINT(xispin); iprec = NINT(xiprec) ! set to integer
                    IF (iprec .EQ. 45210) THEN
                        WRITE (0, *) '*** error - WAVECAR_double requires complex*16'; STOP
                    END IF
                    OPEN (unit=10, file=filename, access='direct', recl=irecl, &
                          iostat=iost, status='old')
                    IF (iost .NE. 0) WRITE (6, *) '1.open error - iostat =', iost
                    READ (10, rec=2) xnk, xnband, ecut,                 !RNKPTS,RNB_TOT,ENCUT&
                    (a1(j), j=1, 3), (a2(j), j=1, 3), (a3(j), j=1, 3)       !A1(3),A2(3),A3(3)
                    nk = NINT(xnk)
                    nband = NINT(xnband)

                    RETURN
                END SUBROUTINE inforead

                SUBROUTINE creditinfo(ver_tag)
                    CHARACTER * 75 ver_tag

                    WRITE (6, *) ver_tag
                    WRITE (6, *) "#This program calculates (1)berry curvature omega(k) "
                    WRITE (6, *) "#for closed loop C on a small patches in k-space,"
                    WRITE (6, *) "#and (2) degree of optical selectivity between "
                    WRITE (6, *) "#two bands specified. "
                    WRITE (6, *) " "
                    WRITE (6, *) "#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                    WRITE (6, *) "#! Copyright 2015. Hyun-Jung Kim All rights reserved.!"
                    WRITE (6, *) "#!           (angpangmokjang@hanmail.net)            !"
                    WRITE (6, *) "#!             Hanyang Univ. 2015.Apr.05.            !"
                    WRITE (6, *) "#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                    WRITE (6, *) " "
                    WRITE (6, *) "#*Ref1: T. Fukui, Y. Hatsugai, and H. Suzuki, "
                    WRITE (6, *) "#       J. J. Phys. Soc. Jap. 74, 1674 (2005),"
                    WRITE (6, *) "#     'Chern Numbers in Discretized Brillouin Zone: "
                    WRITE (6, *) "#  Efficient Method of Computing (Spin) Hall &
              &     Conductances'"
                    WRITE (6, *) "#*Ref2: R. Resta, J. Phys.: Condens. Matter. 12, R107 &
              &     (2000)"
                    WRITE (6, *) "#     'Menifestations of Berry's phase in molecules "
                    WRITE (6, *) "#      and condensed matter'"
                    WRITE (6, *) "#*Ref3: W. Yao, D. Xiao, and Q. Niu, PRB 77, 235406 &
              &     (2008)"
                    WRITE (6, *) "#     'Valley-dependent optoelectronics from inversion"
                    WRITE (6, *) "#      symmetry breaking'"
                    WRITE (6, *) "#*Ref4: http://www.andrew.cmu.edu/user/feenstra/&
              &     wavetrans/"
                    WRITE (6, *) "#       R. M. Feenstra and M. Widom                  "
                    WRITE (6, *) "#       some routines has been adopted from WAVETRANS"
                    WRITE (6, *) "#*      routines for Wavefunction reading and G "
                    WRITE (6, *) "#       matching "
                    WRITE (6, *) " "
                    WRITE (6, *) "#*Syntax:"
                    WRITE (6, *) "#       berry -f file -kx nx -ky ny -s 2 -ii ni -if nf"
                    WRITE (6, *) "#*  (or)berry -f file -kx nx -ky ny -s 2 -is n "
                    WRITE (6, *) " "
                    WRITE (6, *) "#*For the detailed help: ./berry -h"
                    WRITE (6, *) " "
                    WRITE (6, *) " "
                END SUBROUTINE creditinfo

                SUBROUTINE write_special_kpoint(b1, b2, b3)
                    IMPLICIT REAL * 8(a - h, o - z)
                    DIMENSION SKP(3), SP(3), b1(3), b2(3), b3(3)
                    OPEN (41, file='SKP.dat', status='unknown')
                    DO i = 1, 3; SKP(i) = (0.*b1(i)) + &
 (0.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # G"

                    DO i = 1, 3; SKP(i) = (2./3.*b1(i)) + &
 (1./3.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # K1"
                    DO i = 1, 3; SKP(i) = (-1./3.*b1(i)) + &
 (1./3.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # K1"

                    DO i = 1, 3; SKP(i) = (1./3.*b1(i)) + &
 (2./3.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # K2"
                    DO i = 1, 3; SKP(i) = (1./3.*b1(i)) + &
 (-1./3.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # K2"

                    DO i = 1, 3; SKP(i) = (1./2.*b1(i)) + &
 (0.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # M1"

                    DO i = 1, 3; SKP(i) = (1./2.*b1(i)) + &
 (1./2.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # M2"

                    DO i = 1, 3; SKP(i) = (0.*b1(i)) + &
 (1./2.*b2(i)) + &
 (0.*b3(i)); END DO
                    WRITE (41, '(3F12.6,A)') (SKP(i), i=1, 3), "  # M3"

                    CLOSE (41)
                END SUBROUTINE write_special_kpoint

                SUBROUTINE extendingBZ(fbz, ixt)
                    IMPLICIT REAL * 8(a - h, o - z)
                    CHARACTER * 75 fbz
                    CHARACTER * 200 A, S, P
                    DIMENSION b1(3), b2(3), b3(3), xb(3)
                    REAL*8, ALLOCATABLE :: xrecivec(:, :), xrecilat(:, :), xdata_(:)
                    REAL*8, ALLOCATABLE :: recivec(:, :), recilat(:, :), data_(:)
                    INTEGER IOstatus

                    IOstatus = 0; II = 1; ik = 0
                    OPEN (51, file=fbz, status='unknown')
                    OPEN (61, file='EXT.dat', status='unknown')
                    READ (51, '(A)', IOSTAT=IOstatus) A
                    S = A(1:1)
                    DO WHILE (TRIM(S) .EQ. '#' .OR. II .EQ. 35)
                        WRITE (6, '(A)') TRIM(A)
                        WRITE (61, '(A)') TRIM(A)
                        READ (51, '(A)', IOSTAT=IOstatus) A
                        IF (TRIM(A(3:9)) .EQ. 'NKPOINT') THEN
                            P = A(23:27); READ (P, *) nk
                        ELSE IF (TRIM(A(3:12)) .EQ. 'RECIVEC B1') THEN
                            P = A(26:34); READ (P, *) b1(1)
                            P = A(39:47); READ (P, *) b1(2)
                            P = A(52:60); READ (P, *) b1(3)
                        ELSE IF (TRIM(A(3:12)) .EQ. 'RECIVEC B2') THEN
                            P = A(26:34); READ (P, *) b2(1)
                            P = A(39:47); READ (P, *) b2(2)
                            P = A(52:60); READ (P, *) b2(3)
                        ELSE IF (TRIM(A(3:12)) .EQ. 'RECIVEC B3') THEN
                            P = A(26:34); READ (P, *) b3(1)
                            P = A(39:47); READ (P, *) b3(2)
                            P = A(52:60); READ (P, *) b3(3)
                        END IF
                        S = TRIM(A(1:1))
                        II = ICHAR(TRIM(A))
                    END DO
                    ALLOCATE (recivec(3, nk))
                    ALLOCATE (recilat(3, nk))
                    ALLOCATE (data_(nk))
                    ALLOCATE (xrecivec(3, ixt * 2 * ixt * 2 * nk))
                    ALLOCATE (xrecilat(3, ixt * 2 * ixt * 2 * nk))
                    ALLOCATE (xdata_(ixt * 2 * ixt * 2 * nk))

                    BACKSPACE (51)
                    DO WHILE (IOstatus .EQ. 0)
                        ik = ik + 1
                        READ (51, *, IOSTAT=IOstatus) (recivec(i, ik), i=1, 3), data_(ik), &
                            (recilat(i, ik), i=1, 3)
                    END DO

                    DO ik = 1, nk
                        WRITE (6, '(3F11.6,A,F16.6,A,3F11.6)') (recivec(i, ik), i=1, 3), &
                            "     ", data_(ik), &
                            "                  ", (recilat(i, ik), i=1, 3)
                    END DO

                    kk = 0   ! extend
                    DO ib2 = -1 * (ixt - 1) + 1, ixt
                        DO ib1 = -1 * (ixt - 1) + 1, ixt
                            DO ik = 1, nk
                                kk = kk + 1
                                xrecivec(1, kk) = recivec(1, ik) + (ib1 - 1) * b1(1) + (ib2 - 1) * b2(1)
                                xrecivec(2, kk) = recivec(2, ik) + (ib1 - 1) * b1(2) + (ib2 - 1) * b2(2)
                                xrecivec(3, kk) = recivec(3, ik) + (ib1 - 1) * b1(3) + (ib2 - 1) * b2(3)
                                xrecilat(1, kk) = recilat(1, ik) + (ib1 - 1)
                                xrecilat(2, kk) = recilat(2, ik) + (ib2 - 1)
                                xrecilat(3, kk) = recilat(3, ik)
                                xdata_(kk) = data_(ik)
                                kext = kk
                            END DO
                        END DO
                    END DO

                    WRITE (6, *) " "
                    WRITE (6, '(A,I1,A,I1)') "# EXTENDING DATA GRID by : ", ixt, ' x ', ixt
                    WRITE (6, '(A)') "# SORTING K-grids..."
                    DO k = ixt - 1, 1, -1  ! sorting kx
                        DO j = 1, k
                            IF (xrecivec(1, j + 1) .GT. xrecivec(1, j)) THEN
                                xb(:) = xrecivec(:, j)
                                xrecivec(:, j) = xrecivec(:, j + 1)
                                xrecivec(:, j + 1) = xb(:)
                                xb(:) = xrecilat(:, j)
                                xrecilat(:, j) = xrecilat(:, j + 1)
                                xrecilat(:, j + 1) = xb(:)
                                xtemp = xdata_(j)
                                xdata_(j) = xdata_(j + 1)
                                xdata_(j + 1) = xtemp
                            END IF
                        END DO
                    END DO
                    DO k = ixt - 1, 1, -1  ! sorting ky
                        DO j = 1, k
                            IF (xrecivec(1, j + 1) .EQ. xrecivec(1, j)) THEN
                                IF (xrecivec(2, j + 1) .GT. xrecivec(2, j)) THEN
                                    xb(:) = xrecivec(:, j)
                                    xrecivec(:, j) = xrecivec(:, j + 1)
                                    xrecivec(:, j + 1) = xb(:)
                                    xb(:) = xrecilat(:, j)
                                    xrecilat(:, j) = xrecilat(:, j + 1)
                                    xrecilat(:, j + 1) = xb(:)
                                    xtemp = xdata_(j)
                                    xdata_(j) = xdata_(j + 1)
                                    xdata_(j + 1) = xtemp
                                END IF
                            END IF
                        END DO
                    END DO

                    DO ik = 1, kext
                        WRITE (61, '(3F11.6,A,F16.6,A,3F11.6)') (xrecivec(i, ik), i=1, 3), &
                            "     ", xdata_(ik), "                  ", &
                            (xrecilat(i, ik), i=1, 3)
                    END DO
                    WRITE (6, '(A)') "# DONE! result is in 'EXT.dat' "

                    STOP
                END SUBROUTINE extendingBZ

                SUBROUTINE test
                    IMPLICIT REAL * 8(a - h, o - z)
                    COMPLEX * 8 a, b, c, d
                    CHARACTER * 75 foname
                    a = (2.2, -1.3)
                    b = (3.2, -5.4)

                    c = a + (0., 1.) * b
                    d = a - (0., 1.) * b

                    WRITE (6, *) "AAA", a
                    WRITE (6, *) "BBB", b
                    WRITE (6, *) "CCC", c, ABS(c)
                    WRITE (6, *) "DDD", d, ABS(d)

                    STOP
                END SUBROUTINE test

                SUBROUTINE help(ver_tag)
                    CHARACTER * 75 ver_tag
                    WRITE (6, *) "          **** PROGRAM INSTRUCTION ***"
                    WRITE (6, *) " "
                    WRITE (6, *) ver_tag
                    WRITE (6, *) " "
                    WRITE (6, *) "*LIMITATION : -This program is ONLY for 2D system."
                    WRITE (6, *) "            : -It tries to find the VBM for first KPT,"
                    WRITE (6, *) "            :  it may result in some problem when you"
                    WRITE (6, *) "            :  dealing with metallic states. "
                    WRITE (6, *) "            :  Be careful!"
                    WRITE (6, *) "*NOTE1 The x,y,z components of each G value are given"
                    WRITE (6, *) "       in terms of the ig values and the components "
                    WRITE (6, *) "       of the recip. lattice vectors according to:"
                    WRITE (6, *) "        ig1*b1_x + ig2*b2_x + ig3*b3_x,"
                    WRITE (6, *) "        ig1*b1_y + ig2*b2_y + ig3*b3_y, and"
                    WRITE (6, *) "        ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively,"
                    WRITE (6, *) "       with"
                    WRITE (6, *) "        ig1=ig(1,iplane),"
                    WRITE (6, *) "        ig2=ig(2,iplane), and"
                    WRITE (6, *) "        ig3=ig(3,iplane),"
                    WRITE (6, *) "       where iplane=1,2,...,nplane(k,ispin) is an "
                    WRITE (6, *) "       index incrementing the plane waves for specific"
                    WRITE (6, *) "       k and spin values"
                    WRITE (6, *) " "
                    WRITE (6, *) "*NOTE2 The energy eigenvalues are complex, as provided"
                    WRITE (6, *) "       in the WAVECAR file, but the imaginary part is "
                    WRITE (6, *) "       zero (at least for cases investigated thus far)"
                    WRITE (6, *) " "
                    WRITE (6, *) "*Syntax:berry -f file -kx nx -ky ny -s 2 -ii ni -if nf"
                    WRITE (6, *) "*   (or)berry -f file -kx nx -ky ny -s 2 -is n "
                    WRITE (6, *) " "
                    WRITE (6, *) "             ### POSSIBLE OPTIONS ###"
                    WRITE (6, *) " -f filename      : File name to be read"
                    WRITE (6, *) "                  : Default: WAVECAR"
                    WRITE (6, *) " -kx(ky) kx(ky)   : k-point grid of your system"
                    WRITE (6, *) " -s 2 or 1        : for the noncollinear case, -s 2"
                    WRITE (6, *) "                  : for the collinear or NM,   -s 1"
                    WRITE (6, *) "                  :  Default : 2 if ISPIN 1"
                    WRITE (6, *) "                  :          : 1 if ISPIN 2"
                    WRITE (6, *) " -ii(if) ni(nf)   : Once specified, berry curvature "
                    WRITE (6, *) "                  : for multiband (from ni-th to nf-th"
                    WRITE (6, *) "                  : state) states will be evaluated,"
                    WRITE (6, *) "                  : unless '-cd' is not 1."
                    WRITE (6, *) "                  :  Default -cd 0 -> -ii  1  -if VBM "
                    WRITE (6, *) "                  :          -cd 1 -> -ii VBM -if CBM"
                    WRITE (6, *) " -is n            : if specified, berry curvature for"
                    WRITE (6, *) "                  : single band (n-th band) will be "
                    WRITE (6, *) "                  : evaluated."
                    WRITE (6, *) " -o filename      : Print Berry curvature distribution"
                    WRITE (6, *) "                  : to the 'filename.dat'             "
                    WRITE (6, *) "                  :  Default : BERRYCURV.dat          "
                    WRITE (6, *) " -kp np           : Print Berry curvature distribution"
                    WRITE (6, *) "                  : with extended BZ with 'np x np'-BZ"
                    WRITE (6, *) "                  :  Default : 2              "
                    WRITE (6, *) " -skp 1(or 0)     : Specify whether special K-points "
                    WRITE (6, *) "                  : will be printed in the SKP.dat "
                    WRITE (6, *) "                  : You may specify the points by hand"
                    WRITE (6, *) "                  : into the source code, the routine,"
                    WRITE (6, *) "                  : 'write_special_kpoint'"
                    WRITE (6, *) "                  :  Default : 0"
                    WRITE (6, *) " -ne  ne          : Specify total number of electrons"
                    WRITE (6, *) "                  : Usual insulating cases, you don't"
                    WRITE (6, *) "                  : need to specify it, program will "
                    WRITE (6, *) "                  : find total number of electrons,"
                    WRITE (6, *) "                  : but for semimetallic system, it "
                    WRITE (6, *) "                  : may be more clear to put by hand"
                    WRITE (6, *) " -cd  1(or 0)     : Calculate spin- and k-resolved"
                    WRITE (6, *) "                  : degree of optical polarization,"
                    WRITE (6, *) "                  : degree of circular polarization,"
                    WRITE (6, *) "                  : between valence & conduction band"
                    WRITE (6, *) "                  : If set to '1', Berry cuvature will"
                    WRITE (6, *) "                  : not be evaluated"
                    WRITE (6, *) "                  : You may set -ii and -if together,"
                    WRITE (6, *) "                  : defining VBM & CBM, respectively"
                    WRITE (6, *) "                  :  Default : 0"
                    WRITE (6, *) " -ixt np -fbz f   : **For the special purpose"
                    WRITE (6, *) "                  : Read file 'f' and extend data to"
                    WRITE (6, *) "                  : np x np periodic field. Note that"
                    WRITE (6, *) "                  : -ixt option should be used with "
                    WRITE (6, *) "                  :'-fbz filename_of_.dat_file' option"
                    WRITE (6, *) " -vel 1 -is n     : **For the special purpose"
                    WRITE (6, *) "                  : Calculate velocity expectation "
                    WRITE (6, *) "                  : vaule (v_x, v_y) of n-th state"
                    WRITE (6, *) " -z2  1           : Calculate Z2 invariant   "
                    WRITE (6, *) "                  : Using Fukui's method (JSPJ 76,"
                    WRITE (6, *) "                  : 053702 2007) which is lattice "
                    WRITE (6, *) "                  : version of Fu and Kane method"
                    WRITE (6, *) "                  : (PRB 74, 195312 2006) "
                    WRITE (6, *) " -hf  1           : Use half BZ, ISYM=1 or 2, "
                    WRITE (6, *) "                  : currently only works with -z2 1"
                    WRITE (6, *) " -wf nb -k nk     : **For the special purpose"
                    WRITE (6, *) "  -ng nx,ny,nz    : Calculate real-space wavefunction"
                    WRITE (6, *) "  -im 1           : output will be written in "
                    WRITE (6, *) "  -ishift rx,ry,rz: PARCHG-W-K$k-E$wf-(IM)-SPIN$ispin"
                    WRITE (6, *) "                  : nb and nk : index of band & kpoint"
                    WRITE (6, *) "                  : -ng option determins nx,ny,nz grid"
                    WRITE (6, *) "                  :  for the cube: default=nbmax"
                    WRITE (6, *) "                  : -im option: plot imaginary?"
                    WRITE (6, *) "                  : -ishift: shift origin with respect"
                    WRITE (6, *) "                  :  to the direct coord. rx,ry,rz"
                    WRITE (6, *) " "
                    WRITE (6, *) "* default: -f WAVECAR -kx 2 -ky 2 -s 2 -ii 1 -if VBM &
              &     -kp 1"
                    WRITE (6, *) "* here, VBM is valence band maximum"
                    WRITE (6, *) " "
                    WRITE (6, *) "*Compilation:  gfortran or ifort. "
                    WRITE (6, *) " Flag '-assume byterecl' is required for ifort."
                    WRITE (6, *) " for OSX,  -Wl,-stack_size,0x80000000 may be required"
                    WRITE (6, *) " ex-noMPI)ifort -fpp -assume byterecl -mkl &
              &     -o vaspberry vaspberry.f"
                    WRITE (6, *) " ex-MPI)mpif90 -DMPI_USE -mkl -fpp -assume byterecl &
              &     -o vaspberry vaspberry.f"
                    WRITE (6, *) " ex-noMPI-gfortran) gfortran -I/opt/local/include &
              &     -L/opt/local/lib/lapack/ -l lapack -o vaspberry vaspberry_gfortran&
              &     _serial.f"

!     write(6,*)" ex-MPI) mpif90 -DMPI_USE -mkl -fpp -assume byterecl -o vaspberry vaspberry.f "

                    STOP
                END SUBROUTINE help
