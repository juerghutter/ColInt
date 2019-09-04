! **************************************************************************************************
MODULE task_list_methods
   USE task_list_types,                 ONLY: task_list_type,&
                                              grid_base_type
   USE spmat_types,                     ONLY: spmat_p_type,&
                                              spmat_type

   IMPLICIT NONE

   INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)

   PRIVATE

   PUBLIC :: int2pair, read_task_file, print_task_info

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
   SUBROUTINE int2pair(res, ilevel, image, iatom, jatom, iset, jset, ipgf, jpgf, &
                       nimages, natom, maxset, maxpgf)
      INTEGER(KIND=int_8), INTENT(IN)                    :: res
      INTEGER, INTENT(OUT)                               :: ilevel, image, iatom, jatom, iset, jset, &
                                                            ipgf, jpgf
      INTEGER, INTENT(IN)                                :: nimages, natom, maxset, maxpgf

      INTEGER(KIND=int_8) :: iatom8, ijatom, ijset, img, ipgf8, iset8, jatom8, jpgf8, jset8, &
         maxpgf8, maxset8, natom8, nimages8, nlev1, nlev2, nlev3, nlev4, tmp

      natom8 = natom; maxset8 = maxset; maxpgf8 = maxpgf; nimages8 = nimages
      !
      nlev1 = maxpgf8**2
      nlev2 = maxset8**2*nlev1
      nlev3 = natom8**2*nlev2
      nlev4 = nimages8*nlev3
      !
      ilevel = INT(res/nlev4)
      tmp = MOD(res, nlev4)
      img = tmp/nlev3+1
      tmp = MOD(tmp, nlev3)
      ijatom = tmp/nlev2
      iatom8 = ijatom/natom8+1
      jatom8 = MOD(ijatom, natom8)+1
      tmp = MOD(tmp, nlev2)
      ijset = tmp/nlev1
      iset8 = ijset/maxset8+1
      jset8 = MOD(ijset, maxset8)+1
      tmp = MOD(tmp, nlev1)
      ipgf8 = tmp/maxpgf8+1
      jpgf8 = MOD(tmp, maxpgf8)+1
      !
      image = INT(img)
      iatom = INT(iatom8); jatom = INT(jatom8); iset = INT(iset8); jset = INT(jset8)
      ipgf = INT(ipgf8); jpgf = INT(jpgf8)

   END SUBROUTINE int2pair

! **************************************************************************************************

   SUBROUTINE read_task_file(task_list, spmats, rs_rho, kind_of, posat)
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(spmat_p_type), DIMENSION(:), POINTER          :: spmats
      TYPE(grid_base_type), DIMENSION(:), POINTER        :: rs_rho
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      REAL(KIND=dp), DIMENSION(:,:), POINTER             :: posat

      CHARACTER(LEN=50)                                  :: inline
      INTEGER                                            :: iu, ng, irs, ir, ig, n1, n2, nset
      INTEGER, DIMENSION(2, 3)                           :: gridbounds
      TYPE(spmat_type), POINTER                          :: spmat

      iu = 5
      OPEN (unit=iu, position="rewind")
      !
      READ (iu, *) inline !"S> TASK GENERAL"
      READ (iu, *) task_list%ntasks
      READ (iu, *) task_list%ngrid_levels
      READ (iu, *) task_list%nimages
      READ (iu, *) task_list%natoms
      READ (iu, *) task_list%maxset
      READ (iu, *) task_list%maxpgf
      READ (iu, *) task_list%lmax_all
      READ (iu, *) task_list%maxco
      READ (iu, *) task_list%maxsgf_set
      READ (iu, *) task_list%eps_rho_rspace
      READ (iu, *) inline !"E> TASK GENERAL"
      !
      READ (iu, *) inline !"S> TASK GRID"
      READ (iu, *) ng
      ALLOCATE(task_list%cinfo(ng))
      DO ig=1,ng
         READ(iu, *) task_list%cinfo(ig)%dh
         READ(iu, *) task_list%cinfo(ig)%dh_inv
         READ(iu, *) task_list%cinfo(ig)%perd
         READ(iu, *) task_list%cinfo(ig)%npts
         READ(iu, *) task_list%cinfo(ig)%lb_grid
         READ(iu, *) task_list%cinfo(ig)%drmin
         READ(iu, *) task_list%cinfo(ig)%max_radius
         READ(iu, *) task_list%cinfo(ig)%max_rad_ga
         READ(iu, *) n1, n2
         IF(n1*n2 > 0) THEN
            ALLOCATE(task_list%cinfo(ig)%lb_cube(n1,n2))
            READ(iu, *) task_list%cinfo(ig)%lb_cube
         ENDIF
         READ(iu, *) n1, n2
         IF(n1*n2 > 0) THEN
            ALLOCATE(task_list%cinfo(ig)%ub_cube(n1,n2))
            READ(iu, *) task_list%cinfo(ig)%ub_cube
         ENDIF
         READ(iu, *) n1, n2
         IF(n1*n2 > 0) THEN
            ALLOCATE(task_list%cinfo(ig)%sphere_bounds(n1,n2))
            READ(iu, *) task_list%cinfo(ig)%sphere_bounds
         ENDIF
      END DO
      READ (iu, *) inline !"E> TASK GRID"
      !
      READ (iu, *) inline !"S> TASK BASIS"
      READ (iu, *) ng
      ALLOCATE(task_list%basis(ng))
      DO ig=1,ng
         READ (iu, *) nset
         IF(nset > 0) THEN
            ALLOCATE(task_list%basis(ig)%first_sgf(nset))
            READ (iu, *) task_list%basis(ig)%first_sgf
            ALLOCATE(task_list%basis(ig)%lmax(nset))
            READ (iu, *) task_list%basis(ig)%lmax
            ALLOCATE(task_list%basis(ig)%lmin(nset))
            READ (iu, *) task_list%basis(ig)%lmin
            ALLOCATE(task_list%basis(ig)%npgf(nset))
            READ (iu, *) task_list%basis(ig)%npgf
            ALLOCATE(task_list%basis(ig)%nsgf(nset))
            READ (iu, *) task_list%basis(ig)%nsgf
            READ(iu, *) n1, n2
            ALLOCATE(task_list%basis(ig)%sphi(n1,n2))
            READ(iu, *) task_list%basis(ig)%sphi
            READ(iu, *) n1, n2
            ALLOCATE(task_list%basis(ig)%zet(n1,n2))
            READ(iu, *) task_list%basis(ig)%zet
         END IF
      END DO
      READ (iu, *) inline !"E> TASK BASIS"
      !
      READ (iu, *) inline !"S> TASK DIST"
      READ(iu, *) n1, n2
      ALLOCATE(task_list%dist_ab(n1,n2))
      READ(iu, *) task_list%dist_ab
      READ(iu, *) n1, n2
      ALLOCATE(task_list%taskstart(n1,n2))
      ALLOCATE(task_list%taskstop(n1,n2))
      READ(iu, *) task_list%taskstart
      READ(iu, *) task_list%taskstop
      READ(iu, *) n1
      ALLOCATE(task_list%npairs(n1))
      READ(iu, *) task_list%npairs
      READ (iu, *) inline !"E> TASK DIST"
      !
      READ (iu, *) inline !"S> TASK LIST"
      READ(iu, *) n1, n2
      ALLOCATE(task_list%tasks(n1,n2))
      READ(iu, *) task_list%tasks
      READ (iu, *) inline !"E> TASK LIST"
      !
      READ (iu, *) inline !"S> SPMAT"
      READ (iu, *) ng
      ALLOCATE(spmats(ng))
      DO ig = 1, ng
         spmat => spmats(ig)%spmat
         READ (iu, *) spmat%nrow, spmat%ncol, spmat%nnz
         ALLOCATE(spmat%rind(spmat%nrow))
         DO ir = 1, spmat%nrow
            READ (iu, *) irs
            IF (irs > 0) THEN
               ALLOCATE(spmat%rind(ir)%col(irs))
               ALLOCATE(spmat%rind(ir)%blk(irs))
               READ (iu, *) spmat%rind(ir)%col(1:irs)
               READ (iu, *) spmat%rind(ir)%blk(1:irs)
            ELSE
               READ (iu, *)
               READ (iu, *)
            END IF
         END DO
         !data
         ALLOCATE(spmat%spdata(spmat%nnz))
         DO ir = 1, spmat%nnz
            READ (iu, *) n1, n2
            ALLOCATE(spmat%spdata(ir)%mat(n1,n2))
            READ (iu, *) spmat%spdata(ir)%mat
         END DO
      END DO
      READ (iu, *) inline !"E> SPMAT"
      !
      READ (iu, *) inline !"S> ATOMS"
      READ (iu, *) ng
      NULLIFY(kind_of, posat)
      ALLOCATE(kind_of(ng))
      ALLOCATE(posat(3,ng))
      READ (iu, *) kind_of
      READ (iu, *) posat
      READ (iu, *) inline !"E> ATOMS"
      !
      READ (iu, *) inline !"S> RS GRID"
      READ (iu, *) ng
      ALLOCATE(rs_rho(ng))
      DO ig = 1, ng
         READ (iu, *) gridbounds
         ALLOCATE(rs_rho(ig)%grid(gridbounds(1, 1):gridbounds(2, 1),&
                                  gridbounds(1, 2):gridbounds(2, 2),&
                                  gridbounds(1, 3):gridbounds(2, 3)))
      END DO
      READ (iu, *) !"E> RS GRID"
      !
      CLOSE (unit=iu)

   END SUBROUTINE read_task_file

! **************************************************************************************************

   SUBROUTINE print_task_info(task_list, spmats, rs_rho, kind_of, posat, io)
      TYPE(task_list_type), INTENT(IN)                   :: task_list
      TYPE(spmat_p_type), DIMENSION(:), POINTER          :: spmats
      TYPE(grid_base_type), DIMENSION(:), INTENT(IN)     :: rs_rho
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)          :: posat
      INTEGER, INTENT(IN)                                :: io

      INTEGER                                            :: ng, ig, ip
      INTEGER, DIMENSION(2, 3)                           :: gridbounds
      REAL(KIND=dp), DIMENSION(:,:,:), POINTER           :: grid
      TYPE(spmat_type), POINTER                          :: spmat

      WRITE (UNIT=io, FMT="(T1,A)") "> Tasklist information summary"
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Number of atoms",task_list%natoms
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Number of cell images (k-points)",task_list%nimages
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Number of grid levels",task_list%ngrid_levels
      ng = SIZE(task_list%basis)
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Number of atom types",ng
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Max. l-qn in basis",task_list%lmax_all
      ng = task_list%ngrid_levels
      DO ig=1,ng
         WRITE (UNIT=io, FMT="(T5,A,T70,I10)") "> Grid information",ig
         WRITE (UNIT=io, FMT="(T7,A,T50,3F10.6)") "> Grid box matrix  A",task_list%cinfo(ig)%dh(:,1)
         WRITE (UNIT=io, FMT="(T7,A,T50,3F10.6)") ">                  B",task_list%cinfo(ig)%dh(:,2)
         WRITE (UNIT=io, FMT="(T7,A,T50,3F10.6)") ">                  C",task_list%cinfo(ig)%dh(:,3)
         WRITE (UNIT=io, FMT="(T7,A,T50,3I10)") "> Grid periodicity  ",task_list%cinfo(ig)%perd
         WRITE (UNIT=io, FMT="(T7,A,T50,3I10)") "> Number of Grid pointsi",task_list%cinfo(ig)%npts
         WRITE (UNIT=io, FMT="(T7,A,T70,I10)") "> Max function radius",task_list%cinfo(ig)%max_radius
         grid => rs_rho(ig)%grid(:,:,:)
         WRITE (UNIT=io, FMT="(T7,A,T60,2I10)") "> Realspace Grid bounds X",LBOUND(grid,1),UBOUND(grid,1)
         WRITE (UNIT=io, FMT="(T7,A,T60,2I10)") ">                       Y",LBOUND(grid,2),UBOUND(grid,2)
         WRITE (UNIT=io, FMT="(T7,A,T60,2I10)") ">                       Z",LBOUND(grid,3),UBOUND(grid,3)
      END DO
      WRITE (UNIT=io, FMT="(T3,A,T70,I10)") "> Total number of tasks",task_list%ntasks
      ng = task_list%ngrid_levels
      DO ig=1,ng
         WRITE (UNIT=io, FMT="(T5,A,T70,I10)") "> Number of atom pairs",task_list%npairs(ig)
         ip=task_list%taskstop(task_list%npairs(ig),ig) - task_list%taskstart(1,ig) + 1
         WRITE (UNIT=io, FMT="(T5,A,T70,I10)") "> Number of tasks",ip
      END DO
      WRITE (UNIT=io, FMT="(T3,A)") "> Sparse matrix info"
      ng = SIZE(spmats)
      WRITE (UNIT=io, FMT="(T5,A,T70,I10)") "> Number of sparse matrices ",ng
      DO ig = 1, ng
         spmat => spmats(ig)%spmat
         WRITE (UNIT=io, FMT="(T7,A,T70,I10)") "> Sparse matrix",ig
         WRITE (UNIT=io, FMT="(T7,A,T60,2I10)") "> Number of rows and columns",spmat%nrow, spmat%ncol
         WRITE (UNIT=io, FMT="(T7,A,T70,I10)") "> Number of non-zero blocks stored",spmat%nnz
      END DO

   END SUBROUTINE print_task_info

END MODULE task_list_methods
