PROGRAM colint
   USE spmat_types,                     ONLY: spmat_p_type
   USE task_list_types,                 ONLY: task_list_type,&
                                              grid_base_type,&
                                              allocate_task_list,&
                                              deallocate_task_list
   USE task_list_methods,               ONLY: read_task_file,&
                                              print_task_info
   USE grid_base_ref,                   ONLY: grid_collocate_ref

   IMPLICIT NONE

   INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)

   INTEGER, PARAMETER                          :: io=6
   INTEGER                                     :: lmax
   INTEGER, DIMENSION(:), POINTER              :: kind_of
   REAL                                        :: t1, t2
   REAL(KIND=dp), DIMENSION(:,:), POINTER      :: posat
   TYPE(task_list_type), POINTER               :: task_list
   TYPE(spmat_p_type), DIMENSION(:), POINTER   :: spmats
   TYPE(grid_base_type), DIMENSION(:), POINTER :: rs_rho

   WRITE (UNIT=io, FMT="(T1,A)") REPEAT("=", 80)
   WRITE (UNIT=io, FMT="(T1,A,T34,A,T80,A)") "=","ColInt MiniApp","="
   WRITE (UNIT=io, FMT="(T1,A,T23,A,T80,A)") "=","GPW Collocate and Integrate Routines","="
   WRITE (UNIT=io, FMT="(T1,A)") REPEAT("=", 80)
   CALL allocate_task_list(task_list)
   NULLIFY(spmats, rs_rho)

   WRITE (UNIT=io, FMT="(/,T1,A)") "> Start reading tasklist information"
   CALL read_task_file(task_list, spmats, rs_rho, kind_of, posat)
   WRITE (UNIT=io, FMT="(T1,A)") "> End reading tasklist information"
   CALL print_task_info(task_list, spmats, rs_rho, kind_of, posat, io)

   WRITE (UNIT=io, FMT="(/,T1,A)") "> Start reference implementation"
   CALL cpu_time(t1)
   CALL set_zero(rs_rho)
   lmax = task_list%lmax_all
   CALL grid_collocate_ref(rs_rho, spmats, task_list, .FALSE., kind_of, posat, lmax)
   CALL rs_grid_info(rs_rho, task_list, io)
   CALL set_zero(rs_rho)
   lmax = task_list%lmax_all + 1
   CALL grid_collocate_ref(rs_rho, spmats, task_list, .TRUE., kind_of, posat, lmax)
   CALL rs_grid_info(rs_rho, task_list, io)
   CALL cpu_time(t2)
   WRITE(UNIT=io, FMT="(T3,A,T60,F20.7)") "CPU TIME (seconds)",t2-t1
   WRITE (UNIT=io, FMT="(T1,A)") "> END reference implementation"

   CALL deallocate_task_list(task_list)

   WRITE (UNIT=io, FMT="(/,T1,A)") REPEAT("=", 80)
   WRITE (UNIT=io, FMT="(T1,A,T34,A,T80,A)") "=","ColInt MiniApp","="
   WRITE (UNIT=io, FMT="(T1,A)") REPEAT("=", 80)

   CONTAINS

   SUBROUTINE set_zero(rs_rho)
      TYPE(grid_base_type), DIMENSION(:), POINTER :: rs_rho
      INTEGER                                     :: ig
      DO ig=1,SIZE(rs_rho)
         rs_rho(ig)%grid = 0.0_dp
      END DO
   END SUBROUTINE set_zero

   SUBROUTINE rs_grid_info(rs_rho, task_list, io)
      TYPE(grid_base_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(task_list_type), POINTER               :: task_list
      INTEGER, INTENT(IN)                         :: io
      INTEGER                                     :: ig
      REAL(KIND=dp)                               :: dvol, dsum, dall
      REAL(KIND=dp), DIMENSION(3,3)               :: dh

      WRITE(UNIT=io, FMT="(T3,A,T60,F20.7)") "> Grid values information"
      dall = 0.0_dp
      DO ig=1,SIZE(rs_rho)
         dh = task_list%cinfo(ig)%dh
         dvol = dh(1, 1)*(dh(2, 2)*dh(3, 3)-dh(2, 3)*dh(3, 2))+ &
                dh(1, 2)*(dh(2, 3)*dh(3, 1)-dh(2, 1)*dh(3, 3))+ &
                dh(1, 3)*(dh(2, 1)*dh(3, 2)-dh(2, 2)*dh(3, 1))
         dvol = ABS(dvol)
         dsum = dvol*SUM(rs_rho(ig)%grid)
         WRITE(UNIT=io, FMT="(T5,A,I3,T60,F20.14)") "> Grid integration ",ig,dsum
         dall = dall + dsum
      END DO
      WRITE(UNIT=io, FMT="(T3,A,T60,F20.14)") "> Total Grid integration ",dall

   END SUBROUTINE rs_grid_info

END PROGRAM
