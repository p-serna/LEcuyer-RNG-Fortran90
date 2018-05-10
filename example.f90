!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Example
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
 use rng
 implicit none
!!!!!!!!!!!!! RNG  
 integer(KIND=8)  :: iseed(6)       
 real(KIND=8)     ::RngStream_RandU01,U01d
 character(LEN=2) :: cname
 type(RngStream),allocatable :: gar(:)
 type(RngStream)  :: RngStream_CreateStream,g
!!!!!!!!!!!!! Openmp
 integer*4 :: nth,ith
 integer   :: iproc,omp_get_num_threads,omp_get_thread_num
 integer   :: omp_get_num_procs,omp_get_max_threads
 real*8      :: omp_get_wtime      

!!!!!!!!!!!! This is to initialize openmp: number of threads
 nth = omp_get_num_procs()
 call omp_set_num_threads(nth)
!!!!!! I think this call to env variable:$OMP_NUM_THREADS 


! Seed for the random number generator
!   ISEED MAX
!      iseed = (/ 4294967087, 4294967087, 4294967087, 
! 4294944443,  4294944443,  4294944443 /)

 iseed = (/ 0747946432, 0294967087, 1294967087, &
      00934464363, 2141353523,  1111241368 /)

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialization with seed
 call RngStream_SetPackageSeed(iseed) 
 
! allocating nth threads   
 allocate(gar(nth))
 cname = "gp"
! Main thread
 g = RngStream_CreateStream(cname)
 do ith=1,nth
  write(cname,'(A,I1)') 'g',ith-1        
! ith thread
  gar(ith) = RngStream_CreateStream (cname)
 enddo       
!!!!!!!!!!!! RNG


 write(*,*) RngStream_RandU01(g)


!$omp parallel DEFAULT(shared) private(ith)
 ith=omp_get_thread_num()+1 
 write(*,*) ith,RngStream_RandU01(gar(ith))
!$omp end parallel         



end


