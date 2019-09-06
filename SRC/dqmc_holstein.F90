module DQMC_Hubbard
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_CFG
  use DQMC_PHY0
  use DQMC_PHY2
#ifdef DQMC_CKB
  use DQMC_CheckerBoard
#else
  use DQMC_MATB
#endif
  use DQMC_SEQB
  use DQMC_GFUN
  use DQMC_STRUCT
  
  implicit none 
  ! 
  ! This module contains the data type and subroutines for 
  ! computing DQMC.  This should be the only module that user
  ! program need to include.
  !
  ! The data type is consisted of four parts
  !    1. Parameters of Hubbard's model and Green's function.
  !    2. Parameters for Monte Carlo algorithm.
  !    3. Physical measurements.
  !    4. Working space.
  !
  ! There are only four subroutines for user program to call
  !    1. DQMC_Readin  : read input
  !    2. DQMC_Run     : execute DQMC
  !    3. DQMC_Dump    : write output
  !    4. DQMC_Current_Config : output current Hubbard-Stratonovich
  !                             configuration.
  ! 
  ! References
  ! ==========
  !    [1] Z. Bai, W.Chen, R. Scalettar, I. Yamazaki, "Lecture Notes 
  !        on Advances of Numerical Methods for Hubbard Quantum Monte
  !        Carlo Simulation." 
  !
  ! List of subroutines
  ! ===================
  !    DQMC_Default(Hub) : set the default value of the data type
  !    DQMC_Readin(Hub, IPT, OPT, ReadStruct) : read in parameters
  !    DQMC_Init(Hub) : Initialize the data type.
  !    DQMC_Dump(Hub, OPT) : output the parameters.
  !    DQMC_Sweep(Hub, nMeas0, v1, v2) : Metropolis algorithm.
  !    DQMC_Run(Hub) : the main subroutine of DQMC.  
  !
  !
  ! Data Type
  ! =========
  !
  type Hubbard
     ! Part 1: Parameters of Hubbard's model and Green's function
     ! ==========================================================

     ! Parameters for problem size
     integer           :: n               ! Number of sites
     integer           :: L               ! Number of time slices

     ! Parameters for Hubbard model
     integer           :: n_U
     real(wp), pointer :: U(:)            ! Hubbard U, or el-ph coupling in the case of Holstein model
     integer           :: n_t             ! Number of hoppings 
     real(wp), pointer :: t_up(:)         ! spin-up hopping 
     real(wp), pointer :: t_dn(:)         ! spin-down hopping      
     integer           :: n_mu          
     real(wp), pointer :: mu_up(:)        ! spin-up chemical potential      
     real(wp), pointer :: mu_dn(:)        ! spin-down chemical potential 
     real(wp)          :: dtau            ! imaginary time step size
     integer,  pointer :: HSF (:,:)       ! discrete Hubbard-Stratonovich Fields
     real(wp), pointer :: CHSF (:,:)      ! continuous Hubbard-Stratonovich Fields
     integer           :: HSFtype         ! flag for HSF type
     logical           :: outputHSF       ! flag for output HSF
     logical           :: continuous      ! flag for continuous HSF
     real(wp)          :: delta1          ! continuous field local move step size
     real(wp)          :: delta2          ! continuous field global move step size
     real(wp), pointer :: lambda(:)       ! parameter for contunuous HSF 
     integer           :: n_start, n_end 
     
     real(wp)          :: omega           ! phonon dispersion
     
     
     ! Underline structure  
     type(Struct)      :: S               ! Lattice structure
     

     ! For Green function computation
     type(matB)        :: B_up            ! 
     type(seqB)        :: SB_up           ! Sequential Bs 
     type(matB)        :: B_dn            ! 
     type(seqB)        :: SB_dn           ! Sequential Bs 
     type(G_fun)       :: G_up            ! Green's fun for spin up
     type(G_fun)       :: G_dn            ! Green's fun for spin down
     real(wp), pointer :: V_up(:,:)  
     real(wp), pointer :: V_dn(:,:)   
     
     ! Parameters for random number
     integer           :: idum            ! random seed for ran2
     integer           :: seed(4)         ! random seed for ran1

     ! Auxiliary variables
     real(wp), pointer :: explook(:,:)    ! Lookup table for computing V
     logical           :: comp_dn         ! indicator for wheather computing
                                          ! G_dn or not
     logical           :: neg_u           ! are all U_i < 0 ?
                                          ! CAVEAT: mixed sign may not work

     ! Part 2: Parameters for Monte Carlo algorithm
     ! ============================================
     integer           :: nWarm           ! Number of warm up step
     integer           :: nPass           ! Number of measurement step
     integer           :: nTry            ! Number of global move
     real(wp)          :: gamma           ! Parameters for Metopolis alg
                      
     integer           :: nAccept         ! The following parameters  
     integer           :: nReject         ! are used to dynamically
                                          ! adjust gamma.
                      
     integer           :: nAcceptGlobal   ! global move acceptance
     integer           :: nRejectGlobal
     integer           :: nAcceptGlobal2  ! global move 2 acceptance
     integer           :: nRejectGlobal2
     ! Part 3: Physical measurements
     ! =============================
     type(Phy0)        :: P0              ! Meas0
     type(Phy2)        :: P2              ! MeasPair
     integer           :: nMeas           
     integer           :: tausk           ! Frequency of performing Phy0 measurement 
     logical           :: meas2

     ! Part 4: Working space
     ! =============================
     type(Wspace)      :: WS

     ! Part 5: file units
     ! =============================
     integer           :: OUT_UNIT

  end type Hubbard

  integer, parameter   :: NO_MEAS0     = -1

  ! HSF parameter
  integer, parameter   :: HSF_OUTPUT_UNIT = 28
  integer, parameter   :: HSF_INPUT_UNIT  = 27
  integer, parameter   :: HSF_RANDOM_GEN  = -1 ! generate HS fields randomly
  integer, parameter   :: HSF_FROM_FILE   =  1 ! load HS fields from file
  integer, parameter   :: HSF_RESTORE     =  2 ! load HS field and RNG state from file
  integer, parameter   :: HSF_FROM_MEMORY =  0

  integer, parameter   :: HSF_DISC  =  0
  integer, parameter   :: HSF_CONT  =  1

  integer              :: SimType              ! model selection flag.
                                               !   SimType == 0: Hubbard model
                                               !   SimType == 1: Holstein model
  real(wp)             :: CHSF0                ! phonon field overall scale
  real(wp)             :: norm_phonon(2)       ! normalization factor of phonon terms
                                               !   1 => kinetic energy term
                                               !   2 => potential energy term
contains

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Config(Hub, cfg)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subrotine initializes Hubbard model from the configuration.
    !
    !
    ! Pre-assumption
    ! ==============
    !    DQMC_default should be called before this.
    !    Geometry information should be iniitialized before calling this. 
    !
    ! 
    ! Arguments
    ! =========
    !
    type(config), intent(inout)  :: cfg
    type(Hubbard), intent(inout) :: Hub                   ! Hubbard model

    ! ... Local Variables ...
    integer :: n_t, n_U, n_mu, L, HSF, nWarm, nPass
    integer :: accept, reject, HSFtype, fixw, tdm
    integer :: seed, nOrth, nWrap, nTry, nBin, nMeas, ntausk, ssxx
    character(len=slen) :: HSF_ipt, HSF_opt
    logical :: valid
    real(wp), pointer :: t_up(:)  => null()
    real(wp), pointer :: t_dn(:)  => null()  
    real(wp), pointer :: U(:)     => null()   
    real(wp), pointer :: mu_up(:) => null() 
    real(wp), pointer :: mu_dn(:) => null() 
    real(wp) :: dtau, errrate, difflim, gamma, delta1, delta2, omega 

    ! ... Executable ...

    !Open output file
    !call CFG_Get(cfg, "gfile", fname)
    !outname=trim(adjustl(fname))//".out"
    !call DQMC_open_file(outname, 'unknown', Hub%OUT_UNIT)

    ! integer parameters
    call CFG_Get(cfg, "HSF",      HSF)
    call CFG_Get(cfg, "L",        L)
    call CFG_Get(cfg, "nwarm",    nWarm)
    call CFG_Get(cfg, "npass",    nPass)
    call CFG_Get(cfg, "tdm",      tdm)
    call CFG_Get(cfg, "nbin",     nBin)
    call CFG_Get(cfg, "ntry",     nTry)
    call CFG_Get(cfg, "seed",     seed)
    call CFG_Get(cfg, "nwrap",    nWrap)
    call CFG_Get(cfg, "north",    nOrth)
    call CFG_Get(cfg, "gamma",    gamma)
    call CFG_Get(cfg, "accept",   accept)
    call CFG_Get(cfg, "reject",   reject)
    call CFG_Get(cfg, "HSFtype",  HSFtype)
    call CFG_Get(cfg, "delta1",   delta1)
    call CFG_Get(cfg, "delta2",   delta2)
    call CFG_Get(cfg, "ssxx",     ssxx)
    call CFG_Get(cfg, "fixwrap",  fixw)
    call CFG_Get(cfg, "tausk",    ntausk)
    call CFG_Get(cfg, "SimType",  SimType)
    call CFG_Get(cfg, "omega",    omega)

    ! Array parameters
    call CFG_Get(cfg, "t_up",    n_t,  t_up)
    call CFG_Get(cfg, "t_dn",    n_t,  t_dn)
    call CFG_Get(cfg, "U",       n_U,  U)
    call CFG_Get(cfg, "mu_up",   n_mu, mu_up)
    call CFG_Get(cfg, "mu_dn",   n_mu, mu_dn)

    ! Real parameters
    call CFG_Get(cfg, "dtau",    dtau)
    call CFG_Get(cfg, "difflim", difflim)
    call CFG_Get(cfg, "errrate", errrate)

    !Change nbin to 1 if we are using more than 1 CPU.
    !Results collected on each CPU will be used as bins.
    if (qmc_sim%aggr_size > 1) then
       nBin = 1
       call CFG_set(cfg, "nbin",  nBin)
    endif

    if ( SimType .eq. Holstein_model .and. HSFtype .eq. 0 ) then
      write(*,*)
      write(*,"(' Input parameter error:')")
      write(*,"('   SimType =', i3)") SimType 
      write(*,"('   HSFtype =', i3)") HSFtype
      call DQMC_Error(" Must set HSFtype=1 in Holstein mode (SimType=1).", 1)
    end if

    if (HSF == HSF_FROM_FILE .or. HSF == HSF_RESTORE) then
       ! open input file
       if (DQMC_Config_isSet(cfg, "HSFin")) then
          call CFG_Get(cfg, "HSFin", HSF_ipt)
          inquire(FILE=trim(HSF_ipt), EXIST=valid)
          if (valid) then
             open(HSF_INPUT_UNIT, FILE = trim(HSF_ipt))
          else
             call DQMC_Warning("HSF input file does not exist.", 1)
             HSF = HSF_RANDOM_GEN
          end if
       end if
    elseif (HSF /= HSF_FROM_MEMORY .and. HSF /= HSF_RANDOM_GEN) then
       call DQMC_Warning("Invalid HSF input: Use default", HSF)
       HSF = HSF_RANDOM_GEN
    end if
    
    ! open output file
    Hub%outputHSF =DQMC_Config_isSet(cfg, "HSFout")
    if (Hub%outputHSF) then
       call CFG_Get(cfg, "HSFout", HSF_opt)
       open(HSF_OUTPUT_UNIT, FILE = trim(HSF_opt))
    end if

    nmeas = 1
    ! Deactivate measurements during sweep if tdm is on
    if (tdm > 0) nmeas = 0

    ! call the function
    call DQMC_Hub_Init(Hub, U, t_up, t_dn, mu_up, mu_dn, L, n_t, n_U, n_mu, dtau, &
       HSF, nWarm, nPass, nMeas, nTry, nBin, ntausk, seed, nOrth, nWrap, fixw, &
       errrate, difflim, gamma, accept, reject, delta1, delta2, ssxx, HSFtype, omega)
    
    call CFG_Set(cfg, "n", Hub%n)

    deallocate(t_up, t_dn, mu_up, mu_dn, U)

  end subroutine DQMC_Hub_Config

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Init(Hub, U, t_up, t_dn, mu_up, mu_dn, L, n_t, n_U, n_mu, dtau, &
       HSF_IPT, nWarm, nPass, nMeas, nTry, nBin, ntausk, seed, nOrth, nWrap, fixw, &
       errrate, difflim, gamma, accept, reject, delta1, delta2, ssxx, HSFtype, omega)
    !
    ! Purpose
    ! =======
    !    This subrotine initializes Hubbard model.
    !
    ! Pre-assumption
    ! ==============
    !    DQMC_default should be called before this.
    !    Geometry information should be iniitialized before calling this. 
    !
    ! 
    ! Arguments
    ! =========
    !
    use dqmc_mpi
#ifdef __INTEL_COMPILER
  use IFPORT, only : getpid
#elif __PGI
  integer :: getpid  
#endif

#   ifdef _QMC_MPI
#      define SIMPLE_SPRNG
#      define USE_MPI
#      include "sprng_f.h"
#   endif
    type(Hubbard), intent(inout) :: Hub             ! Hubbard model
    real(wp), intent(in)  :: U(:), t_up(:), t_dn(:) 
    real(wp), intent(in)  :: mu_up(:), mu_dn(:), dtau  ! Parameters
    integer,  intent(in)  :: L, n_t, n_U, n_mu
    integer,  intent(in)  :: HSF_IPT, seed, ssxx
    integer,  intent(in)  :: nWarm, nPass, nOrth, nTry, HSFtype
    integer,  intent(in)  :: nMeas, nBin, ntausk, nWrap, fixw, accept, reject
    real(wp), intent(in)  :: errrate, difflim, gamma, delta1, delta2, omega

    ! ... Local scalar ...
    integer  :: n, i, j, HSF, ilb, val(8)
    real(wp) :: temp, lambda 
    real(wp) :: Map_up(Hub%S%nSite), Map_dn(Hub%S%nSite), Uhubb(Hub%S%nSite)
    logical  :: lex
    character(80) :: msg

    ! phonon field related
    real(wp), allocatable :: CHSFtmp(:)

#   ifdef _QMC_MPI
       SPRNG_POINTER junkPtr
#   endif

    ! ... Executable ...

    if (.not. Hub%S%checklist(STRUCT_INIT)) then
       call DQMC_Error("Must initialize lattice geometry first", 0)
    end if

    Hub%n        = Hub%S%nSite
    n            = Hub%n

    Hub%L        = L
    Hub%dtau     = dtau
    Hub%delta1   = delta1
    Hub%delta2   = delta2
    Hub%omega    = omega

    allocate(CHSFtmp(n))

    ! t parameter
    if (n_t /= Hub%S%n_t) then
       if (n_t == 1 .and. Hub%S%n_t > 1) then
          ! special case for checkerboard method
          Hub%n_t  = Hub%S%n_t
          allocate(Hub%t_up(Hub%S%n_t))
          allocate(Hub%t_dn(Hub%S%n_t))
          Hub%t_up    = t_up(1)
          Hub%t_dn    = t_dn(2)
       else
          write(msg, "(a,i5, a, i5)") "Input lattice requires ", Hub%S%n_t, &
               " hoppings, but only reads ", n_t 
          call DQMC_Error(msg, 0)
       end if
    else
       Hub%n_t = n_t
       allocate(Hub%t_up(n_t))
       allocate(Hub%t_dn(n_t))
       Hub%t_up = t_up
       Hub%t_dn = t_dn
    end if

    ! U parameter
    Hub%n_U      = n_U
    allocate(Hub%U(n_U))
    Hub%U        = U
     
    ! mu parameter
    Hub%n_mu     = n_mu
    allocate(Hub%mu_up(n_mu))
    allocate(Hub%mu_dn(n_mu))

    Hub%mu_up       = mu_up
    Hub%mu_dn       = mu_dn
    do i = 1, n
       Map_up(i) = mu_up(Hub%S%Map(i))
       Map_dn(i) = mu_dn(Hub%S%Map(i))
       Uhubb(i)  = abs(Hub%U(Hub%S%Map(i)))
    end do

    if (SimType .eq. Hubbard_model) then

      Hub%comp_dn = .true.

      if ( all(U < ZERO+1.d-6) )then
         !Negative U and U=0
         Hub%neg_u = .true.

         ! G_up and G_dn are identical. Do not compute G_dn.
         if( maxval(abs(t_up-t_dn)) < 1.d-6 .and. maxval(abs(mu_up-mu_dn)) < 1.d-6 ) then
            Hub%comp_dn = .false.
         endif
      elseif ( all(U > ZERO-1.d-6))then
         !Positive U 
         Hub%neg_u   = .false.

         ! Don't compute G_dn on half-filled (mu=0) bipartite lattices
         if ( all(abs(mu_up) < 1.d-6) .and. all(abs(mu_dn) < 1.d-6) .and. &
              maxval(abs(t_up-t_dn)) < 1.d-6 .and. Hub%S%checklist(STRUCT_PHASE) ) then
            Hub%comp_dn = .false.
         end if
      else
         stop 'All U''s must have the same sign (or be zero)'
      end if

    else if (SimType .eq. Holstein_model) then
      if ( any(U < -1.d-6) )then
        write(*,*)
        write(*,*) ' Electron-phonon couplings must be positive.'
        write(*,*) ' Stop.'
        stop
      end if

      ! No need to compute G_dn for in the Holstein mode.
      Hub%comp_dn = .false.
      Hub%neg_u   = .true.

      ! Normalization factors of the phonon action
      norm_phonon(1) = 0.5_wp/dtau
      norm_phonon(2) = 0.5_wp*omega*omega*dtau

      ! Length scale of phonon displacement
      CHSF0 = 1.0_wp /omega/(exp(float(L)*dtau*omega)-1.0_wp) + 0.5_wp/omega
      CHSF0 = 2.0_wp*sqrt(CHSF0)
    end if

    !write(*,*)
    !write(*,*) "DEBUG info -- In DQMC_Hub_Init():"
    !write(*,*) "  Hub%comp_dn=",Hub%comp_dn
    !write(*,*) "  Hub%neg_u=",Hub%neg_u
    !write(*,*) "  S%P defined? ",Hub%S%checklist(STRUCT_PHASE)
    !write(*,*)

    ! Parameters for MC loop
    Hub%nWarm    = nWarm
    Hub%nPass    = nPass
    Hub%nMeas    = nMeas
    Hub%nTry     = nTry
    Hub%tausk    = ntausk

    ! Initialize random seeds
    Hub%idum     = seed
    if (Hub%idum == 0) then
       call date_and_time(VALUES=val)
       Hub%idum = getpid()+val(8)*val(7)+val(6)**mod(val(5),5)
    end if

    ! LAPACK random variable generation
    Hub%seed = Hub%idum * (/1,2,3,4/)
    Hub%seed = mod(abs(Hub%seed), 4095)
    if (mod(Hub%seed(4),2) == 0) then
       Hub%seed(4) = Hub%seed(4) + 1
    end if
#   ifdef _QMC_MPI
       junkPtr = init_sprng(SPRNG_LCG, Hub%seed(4), SPRNG_DEFAULT)
#   endif

    ! Initialize auxiliary variables
    Hub%gamma   = gamma
    Hub%nAccept = accept
    Hub%nReject = reject

    ! Initialize working space 
    call DQMC_WSpace_Allocate(n, Hub%S%n_b, Hub%WS)

    ! Initialize Hubbard-Stratonovich Field
    HSF = HSF_IPT
    Hub%HSFtype = HSFtype
    if (HSF == HSF_FROM_MEMORY) then
       ilb = Hub%G_up%ilb
       ! discrete case
       if (HSFtype == HSF_DISC) then
          if (.not. associated(Hub%HSF)) then
             call DQMC_Warning("Cannot use current HSF. ", 0)
             HSF = HSF_RANDOM_GEN
          else
             print *, "Read HSF from memory."
          end if
       else
          ! contnuous case
          if (.not. associated(Hub%CHSF)) then
             call DQMC_Warning("Cannot use current HSF. ", 0)
             HSF = HSF_RANDOM_GEN
          else
             print *, "Read HSF from memory."
          end if
       end if
    end if

    if (HSF == HSF_FROM_FILE .or. HSF == HSF_RESTORE) then
       inquire(UNIT=HSF_INPUT_UNIT, EXIST=lex)
       if (lex) then
          ! If a valid input file handle is provided,
          ! read HSF from the file
          if (HSFtype == HSF_DISC) then
             allocate(Hub%HSF(n,L))
             call DQMC_Hub_Input_HSF(Hub, HSF==HSF_RESTORE, ilb, HSF_INPUT_UNIT)
          else
             ! TODO: input continuous HSF from file
             call DQMC_Error("reading continuous HSF from file is not supported", HSF)
          end if
          print *, "Read HSF from a file."
       else
          ! If file does not exist, give a warning message.
          call DQMC_Warning("HSF file does not exist. &
               & Use randomly generated values.", HSF)
          HSF = HSF_RANDOM_GEN
       end if
    end if

    ! generate HSF randomly
    if (HSF == HSF_RANDOM_GEN) then
       ilb = 1

       ! discrete case
       if (HSFtype .eq. HSF_DISC) then
          allocate(Hub%HSF(n,L))
          Hub%HSF = 1
          do i = 1, Hub%L   
             call ran0(n, Hub%WS%R5, Hub%seed)
             where(Hub%WS%R5 > HALF) Hub%HSF(:,i) = -1
          end do
       else if (HSFtype .eq. HSF_CONT) then
          ! continuous case
          allocate(Hub%CHSF(n,L))
          if ( SimType .eq. Hubbard_model ) then
            do i = 1, Hub%L   
               call ran1(n, Hub%CHSF(:,i), Hub%seed)
            end do          
          else if ( SimType .eq. Holstein_model ) then
            call ran1(n, CHSFtmp, Hub%seed)
            do i = 1, Hub%L
              Hub%CHSF(:,i) = CHSF0*CHSFtmp
            end do
          end if
       end if
    end if

    ! Initialize lookup table
    if  (HSFtype .eq. HSF_DISC) then
       nullify(Hub%explook)
       allocate(Hub%explook(-2:2,1:n_U))
       do j = 1, n_U
          temp = exp(dtau*abs(U(j))*HALF)    
          lambda = log(temp+sqrt(temp*temp-ONE))
          do i = -2, 2
             Hub%explook(i,j)=exp(float(i)*lambda)
          end do
          ! for use by U<0, save lambda in explook(0, j)
          Hub%explook(0, j) = lambda
       end do
    else if (HSFtype .eq. HSF_CONT) then
       allocate(Hub%lambda(n_U))
       ! Hubbard model mode.
       if (SimType .eq. Hubbard_model) then
         do j = 1, n_U
            Hub%lambda(j) = sqrt(dtau*abs(U(j)))
         end do
       ! Holstein model mode. Here U should be understood as the electron-phonon couplings.
       elseif (simType .eq. Holstein_model) then
         do j = 1, n_U
            Hub%lambda(j) = U(j)*dtau
         end do
       end if
    end if

    ! Initialize V matrices
    call DQMC_Hub_Init_Vmat(Hub)

    ! Initialize Green functions
    if ((SimType .eq. Hubbard_model) .and. Hub%neg_u .and. (HSFType .eq. HSF_CONT)) then
      call DQMC_B_Init(n, Hub%B_up, Hub%WS, Hub%S%T, Hub%S%ckb, Hub%t_up, Map_up, dtau, Uhubb)
      call DQMC_B_Init(n, Hub%B_dn, Hub%WS, Hub%S%T, Hub%S%ckb, Hub%t_dn, Map_dn, dtau, Uhubb)
    else
      call DQMC_B_Init(n, Hub%B_up, Hub%WS, Hub%S%T, Hub%S%ckb, Hub%t_up, Map_up, dtau)
      call DQMC_B_Init(n, Hub%B_dn, Hub%WS, Hub%S%T, Hub%S%ckb, Hub%t_dn, Map_dn, dtau)
    end if

    call DQMC_SeqB_Init(n, Hub%L, nOrth, Hub%B_up, Hub%SB_up, Hub%WS)
    call DQMC_SeqB_Init(n, Hub%L, nOrth, Hub%B_dn, Hub%SB_dn, Hub%WS)

    ! Initialize G
    call DQMC_GFun_Init(n, L, Hub%G_up, Hub%V_up,  Hub%WS, &
       nWrap, difflim, errrate, GMAT_UP, ssxx, fixw)

    ! for positive U or H_dn/=H_up we need to construct G_dn implicitly
    if (Hub%comp_dn .or. .not. Hub%neg_u) then
       call DQMC_GFun_Init(n, L, Hub%G_dn, Hub%V_dn,  Hub%WS, &
          nWrap, difflim, errrate, GMAT_DN, ssxx, fixw)
    else
       ! Negative U or U=0, G_dn is a clone of G_up
       call DQMC_Gfun_Clone(Hub%G_dn, Hub%G_up)
    end if

    !Fill G
    call DQMC_GetG(ilb, Hub%G_up, Hub%SB_up)
    if (Hub%comp_dn .or. .not.Hub%neg_u) then
       call DQMC_GetG(ilb, Hub%G_dn, Hub%SB_dn)
    end if

    ! Initialize measurements
    temp = Hub%dtau * Hub%L
    call DQMC_Phy0_Init(Hub%P0, Hub%S, temp, nBin, Hub%WS)
    call DQMC_Phy2_Init(Hub%P2, nBin, Hub%S, Hub%WS, Hub%meas2)

    ! Initialize simulation range
    Hub%n_start = 1
    Hub%n_end   = n

    deallocate(CHSFtmp)

  end subroutine DQMC_Hub_Init

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Free(Hub)
    !
    ! Purpose
    ! =======
    !    This subrotine deallocate variables in Hubbard
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout) :: Hub  ! Hubbard model

    ! ... Executable ...

    call DQMC_B_Free(Hub%B_up)
    call DQMC_B_Free(Hub%B_dn)
    call DQMC_Gfun_Free(Hub%G_up)
    call DQMC_Gfun_Free(Hub%G_dn)

    call DQMC_Phy0_Free(Hub%P0)
    call DQMC_Phy2_Free(Hub%P2)

    if (associated(Hub%V_up)) deallocate(Hub%V_up)
    if (.not.Hub%neg_u .and. associated(Hub%V_dn)) deallocate(Hub%V_dn)
    deallocate(Hub%t_up, Hub%t_dn, Hub%mu_up, Hub%mu_dn, Hub%U)

    if (Hub%HSFtype == HSF_DISC) then
       deallocate(Hub%HSF)
       deallocate(Hub%explook)
    else
       deallocate(Hub%CHSF)
       deallocate(Hub%lambda)
    end if

    call DQMC_WSpace_Free(Hub%WS)
    call DQMC_SeqB_Free(Hub%SB_up)
    call DQMC_SeqB_Free(Hub%SB_dn)
    call DQMC_Struct_Free(Hub%S)

  end subroutine DQMC_Hub_Free

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Output_HSF(Hub, restore, slice, OPT)
    !
    ! Purpose
    ! =======
    !    This subrotine outputs Hubbard-Stratonovich Field to a
    !    output file OPT.
    !
    ! Arguments
    ! =========
    !
    use dqmc_mpi
#   ifdef _QMC_MPI
#      include "sprng_f.h"
#   endif
    type(Hubbard), intent(in) :: Hub
!     integer, intent(in) :: n, L         ! dim of HSF
!     integer, intent(in) :: HSF(n,L)     ! Hubbard-Stratonovich Field
!     integer, intent(in) :: seed(4)      ! RNG state
    integer, intent(in) :: OPT          ! output handle
    logical, intent(in) :: restore
    integer, intent(in) :: slice

    ! ... local variables ...
    integer :: i, j, k, h, nproc, nL
    integer :: HSF(Hub%n*Hub%L+1)
#   ifdef _QMC_MPI
       character :: rndbuf(MAX_PACKED_LENGTH) 
       character(len=60) :: fname
       integer   :: stat
#   endif

    ! ... Executable ....
    nL = Hub%n * Hub%L
    nproc = qmc_sim%aggr_size
   

    if (qmc_sim%rank == qmc_sim%aggr_root) then
       !loop over processors
       do j = 0, nproc-1
          ! Collect data from processor "j"
          if (j /= qmc_sim%aggr_root) then
#            ifdef _QMC_MPI
                call mpi_recv(HSF, nL+1, MPI_INT, j, j, MPI_COMM_WORLD, stat, k)
#            endif
          else 
             call pack_fields
          endif
          !Write data to disk
          do k = 1, nL
             write(OPT,'(i1)',advance='no') HSF(k)
          enddo
          write(OPT,*) HSF(nL+1)
       enddo
    else
      !Send data to root node
      call pack_fields
#     ifdef _QMC_MPI
         ! Send field to root processor
         call mpi_send(HSF, nL+1, MPI_INT, qmc_sim%aggr_root, qmc_sim%rank, MPI_COMM_WORLD, k)
#     endif
    endif
    
    !If restore is true, write the random number generator status
    if (restore) then
#      ifdef _QMC_MPI
          k = pack_sprng(rndbuf)
          inquire(unit=OPT, name=fname)
          close(OPT)
          !This would be better done with send/recv (see commented code below)
          ! but there appears to be a problem with the written rndbuf (!?)
          do i = 0, nproc - 1
             if (i == qmc_sim%aggr_rank) then
                open(file=fname, form='formatted', unit=OPT, position='append')
                write(OPT,'(i8,1x)',advance='no') k
                do j = 1, k
                   write(OPT,'(A1)',advance='no') rndbuf(j)
                enddo
                write(OPT,*)
                close(OPT)
             endif
             call mpi_barrier(MPI_COMM_WORLD,j)
          enddo
          open(file=fname, form='formatted', unit=OPT, position='append')
#      else
          !Simply write the seed if MPI is off
          write(OPT, *) Hub%seed
#      endif
       write(OPT, *) Hub%gamma, Hub%naccept, Hub%nreject
    endif


!    if (DQMC_MPI_Is_Root(qmc_sim, CHANNEL_AGGR) .and. restore) then
!#      ifdef _QMC_MPI
!          do j = 0, nproc - 1
!             ! Collect data from processor "j"
!             if (j /= qmc_sim%aggr_root) then
!                call mpi_recv(k, 1, MPI_INT, j, j, MPI_COMM_WORLD, stat, i)
!                call mpi_recv(rndbuf, k, MPI_CHARACTER, j, j, MPI_COMM_WORLD, stat, i)
!             else
!                k = pack_sprng(rndbuf)
!             endif
!             !Write data to disk
!             write(OPT,'(i8,1x)',advance='no') k
!             do i = 1, k
!                write(OPT,'(A1)',advance='no') rndbuf(i)
!             enddo
!             write(OPT,*)
!          enddo
!#      else
!          !Simply write the seed if MPI is off
!          write(OPT, *) Hub%seed
!#      endif
!       !Write info to reastablish acceptance rate
!       write(OPT, *) Hub%gamma, Hub%naccept, Hub%nreject
!    elseif (restore) then
!#      ifdef _QMC_MPI
!          k = pack_sprng(rndbuf)
!          call mpi_send(k, 1, MPI_INT, qmc_sim%aggr_root, qmc_sim%rank, MPI_COMM_WORLD, k)
!          call mpi_send(rndbuf, k, MPI_CHARACTER, qmc_sim%aggr_root, qmc_sim%rank, MPI_COMM_WORLD, k)
!#      endif
!    endif

#   ifdef _QMC_MPI
       call mpi_barrier(MPI_COMM_WORLD, i)
#   endif
 
  contains 

     subroutine pack_fields
        !Fill the vector of fields to write
        k = 1
        do i = 1, Hub%L
           do h = 1, Hub%n
              HSF(k) = (Hub%HSF(h,i)+1)/2
              k = k + 1
           enddo
        enddo
        !Last element is the time slice
        if (slice <= 0) then
           HSF(k) = Hub%G_up%ilb
        else
           HSF(k) = slice
        endif
     end subroutine pack_fields

  end subroutine DQMC_Hub_Output_HSF

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Input_HSF(Hub, restore, slice, IPT)
    !
    ! Purpose
    ! =======
    !    This subrotine reads Hubbard-Stratonovich Field from a
    !    file OPT.
    !
    ! Arguments
    ! =========
    !
    use dqmc_mpi
#   ifdef _QMC_MPI
#      include "sprng_f.h"
#   endif
    type(Hubbard), intent(inout) :: Hub
    logical, intent(in)  :: restore      ! restore RNG?
    integer, intent(in)  :: IPT          ! input handle
    integer, intent(out)  :: slice
    ! ... local varaible ...
    Integer :: i, j, k, ip, nproc, nL
    integer :: HSF(Hub%n*Hub%L+1)
#   ifdef _QMC_MPI
       character :: rndbuf(MAX_PACKED_LENGTH)
       SPRNG_POINTER :: junkPtr
#   endif

    ! ... Executable ....
    nL = Hub%n * Hub%L
    nproc = qmc_sim%size

    do ip = 0, nproc-1
       if (qmc_sim%aggr_rank == ip) then
          k = 1
          do i = 1, Hub%L
             do j = 1, Hub%n
                read(IPT,'(i1)',advance='no',ERR=100) HSF(k)
                k = k + 1
             enddo
             write(98,*)
          enddo
          read(IPT,*) HSF(k)
       else
          read(IPT,*)
       endif
       k = 1
       do i = 1, Hub%L
          do j = 1, Hub%n
             Hub%HSF(j,i) = 2*HSF(k) - 1
             k = k + 1
          enddo
       enddo
       slice = HSF(k)
    enddo
    
    if (restore) then
       do ip = 0, nproc-1
          if (qmc_sim%aggr_rank == ip) then
#            ifdef _QMC_MPI
                read(IPT,'(i8,1x)',advance='no') k
                do j = 1, k
                   read(IPT,'(A1)',advance='no') rndbuf(j)
                enddo
                read(IPT,*)
#            else 
                read(IPT, *) Hub%seed
#            endif
          else
             read(IPT,*)
          endif
       enddo
       read(IPT, *,ERR=102) Hub%gamma, Hub%naccept, Hub%nreject
#      ifdef _QMC_MPI
          junkPtr = unpack_sprng(rndbuf)
#      endif

    endif

    return

100 call DQMC_Error("cannot read HSF input file:", HSF_INPUT_UNIT)
102 call DQMC_Error("cannot read gamma/naccept/nreject from HSF input file", HSF_INPUT_UNIT)

  end subroutine DQMC_Hub_Input_HSF

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_OutputParam(Hub, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subrotine outputs parameters of Hubbard model and
    !    computed results.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(in) :: Hub     ! Hubbard model
    integer, intent(in)       :: OPT     ! output handle
    
    ! ... Local ...
    character(35) :: FMT
    logical       :: lex
    integer, parameter :: slice = 0
    logical, parameter :: restore = .true.

    ! ... Executable ....

    if (qmc_sim%rank == qmc_sim%aggr_root) then

       write(OPT,*)  Hub%S%Name(:)
       if (Hub%n_U == 1) then
          FMT = FMT_STRDBL
       else
          write(FMT, "('(a30,f19.6,(',I3,'(f12.6)))')") Hub%n_U - 1
       end if
       if (SimType .eq. Hubbard_model) then
         write(OPT,FMT)         "                          U : ", Hub%U
       elseif (SimType .eq. Holstein_model) then
         write(OPT,FMT)         "      el-ph coupling lambda : ", Hub%U
         write(OPT,FMT)         "     phonon frequency omega : ", Hub%omega
       end if
       if (Hub%n_t == 1) then
          FMT = FMT_STRDBL
       else
          write(FMT, "('(a30,f19.6,(',I3,'(f12.6)))')") Hub%n_t - 1
       end if
       write(OPT,FMT)         "                       t_up : ", Hub%t_up
       write(OPT,FMT)         "                       t_dn : ", Hub%t_dn
       if (Hub%n_mu == 1) then
          FMT = FMT_STRDBL
       else
          write(FMT, "('(a30,f19.6,(',I3,'(f12.6)))')") Hub%n_mu - 1
       end if
       write(OPT,FMT)         "                      mu_up : ", Hub%mu_up
       write(OPT,FMT)         "                      mu_dn : ", Hub%mu_dn
       write(OPT,FMT_STRINT)  "             Time slice - L : ", Hub%L
       write(OPT,FMT_STRINT)  "            Number of sites : ", Hub%n
       write(OPT,FMT_STRDBL)  "                       dtau : ", Hub%dtau
       write(OPT,FMT_STRDBL)  "                       beta : ", Hub%dtau * Hub%L
       write(OPT,FMT_STRINT)  "     Number of warmup sweep : ", Hub%nWarm
       write(OPT,FMT_STRINT)  "Number of measurement sweep : ", Hub%nPass
       write(OPT,FMT_STRINT)  "   Frequency of measurement : ", Hub%nMeas
       write(OPT,FMT_STRINT)  "                Random seed : ", Hub%idum
       write(OPT,FMT_STRINT)  " Frequency of recomputing G : ", Hub%G_up%nWrap
       if (Hub%nTry .gt. 0) then
         write(OPT,FMT_STRINT)  "Global move number of sites : ", Hub%nTry
       end if
       write(OPT,FMT_STRINT)  "               Accept count : ", Hub%naccept
       write(OPT,FMT_STRINT)  "               Reject count : ", Hub%nreject
       write(OPT,FMT_STRDBL)  "    Approximate accept rate : ", &
            dble(Hub%naccept)/dble(Hub%naccept+Hub%nreject)
       write(OPT,FMT_STRDBL)  "                      gamma : ", Hub%gamma
       if (Hub%nTry > 0) then         
         write(OPT,FMT_STRINT)  "   Global move accept count : ", Hub%nAcceptGlobal
         write(OPT,FMT_STRINT)  "   Global move reject count : ", Hub%nRejectGlobal
         write(OPT,FMT_STRDBL)  "    Global move accept rate : ", &
              dble(Hub%nAcceptGlobal)/dble(Hub%nAcceptGlobal + Hub%nRejectGlobal)
       end if
       if (SimType .eq. Holstein_model) then
         write(OPT,FMT_STRINT)  " Global move 2 accept count : ", Hub%nAcceptGlobal2
         write(OPT,FMT_STRINT)  " Global move 2 reject count : ", Hub%nRejectGlobal2
         write(OPT,FMT_STRDBL)  "  Global move 2 accept rate : ", &
              dble(Hub%nAcceptGlobal2)/dble(Hub%nAcceptGlobal2 + Hub%nRejectGlobal2)
       end if
       write(OPT,*)           "          Type of matrix B : ", Hub%B_up%name
       if (Hub%HSFtype == HSF_DISC) then
          write(OPT,*)        "        Type of matrix HSF : discrete"
       else
          write(OPT,*)        "        Type of matrix HSF : continuous"
          write(OPT,*)        "                   delta 1 : ", Hub%delta1
          write(OPT,*)        "                   delta 2 : ", Hub%delta2
       end if
       if ( SimType .eq. Holstein_model ) then
          write(OPT,*)        "        phonon field scale : ", CHSF0
       end if

    endif

    ! Check if the file is valid.
    if (Hub%outputHSF) then
       inquire(UNIT=HSF_OUTPUT_UNIT, EXIST=lex)
       if (lex) then
          call DQMC_Hub_Output_HSF(Hub, restore, slice, HSF_OUTPUT_UNIT)
       else
          if (qmc_sim%rank == qmc_sim%aggr_root) &
             call DQMC_Warning("HSF output file is not initialized.", 1)
       end if
    end if

  end subroutine DQMC_Hub_OutputParam
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Print(Hub, OPT)
    implicit none
    !
    ! Purpose
    ! =======
    !    This subrotine outputs parameters of Hubbard model and
    !    computed results.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(in) :: Hub     ! Hubbard model
    integer, intent(in)       :: OPT     ! output handle
    
    ! ... Executable ....

    call DQMC_Hub_OutputParam(Hub, OPT)
    write(OPT, FMT_DBLINE)

    call DQMC_Phy0_Print(Hub%P0, Hub%S, SimType, OPT)
    call DQMC_Phy2_Print(Hub%P2, Hub%S%wlabel, OPT)

  end subroutine DQMC_Hub_Print

  ! --------------------------------------------------------------------!

  subroutine DQMC_Hub_Sweep(Hub, nMeas0)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !   This subroutine performs the DQMC sweep, which is consisted of 
    !   four steps. (See [1] for more details.)
    !
    !      1. Swap the slice of G and recompute G if necessary.
    !      2. Metropolis Algorithm
    !      3. Update the model and perform physical measurement.
    !      4. Adjust parameters.
    !
    !   The first three steps are within a big loop, which run
    !   through each time slice of G. The major part is the second
    !   step, which is explained below.
    !
    !      1. Try the new configuration by single spin-flip sampling 
    !         at site j at time slice i.
    !      2. Compute the probability of this new configuration.
    !         
    !             p =  r/(1+gamma*r)    if r < 1
    !             p =  r/(gamma+r)      if r >=1
    !        
    !         where r is the ratio of determinants of Green's function
    !         of spin up and spin down.
    !      3. If p > ran, a uniform random number in [0,1], then change
    !         the configuration and update the Green's function.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout),target  :: Hub      ! Hubbard class
    integer, intent(in)                  :: nMeas0   ! Duration of measurement

    ! ... paremeters ...
    integer, parameter  :: DQMC_CHECK_ITER  = 10000
    integer, parameter  :: DQMC_ADJUST      = 100
    real(wp), parameter :: DQMC_ACC_UP      = 0.52_wp
    real(wp), parameter :: DQMC_ACC_LO      = 0.48_wp

    ! ... local scalar ...

    integer  :: i, j, n, L, m      ! Loop iterator
    integer  :: accept_cnt            ! Counter for accept in Met-alg
    integer  :: reject_cnt            ! Counter for accept in Met-alg
    real(wp) :: accrat
    real(wp) :: alpha_up, alpha_dn=0  ! Change of configuration
    real(wp) :: p, randn              ! Probability of changing
    real(wp) :: r_up, r_dn, r         ! Ratio of determinant
    real(wp) :: gjj                   ! (j,j) element of G_up or G_dn
#   ifdef _QMC_MPI
       integer  :: send_cnt(2) 
       integer  :: recv_cnt(2)           
#   endif
    
    ! To speed up the computation
    real(wp), pointer :: G_up(:,:) 
    real(wp), pointer :: G_dn(:,:) 
    real(wp), pointer :: U_up(:,:) 
    real(wp), pointer :: U_dn(:,:) 
    real(wp), pointer :: W_up(:,:) 
    real(wp), pointer :: W_dn(:,:) 
    real(wp), pointer :: V_up(:,:) 
    real(wp), pointer :: V_dn(:,:) 
    integer,  pointer :: blksz_up  
    integer,  pointer :: blksz_dn  
    real(wp), pointer :: sgn_up    
    real(wp), pointer :: sgn_dn    

    real(wp), pointer :: ranlist(:)   
    real(wp), pointer :: explook(:,:) 
    integer,  pointer :: HSF(:,:)     
    integer,  pointer :: map(:)       
    real(wp)  :: gamma
    logical   :: comp_dn, neg_u

    ! ... Executable ...

    !=====================! 
    ! Step 0: Setup alias !
    !=====================!

    G_up     => Hub%G_up%G
    U_up     => Hub%G_up%U
    W_up     => Hub%G_up%W
    V_up     => Hub%G_up%V
    blksz_up => Hub%G_up%blksz
    sgn_up   => Hub%G_up%sgn

    G_dn     => Hub%G_dn%G
    U_dn     => Hub%G_dn%U
    W_dn     => Hub%G_dn%W
    V_dn     => Hub%G_dn%V
    blksz_dn => Hub%G_dn%blksz
    sgn_dn   => Hub%G_dn%sgn

    ranlist  => Hub%WS%R7
    explook  => Hub%explook
    HSF      => Hub%HSF
    map      => Hub%S%map
    gamma    =  Hub%gamma
    comp_dn  =  Hub%comp_dn
    neg_u    =  Hub%neg_u

    n = Hub%n    
    L = Hub%L

    i = Hub%G_up%ilb
    accept_cnt = 0
    reject_cnt = 0

#   ifdef _QMC_MPI
       send_cnt(1) = Hub%naccept
       send_cnt(2) = Hub%nreject
#   endif

    do m = 1, L

       !First thing see if you can make a measurement on a freshly computed G
       if (Hub%G_up%wps==Hub%G_up%nWrap) then

          if(nmeas0>0)then
             ! Construct G_dn for mu = 0 and U > 0 using particle-hole symmetry
             ! Note that ( .not.neg_u .and. not.comp_dn ) implies that S%P is defined, i.e.
             ! S%checklist(STRUCT_PHASE) = 'T'.
             if (.not.neg_u .and. .not.comp_dn) call DQMC_GFun_CopyUp(Hub%G_dn, Hub%G_up, Hub%S%P)

             !Fill GS in Hub%Gfun
             call DQMC_GetG_2nd_order(Hub%G_up, Hub%B_up)
             if (comp_dn .or. .not.neg_u) then
                call DQMC_GetG_2nd_order(Hub%G_dn, Hub%B_dn)
             endif

             ! Basic measurement
             call DQMC_Phy0_Meas(Hub%n, Hub%P0, Hub%G_up%GS, Hub%G_dn%GS, Hub%U,   &
                  Hub%mu_up, Hub%mu_dn, Hub%t_up, Hub%t_dn, sgn_up, sgn_dn, Hub%S)

             if (Hub%meas2) then
                ! Pair measurement
                r = sgn_up*sgn_dn
                call DQMC_Phy2_Meas(n, Hub%P2%M1, Hub%P2%M2, &
                     Hub%P2, Hub%S%B, Hub%G_up%GS, Hub%G_dn%GS, r)
             end if

          endif

       end if

       i = i + 1
       if( i>L .or. i<0 ) i=1

       !==============================! 
       ! Step 1: Swap the slice of G  !
       !==============================!

       call DQMC_GetG(i, Hub%G_up, Hub%SB_up)
       if (comp_dn) then
          call DQMC_GetG(i, Hub%G_dn, Hub%SB_dn)
       else
          sgn_dn = sgn_up
       end if

       !==============================!
       ! Step 2: Metropolis Algorithm !
       !==============================!
       call ran0(n, ranlist, Hub%seed)

       ! Loop over lattice sites
       do j = Hub%n_start, Hub%n_end
          ! Try the new configuration by single spin-flip sampling 
          ! at site j at time slice i.
          ! See reference [1] for more detail for these formula
          if (neg_u) then
             alpha_up = explook(-2*HSF(j,i), map(j)) - ONE
             alpha_dn = alpha_up
          else
             alpha_up = explook(-2*HSF(j,i), map(j)) - ONE
             alpha_dn = explook( 2*HSF(j,i), map(j)) - ONE
          end if

          gjj = DQMC_Gfun_Getjj(n, j, blksz_up, G_up, U_up, W_up)

          r_up = ONE + (ONE - gjj)*alpha_up

          if (comp_dn) then
             gjj = DQMC_Gfun_Getjj(n, j, blksz_dn, G_dn, U_dn, W_dn)
             r_dn = ONE + (ONE - gjj)*alpha_dn
          elseif (neg_u) then
             r_dn = r_up
          else
             r_dn = ONE + gjj*alpha_dn
          end if

          r = abs(r_up * r_dn)

          if (neg_u) then
             r = r * explook(+2*HSF(j, i), map(j))
          end if

          ! Compute the probability
          if(r <= ONE) then
             p = r/(ONE+gamma*r)
          else
             p = r/(gamma+r)
          end if

          randn = ranlist(j)

          ! Accept 
          if (p > randn) then
             accept_cnt = accept_cnt + 1

             if(r_up < ZERO) sgn_up = -sgn_up
             if(r_dn < ZERO) sgn_dn = -sgn_dn
             HSF(j,i) = -HSF(j,i)

             ! Update G_up
             call DQMC_UpdateG(j, alpha_up/r_up, Hub%G_up)
             ! invalidate the cache for B_i
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, i)
#endif
             V_up(j,i) = V_up(j,i) * (alpha_up + ONE)
             Hub%G_up%nModify = i
             Hub%G_up%det = Hub%G_up%det - log(abs(r_up))

             ! If mu /= zero, then update G_dn as well.
             if (comp_dn) then
                ! Update G_dn
                call DQMC_UpdateG(j,  alpha_dn/r_dn, Hub%G_dn)
                Hub%G_dn%det = Hub%G_dn%det - log(abs(r_dn))
             end if
             if (.not. neg_u) then
#if defined(DQMC_ASQRD)
                call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, i)
#endif
                V_dn(j,i) = V_dn(j,i) * (alpha_dn + ONE)
             end if
             Hub%G_dn%nModify = i

          else 
             reject_cnt = reject_cnt + 1

          endif

       end do


       !============================!
       ! Step 3: Update and Measure !
       !============================!
       ! update G_up/G_dn if there are some updates not applied.

       call DQMC_ApplyUpdate(Hub%G_up, forced = .true.)
       if (comp_dn) then
          call DQMC_ApplyUpdate(Hub%G_dn, forced = .true.)
       end if

    end do

#   ifdef _QMC_MPI
       send_cnt(1) = accept_cnt
       send_cnt(2) = reject_cnt
       call mpi_allreduce(send_cnt, recv_cnt, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD, m)
       accept_cnt = recv_cnt(1)
       reject_cnt = recv_cnt(2)
#   endif

    ! update accept and reject counts
    Hub%naccept = Hub%naccept + accept_cnt
    Hub%nreject = Hub%nreject + reject_cnt

    !===========================!
    ! Step 4: Adjust parameters !
    !===========================!
    if(Hub%naccept+Hub%nreject > DQMC_CHECK_ITER) then
       accrat = dble(Hub%naccept)/dble(Hub%naccept+Hub%nreject)
       if(accrat > DQMC_ACC_UP .or. accrat < DQMC_ACC_LO)then
          Hub%gamma = Hub%gamma + (accrat - HALF)
          Hub%gamma = dmax1(ZERO,Hub%gamma)
          Hub%gamma = dmin1(ONE, Hub%gamma)
          Hub%naccept = int(DQMC_ADJUST*accrat)
          Hub%nreject = int(DQMC_ADJUST*(ONE-accrat))
       endif
    endif

    !Update nwrap. Stop doing it if measurements are performed 
    !every nwrap to ensure bins contains same number of measurements.
    if(nmeas0 <= 0) then
      call DQMC_UpdateWraps(Hub%G_up)
      if (comp_dn) then
         call DQMC_UpdateWraps(Hub%G_dn)
         call DQMC_SyncWraps(Hub%G_up, Hub%G_dn)
      end if
    endif
    

  end subroutine DQMC_Hub_Sweep

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Sweep2(Hub, numTry)
    !
    ! Purpose
    ! =======
    !   This subroutine performs the global moves of DQMC sweep, in which 
    !   all the Hub(i) on some selected sites are flipped for all slice.
    !
    !      1. Try the new configuration.
    !      2. Compute the probability of this new configuration.
    !         
    !             p =  r/(1+gamma*r)    if r < 1
    !             p =  r/(gamma+r)      if r >=1
    !        
    !         where r is the ratio of determinants of Green's function
    !         of spin up and spin down.
    !      3. If p > ran, a uniform random number in [0,1], then change
    !         the configuration and update the Green's function.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout)  :: Hub      ! Hubbard model
    integer, intent(in)           :: numTry   ! Number of Try

    ! ... Local Variables ...
    real(wp) :: ranList(numTry), rat, ratexp, ranSlice(1)
    integer  :: i, j, n, L, accept, tmp, nSite, tslice
    real(wp) :: det_up, det_dn, new_up, new_dn
    real(wp) :: copy_sgn_up, copy_sgn_dn
    integer, pointer :: map(:) 
    integer  :: siteList(Hub%n), site(numTry)
    integer  :: hs_sum          ! sum over HS fields at site, for U<0
    
    ! ... Executable ...

    n = Hub%n
    L = Hub%L
    accept = 0
    if (numTry <= 0) return
    Map=> Hub%S%Map

    call ran0(1, ranSlice, Hub%seed)
    tslice=ceiling(ranSlice(1)*L)

    ! Compute the Green's matrix and the sign
    Hub%G_up%ilb  = -1
    call DQMC_GetG(tslice, Hub%G_up, Hub%SB_up)
    det_up = Hub%G_up%det
    if ( Hub%comp_dn ) then
       Hub%G_dn%ilb  = -1
       call DQMC_GetG(tslice, Hub%G_dn, Hub%SB_dn)
       det_dn = Hub%G_dn%det
    elseif ( Hub%neg_u ) then
       det_dn = det_up
       Hub%G_dn%sgn = Hub%G_up%sgn
    else
       ! Note that here we have .not.neg_u and not.comp_dn.
       ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
       ! So we can safely call DQMC_Gfun_CopyUp() which uses particle-hole symmetry.
       call DQMC_Gfun_CopyUp(Hub%G_dn, Hub%G_up, Hub%S%P)
       det_dn = Hub%G_dn%det
    end if

    ! get random numbers
    call ran0(numTry, ranList, Hub%seed)

    nsite = Hub%n_End - Hub%n_start + 1
    siteList(1:nSite) = Hub%n_Start+(/(i,i=0,nSite-1)/)

    ! generate sites
    do i = 1, numtry
       tmp = int(ranList(i)*nSite) + 1
       site(i) = siteList(tmp)
       ! compress the list
       do j = tmp+1, nSite
          siteList(j-1) = siteList(j) 
       end do
       nSite = nSite - 1
    end do

    call ran0(numTry, ranList, Hub%seed)

    ! Global move
    do i = 1, numTry
       ! Flip its HS field for all the slices
       do j = 1, L
          tmp = -Hub%HSF(site(i),j)
          Hub%HSF (site(i),j) = tmp
#if defined(DQMC_ASQRD)
          call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, j)
#endif
          Hub%V_up(site(i),j) = Hub%explook( tmp, map(site(i)))

          if (.not. Hub%neg_u) then
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, j)
#endif
             Hub%V_dn(site(i),j) = Hub%explook(-tmp, map(site(i)))
          end if
       end do
       
       ! Store the value of G first
       Hub%G_up%tmp = Hub%G_up%G
       if (Hub%comp_dn) then
          Hub%G_dn%tmp = Hub%G_dn%G
       end if
       copy_sgn_up = Hub%G_up%sgn
       copy_sgn_dn = Hub%G_dn%sgn

       ! Compute G with new configuration
       Hub%G_up%ilb  = -1
       call DQMC_GetG(tslice, Hub%G_up, Hub%SB_up)
       new_up = Hub%G_up%det
       if ( Hub%comp_dn ) then
          Hub%G_dn%ilb  = -1
          call DQMC_GetG(tslice, Hub%G_dn, Hub%SB_dn)
          new_dn = Hub%G_dn%det
       elseif ( Hub%neg_u ) then
          Hub%G_dn%sgn = Hub%G_up%sgn
          new_dn = new_up
       else
          ! Note that here we have (.not.neg_u) and (.not.comp_dn).
          ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
          ! So we can safely call DQMC_Gfun_CopyUp() which uses particle-hole symmetry.
          call DQMC_Gfun_CopyUp(Hub%G_dn, Hub%G_up, Hub%S%P)
          new_dn = Hub%G_dn%det
       end if

       ! Compute the Det ratio
       ! NB: the determinant computed by GetG is log(abs(det(G)))
       !     Here we need log(abs(Z))= -log(abs(det(G)))
       rat = det_up + det_dn - new_up - new_dn

       if (Hub%neg_u) then
          hs_sum=0
          do j = 1, L
             hs_sum = hs_sum + Hub%HSF(site(i), j)
          end do

          ! extra factor for U<0: exp(2 * lambda * Sum(l, s_new(i, l)))
          rat = rat - 2*Hub%explook(0, map(site(i))) * hs_sum
       end if

       if (rat > ZERO) then
          ratexp = ONE
       else
          ratexp = exp(rat)
       end if

       ! Compare the ratio to a random number
       if (ratexp >= ranList(i)) then
          ! accept
          det_up = new_up
          det_dn = new_dn
          accept = accept + 1
       else
          ! reject
          ! recover the old values
          Hub%G_up%G = Hub%G_up%tmp
          if (Hub%comp_dn) then
             Hub%G_dn%G = Hub%G_dn%tmp
          end if
          Hub%G_up%sgn = copy_sgn_up
          Hub%G_dn%sgn = copy_sgn_dn

          do j = 1, L
             tmp = -Hub%HSF(site(i),j)
             Hub%HSF (site(i),j) = tmp
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, j)
#endif
             Hub%V_up(site(i),j) = Hub%explook( tmp, map(site(i)))

             if (.not. Hub%neg_u) then
#if defined(DQMC_ASQRD)
                call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, j)
#endif
                Hub%V_dn(site(i),j) = Hub%explook(-tmp, map(site(i)))
             end if
          end do
       end if
    end do

    !Update determinant value
    Hub%G_up%det = det_up
    if (Hub%comp_dn) then
       Hub%G_dn%det = det_dn
    endif

    ! update G's counter
    Hub%G_up%wps = Hub%G_up%nWrap
    Hub%G_dn%wps = Hub%G_dn%nWrap

    ! update accept and reject counts
    Hub%nAcceptGlobal = Hub%nAcceptGlobal + accept
    Hub%nRejectGlobal = Hub%nRejectGlobal + (numTry-accept)

  end subroutine DQMC_Hub_Sweep2

  !---------------------------------------------------------------------!
  ! sweep for continuous HSF
  ! --------------------------------------------------------------------!

  subroutine DQMC_Hub_Sweep_Cont(Hub, nMeas0)
    !
    ! Purpose
    ! =======
    !   This subroutine performs the DQMC sweep, which is consisted of 
    !   four steps. (See [1] for more details.)
    !
    !      1. Swap the slice of G and recompute G if necessary.
    !      2. Metropolis Algorithm
    !      3. Update the model and perform physical measurement.
    !      4. Adjust parameters.
    !
    !   The first three steps are within a big loop, which run
    !   through each time slice of G. The major part is the second
    !   step, which is explained below.
    !
    !      1. Try the new configuration by single spin-flip sampling 
    !         at site j at time slice i.
    !      2. Compute the probability of this new configuration.
    !         
    !             p =  r/(1+gamma*r)    if r < 1
    !             p =  r/(gamma+r)      if r >=1
    !        
    !         where r is the ratio of determinants of Green's function
    !         of spin up and spin down.
    !      3. If p > ran, a uniform random number in [0,1], then change
    !         the configuration and update the Green's function.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout),target  :: Hub        ! Hubbard model
    integer, intent(in)                  :: nMeas0     ! Duration of measurement

    ! ... paremeters ...
    integer, parameter  :: DQMC_CHECK_ITER  = 10000
    integer, parameter  :: DQMC_ADJUST      = 100
    real(wp), parameter :: DQMC_ACC_UP      = 0.52_wp
    real(wp), parameter :: DQMC_ACC_LO      = 0.48_wp

    ! ... local scalar ...

    integer  :: i, j, k, n, L         ! Loop iterator
    integer  :: cnt                   ! Counter for measurement
    integer  :: accept, reject        ! Counter for accept in Met-alg
    real(wp) :: accrat
    real(wp) :: alpha_up, alpha_dn    ! Change of configuration
    real(wp) :: p, randn              ! Probability of changing
    real(wp) :: r_up, r_dn, r         ! Ratio of determinant
    real(wp) :: gjj_up, gjj_dn        ! (j,j) element of G_up or G_dn
    
    ! To speed up the computation
    real(wp), pointer :: G_up(:,:) 
    real(wp), pointer :: G_dn(:,:) 
    real(wp), pointer :: U_up(:,:) 
    real(wp), pointer :: U_dn(:,:) 
    real(wp), pointer :: W_up(:,:) 
    real(wp), pointer :: W_dn(:,:) 
    real(wp), pointer :: V_up(:,:) 
    real(wp), pointer :: V_dn(:,:) 
    integer,  pointer :: blksz_up  
    integer,  pointer :: blksz_dn  
    real(wp), pointer :: sgn_up    
    real(wp), pointer :: sgn_dn    

    real(wp), pointer :: ranlist(:)
    real(wp), pointer :: CHSF(:,:) 
    real(wp), pointer :: lambda(:) 
    integer,  pointer :: map(:)    
    real(wp)          :: gamma, edx, delta1, dx, dE, dPE, dKE
    logical           :: comp_dn, neg_u

    ! ... Executable ...

    !=====================! 
    ! Step 0: Setup alias !
    !=====================!

    G_up     => Hub%G_up%G
    U_up     => Hub%G_up%U
    W_up     => Hub%G_up%W
    V_up     => Hub%G_up%V
    blksz_up => Hub%G_up%blksz
    sgn_up   => Hub%G_up%sgn

    G_dn     => Hub%G_dn%G
    U_dn     => Hub%G_dn%U
    W_dn     => Hub%G_dn%W
    V_dn     => Hub%G_dn%V
    blksz_dn => Hub%G_dn%blksz
    sgn_dn   => Hub%G_dn%sgn

    ranlist  => Hub%WS%R7
    CHSF     => Hub%CHSF
    map      => Hub%S%map
    lambda   => Hub%lambda

    comp_dn  =  Hub%comp_dn
    neg_u    =  Hub%neg_u
    gamma    =  Hub%gamma
    n        =  Hub%n    
    cnt      =  nMeas0
    L        =  Hub%L
    delta1   =  Hub%delta1

    ! This is a reminder:
    ! comp_dn = .false. if
    !   1) U <= 0
    !   2) U > 0 and mu = 0 and bipartitle lattice
    !   3) Holstein model

    ! Loop over time slices
    do i = 1, L
       !==============================! 
       ! Step 1: Swap the slice of G  !
       !==============================!
       call DQMC_GetG(i, Hub%G_up, Hub%SB_up)
       if (comp_dn) then
          call DQMC_GetG(i, Hub%G_dn, Hub%SB_dn)
       else
          sgn_dn = sgn_up
       end if

       !==============================!
       ! Step 2: Metropolis Algorithm !
       !==============================!
       accept = 0
       reject = 0

       call ran0(2*n, ranlist, Hub%seed)

       ! Loop over lattice sites
       do j = Hub%n_start, Hub%n_end
          ! propose a new move
          dx = delta1*(ranlist(j+n)-HALF)

          ! Hubbard model mode
          if ( SimType .eq. Hubbard_model ) then
             ! Remember positive and negative U use different HS transfoamtaions
             if (neg_u) then
               alpha_up = exp( lambda(map(j))*dx ) - ONE
               alpha_dn = alpha_up
             else
               alpha_up = exp( lambda(map(j))*dx ) - ONE
               alpha_dn = exp(-lambda(map(j))*dx ) - ONE
               !alpha_dn = ONE/edx - ONE
             end if

             ! Compute the ratio of new and old partition functions
             gjj_up = DQMC_Gfun_Getjj(n, j, blksz_up, G_up, U_up, W_up)
             r_up   = ONE + (ONE - gjj_up)*alpha_up

             ! Need to compute r_dn explicitly if not on the half-filled bipartite lattices
             if (comp_dn) then
                gjj_dn = DQMC_Gfun_Getjj(n, j, blksz_dn, G_dn, U_dn, W_dn)
                r_dn = ONE + (ONE - gjj_dn)*alpha_dn
             ! Spin up and down ratio are identical for negative U
             else if (neg_u) then
                r_dn = r_up
             ! For positive U and on half-filled bipartite lattices, r_dn can be computed from gjj_up
             else
                r_dn = ONE + gjj_up*alpha_dn
             end if
            
             ! dE = [(x+dx)^2-x^2]/2 = x*dx + dx*dx/2
             dE = CHSF(j,i)*dx + 0.5_wp*dx*dx 

          ! Holstein model
          else if ( SimType .eq. Holstein_model) then
             alpha_up = exp( lambda(map(j))*dx ) - ONE
             gjj_up = DQMC_Gfun_Getjj(n, j, blksz_up, G_up, U_up, W_up)
             r_up   = ONE + (ONE - gjj_up)*alpha_up

             alpha_dn = alpha_up
             r_dn     = r_up
            
             ! phonon potential energy difference
             dPE = 2.0_wp*CHSF(j,i)*dx + dx*dx
            
             ! phonon kinetic energy difference
             if ( i .eq. 1 ) then
               !dKE = 2.0_wp*( 2.0_wp*CHSF(j,i)*dx + dx*dx - dx*( CHSF(j, L)   + CHSF(j,i+1) ) )
               dKE = 2.0_wp*( dPE - dx*( CHSF(j,  L) + CHSF(j,i+1) ) )
             elseif ( i .eq. L ) then
               !dKE = 2.0_wp*( 2.0_wp*CHSF(j,i)*dx + dx*dx - dx*( CHSF(j, i-1) + CHSF(j,1)   ) )
               dKE = 2.0_wp*( dPE - dx*( CHSF(j,i-1) + CHSF(j,  1) ) )
             else
               !dKE = 2.0_wp*( 2.0_wp*CHSF(j,i)*dx + dx*dx - dx*( CHSF(j, i-1) + CHSF(j,i+1) ) )
               dKE = 2.0_wp*( dPE - dx*( CHSF(j,i-1) + CHSF(j,i+1) ) )
             end if
             dE = dKE*norm_phonon(1) + dPE*norm_phonon(2)
          end if

          r  = abs(r_up*r_dn)*exp(-dE)

          ! Compute the probability
          if(r <= ONE) then
             p = r/(ONE+gamma*r)
          else
             p = r/(gamma+r)
          end if

          randn = ranlist(j)

          ! Accept 
          if (p > randn) then
             accept = accept + 1

             if(r_up < ZERO) sgn_up = -sgn_up
             if(r_dn < ZERO) sgn_dn = -sgn_dn

             CHSF(j,i) = CHSF(j,i) + dx

             ! Update G_up
             call DQMC_UpdateG(j, alpha_up/r_up, Hub%G_up)
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, i)
#endif
             V_up(j,i) = V_up(j,i)*(alpha_up + ONE)
             Hub%G_up%nModify = i
             Hub%G_up%det = Hub%G_up%det - log(abs(r_up))

             ! Update G_dn when it is necessary
             if (comp_dn) then
                call DQMC_UpdateG(j,  alpha_dn/r_dn, Hub%G_dn)
                Hub%G_dn%det = Hub%G_dn%det - log(abs(r_dn))
             end if

             if (.not. neg_u) then
#if defined(DQMC_ASQRD)
               call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, i)
#endif
               V_dn(j,i) = V_dn(j,i)*(alpha_dn + ONE)
               Hub%G_dn%det = Hub%G_dn%det - log(abs(r_dn))
             end if
             Hub%G_dn%nModify = i

          else
            ! If reject, advance the counter then move on.          
            reject = reject + 1
          end if
       end do

       !============================!
       ! Step 3: Update and Measure !
       !============================!
       ! update G_up/G_dn if there are some updates not applied.
      
       call DQMC_ApplyUpdate(Hub%G_up, forced = .true.)
       if (comp_dn) then
          call DQMC_ApplyUpdate(Hub%G_dn, forced = .true.)
       end if

       ! update accept and reject counts
       Hub%naccept = Hub%naccept + accept
       Hub%nreject = Hub%nreject + reject

       !cnt = cnt - 1
       !if (cnt == 0) then
       !   ! construct G_dn for mu = 0
       !   if (.not.Hub%neg_u .and. .not.Hub%comp_dn) then
       !      do k = 1,n
       !         do j = 1,n
       !            G_dn(k,j) = -Hub%S%P(k)*Hub%S%P(j)*G_up(j,k)
       !         end do
       !         G_dn(k,k) = G_dn(k,k) + ONE 
       !      end do
       !   end if
       !   ! Basic measurement
       !   !call DQMC_Phy0_Meas(Hub%n, Hub%P0, G_up, G_dn, &
       !   !     Hub%U, Hub%mu_up, Hub%mu_dn, Hub%t_up, Hub%t_dn, sgn_up, sgn_dn, Hub%S)
       !   !if (Hub%meas2) then
       !   !   ! Pair measurement
       !   !   r = sgn_up*sgn_dn
       !   !   call DQMC_Phy2_Meas(n, Hub%P2%M1, Hub%P2%M2, &
       !   !        Hub%P2, Hub%S%B, G_up, G_dn, r)
       !   !   
       !   !end if
       !   ! Reset the counter
       !   cnt = nMeas0
       !end if
       
    end do
    
    !===========================!
    ! Step 4: Adjust parameters !
    !===========================!
    if(Hub%naccept+Hub%nreject > DQMC_CHECK_ITER) then
       accrat = float(Hub%naccept)/float(Hub%naccept + Hub%nreject)
       if(accrat > DQMC_ACC_UP .or. accrat < DQMC_ACC_LO)then
          Hub%gamma = Hub%gamma + (accrat - HALF)
          Hub%gamma = dmax1(ZERO,Hub%gamma)
          Hub%gamma = dmin1(ONE, Hub%gamma)
          Hub%naccept = int(DQMC_ADJUST*accrat)
          Hub%nreject = int(DQMC_ADJUST*(ONE-accrat))
       endif
    endif

    call DQMC_UpdateWraps(Hub%G_up)
    if (comp_dn) then
      call DQMC_UpdateWraps(Hub%G_dn)
    end if

  end subroutine DQMC_Hub_Sweep_Cont

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Sweep2_Cont(Hub, numTry)
    !
    ! Purpose
    ! =======
    !   This subroutine performs the global moves of DQMC sweep, in which 
    !   all the Hub(i) on some selected sites are flipped for all slice.
    !
    !      1. Try the new configuration.
    !      2. Compute the probability of this new configuration.
    !         
    !             p =  r/(1+gamma*r)    if r < 1
    !             p =  r/(gamma+r)      if r >=1
    !        
    !         where r is the ratio of determinants of Green's function
    !         of spin up and spin down.
    !      3. If p > ran, a uniform random number in [0,1], then change
    !         the configuration and update the Green's function.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout)  :: Hub      ! Hubbard model
    integer, intent(in)           :: numTry   ! Number of Try

    ! ... Local Variables ...
    integer, pointer       :: map(:) 
    real(wp), pointer      :: CHSF(:,:) 
    integer                :: i, j, k, n, L, si, sj, tmp, nSite, accept, reject
    integer                :: siteList(Hub%n)
    integer                :: slice(numTry), site(numTry)
    real(wp)               :: ranList(2*numTry), rat, ratexp
    real(wp)               :: det_up, det_dn, new_up, new_dn
    real(wp)               :: copy_sgn_up, copy_sgn_dn, delta2, dx
    real(wp)               :: G_dn_tmp(Hub%n,Hub%n)
    real(wp)               :: E_old, E_new
    logical                :: compute_dn, neg_u

    ! ... Executable ...

    if (numTry <= 0) return

    compute_dn = Hub%comp_dn .or. .not.Hub%neg_u
    neg_u      = Hub%neg_u
    n          = Hub%n
    L          = Hub%L
    accept     = 0
    reject     = 0
    delta2     = Hub%delta2
    map        => Hub%S%map
    CHSF       => Hub%CHSF
    
    ! Compute the Green's matrix and the sign
#if defined(DQMC_ASQRD)
    call cpp_gfun_computeg(Hub%G_up%cpp_data, L, Hub%G_up%sgn, &
         Hub%G_up%G, Hub%V_up, Hub%SB_up%B%B, &
         Hub%SB_up%nOrth, det_up)
#else
    call DQMC_ComputeG(L, n, Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, &
         Hub%SB_up, Hub%G_up%pvt, .true., det_up, HUb%G_up%sxx)
#endif
    if (compute_dn) then
#if defined(DQMC_ASQRD)
       call cpp_gfun_computeg(Hub%G_dn%cpp_data, L, Hub%G_dn%sgn, &
            Hub%G_dn%G, Hub%V_dn, Hub%SB_dn%B%B, &
            Hub%SB_dn%nOrth, det_dn)
#else
       call DQMC_ComputeG(L, n, Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, &
            Hub%SB_dn, Hub%G_dn%pvt, .true., det_dn, Hub%G_dn%sxx)
#endif
    else
       det_dn = det_up
       Hub%G_dn%sgn = Hub%G_up%sgn
    end if

    ! get random numbers
    call ran0(2*numTry, ranList, Hub%seed)

    nsite = Hub%n_End - Hub%n_start + 1
    siteList(1:nSite) = Hub%n_Start+(/(i,i=0,nSite-1)/)

    ! generate sites
    do i = 1, numtry
       tmp = int(ranList(i)*nSite) + 1
       site(i) = siteList(tmp)
       ! compress the list
       do j = tmp+1, nSite
          siteList(j-1) = siteList(j) 
       end do
       nSite = nSite - 1
    end do
    
    ! generate slice
    do i = 1, numtry
       tmp = int(ranList(i+numTry)*L) + 1
       slice(i) = tmp
    end do

    call ran0(numTry, ranList, Hub%seed)

    ! Global move
    do i = 1, numTry
       si = site(i)
       !sj = slice(i)

       if (SimType .eq. Hubbard_model) then
          E_old = 0.0_wp
          E_new = 0.0_wp
          do j = 1, L
             !CHSF (si,sj) = CHSF(site(i),sj) + dx
             !CHSF (si,sj) = CHSF(si,sj) + dx

             dx = delta2*(ranlist(i)-0.5d0)
             ! compute the old Gaussian exponent
             E_old = E_old - 0.5_wp*CHSF(si,j)*CHSF(si,j)

             ! update fields at all time slices at site index si
             CHSF (si,j) = CHSF(si,j) + dx

             ! compute the new Gaussian exponent
             E_new = E_new - 0.5_wp*CHSF(si,j)*CHSF(si,j)
          
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, sj)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, sj)
#endif
             !Hub%V_up(si,sj) = exp(Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             !Hub%V_dn(si,sj) = exp(-Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             Hub%V_up(si,j) = exp(Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             if (neg_u) then
               Hub%V_dn(si,j) = Hub%V_up(si,j)
             else
               Hub%V_dn(si,j) = exp(-Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             end if
          end do

       else if (SimType .eq. Holstein_model) then
          E_old = 0.0_wp
          E_new = 0.0_wp

          !subroutine PHONON_ACTION(L, norm, x, Sb)
          call Phonon_Action(L, norm_phonon, CHSF(si,:), E_old)
          dx = delta2*(ranlist(i)-0.5d0)
          do j = 1, L
             ! update fields at all time slices at site index si
             CHSF (si,j) = CHSF(si,j) + dx

#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, sj)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, sj)
#endif
             Hub%V_up(si,j) = exp(Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             if (neg_u) then
               Hub%V_dn(si,j) = Hub%V_up(si,j)
             else
               Hub%V_dn(si,j) = exp(-Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             end if
          end do

          call Phonon_Action(L, norm_phonon, CHSF(si,:), E_new)
       end if
       
       ! Store the value of G first
       Hub%G_up%tmp = Hub%G_up%G
       if (compute_dn) then
          G_dn_tmp = Hub%G_dn%G
       end if
       copy_sgn_up = Hub%G_up%sgn
       copy_sgn_dn = Hub%G_dn%sgn

       ! Compute G with new configuration
#if defined(DQMC_ASQRD)
       call cpp_gfun_computeg(Hub%G_up%cpp_data, L, n, &
            Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, Hub%SB_up%B%B, &
            Hub%SB_up%L, Hub%SB_up%nOrth, new_up)
#else
       call DQMC_ComputeG(L, n, Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, &
            Hub%SB_up, Hub%G_up%pvt, .true., new_up, Hub%G_up%sxx)
#endif
       if (compute_dn) then
#if defined(DQMC_ASQRD)
          call cpp_gfun_computeg(Hub%G_up%cpp_data, L, n, &
               Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, Hub%SB_dn%B%B, &
               Hub%SB_dn%L, Hub%SB_dn%nOrth, new_dn)
#else
          call DQMC_ComputeG(L, n, Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, &
               Hub%SB_dn, Hub%G_dn%pvt, .true., new_dn, Hub%G_up%sxx)
#endif
       else
          new_dn = new_up
          Hub%G_dn%sgn =  Hub%G_up%sgn
       end if

       ! Compute the Det ratio
       ! NB: the determinant computed by GetG is log(abs(det(G)))
       !     Here we need log(abs(Z))= -log(abs(det(G)))
       rat = det_up + det_dn - new_up - new_dn

       ratexp = exp(rat)
       ratexp = ratexp * exp(E_new - E_old)

       ! Compare the ratio to a random number
       ! add random number
       
       if (ratexp >= ranList(i)) then    
          ! accept
          det_up = new_up
          det_dn = new_dn
          accept = accept + 1

          ! update G's counter
          !Hub%G_up%wps = Hub%G_up%nWrap
          !Hub%G_dn%wps = Hub%G_dn%nWrap
       else                  
          ! reject
          ! recover the old values
          Hub%G_up%G = Hub%G_up%tmp
          if (compute_dn) then
             Hub%G_dn%G = G_dn_tmp
          end if
          Hub%G_up%sgn = copy_sgn_up
          Hub%G_dn%sgn = copy_sgn_dn
          reject = reject + 1

          !sj = slice(i)
          do j = 1, L
             CHSF (si,j) = CHSF (si,j) - dx
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, sj)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, sj)
#endif
             !Hub%V_up(si,sj) = exp(Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             !Hub%V_dn(si,sj) = exp(-Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             Hub%V_up(si,j) = exp(Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             if (neg_u) then
               Hub%V_dn(si,j) = Hub%V_up(si,j)
             else
               Hub%V_dn(si,j) = exp(-Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             end if

          end do
       end if
    end do

    !Update determinant value
    Hub%G_up%det = det_up
    if (Hub%comp_dn) then
       Hub%G_dn%det = det_dn
    endif

    ! update G's counter
    Hub%G_up%wps = Hub%G_up%nWrap
    Hub%G_dn%wps = Hub%G_dn%nWrap

    ! update accept and reject counts
    !Hub%naccept = Hub%naccept + accept
    !Hub%nreject = Hub%nreject + (numTry-accept)
    Hub%nAcceptGlobal = Hub%nAcceptGlobal + accept
    Hub%nRejectGlobal = Hub%nRejectGlobal + (numTry-accept)

    do k = 1, L
       write (*,'(A,f6.2,1x,f6.2,1x,f6.2,1x,f6.2)') 'Sweep2', CHSF(1,k), CHSF(2,k), CHSF(3,k), CHSF(4,k)
    end do

    contains
      subroutine PHONON_ACTION(L, norm, x, Sb)
        implicit none
        integer, intent(in)     :: L
        real(wp), intent(in)    :: norm(2)
        real(wp), intent(in)    :: x(L)
        real(wp), intent(inout) :: Sb

        integer  :: iL
        real(wp) :: Ptmp, Ktmp

        Ptmp = 0.0_wp
        Ktmp = 0.0_wp
        do iL = 1, L-1
          Ptmp = Ptmp - x(iL)*x(iL)
          Ktmp = Ktmp - ( x(iL) - x(iL+1) )*( x(iL) - x(iL+1) )
        end do
        Ptmp = Ptmp - x(L)*x(L)
        Ktmp = Ktmp - ( x(L) - x(1) )*( x(L) - X(1) )
        Sb = Ktmp * norm(1) + Ptmp * norm(2)
      return
      end subroutine PHONON_ACTION
  end subroutine DQMC_Hub_Sweep2_Cont

  !---------------------------------------------------------------------!
  subroutine DQMC_Hol_Sweep3_Cont(Hub, numTry)
    !
    ! Purpose
    ! =======
    !   This subroutine performs the global moves of DQMC sweep, in which 
    !   all the Hub(i) on some selected sites are flipped for all slice.
    !
    !      1. Try the new configuration.
    !      2. Compute the probability of this new configuration.
    !         
    !             p =  r/(1+gamma*r)    if r < 1
    !             p =  r/(gamma+r)      if r >=1
    !        
    !         where r is the ratio of determinants of Green's function
    !         of spin up and spin down.
    !      3. If p > ran, a uniform random number in [0,1], then change
    !         the configuration and update the Green's function.
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout)  :: Hub      ! Hubbard model
    integer, intent(in)           :: numTry   ! Number of Try

    ! ... Local Variables ...
    integer, pointer       :: map(:) 
    real(wp), pointer      :: CHSF(:,:) 
    integer                :: i, j, k, n, L, si, sj, tmp, nSite, accept, reject
    integer                :: siteList(Hub%n)
    integer                :: slice(numTry), site(numTry)
    real(wp)               :: ranList(2*numTry), rat, ratexp
    real(wp)               :: det_up, det_dn, new_up, new_dn
    real(wp)               :: copy_sgn_up, copy_sgn_dn, delta2, dx
    real(wp)               :: G_dn_tmp(Hub%n,Hub%n)
    real(wp)               :: E_old, E_new
    logical                :: compute_dn, neg_u

    ! ... Executable ...

    if (numTry <= 0) return

    compute_dn = Hub%comp_dn .or. .not.Hub%neg_u
    neg_u      = Hub%neg_u
    n          = Hub%n
    L          = Hub%L
    accept     = 0
    reject     = 0
    delta2     = Hub%delta2
    map        => Hub%S%map
    CHSF       => Hub%CHSF
    
    ! Compute the Green's matrix and the sign
#if defined(DQMC_ASQRD)
    call cpp_gfun_computeg(Hub%G_up%cpp_data, L, Hub%G_up%sgn, &
         Hub%G_up%G, Hub%V_up, Hub%SB_up%B%B, &
         Hub%SB_up%nOrth, det_up)
#else
    call DQMC_ComputeG(L, n, Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, &
         Hub%SB_up, Hub%G_up%pvt, .true., det_up, HUb%G_up%sxx)
#endif
    if (compute_dn) then
#if defined(DQMC_ASQRD)
       call cpp_gfun_computeg(Hub%G_dn%cpp_data, L, Hub%G_dn%sgn, &
            Hub%G_dn%G, Hub%V_dn, Hub%SB_dn%B%B, &
            Hub%SB_dn%nOrth, det_dn)
#else
       call DQMC_ComputeG(L, n, Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, &
            Hub%SB_dn, Hub%G_dn%pvt, .true., det_dn, Hub%G_dn%sxx)
#endif
    else
       det_dn = det_up
       Hub%G_dn%sgn = Hub%G_up%sgn
    end if

    ! get random numbers
    call ran0(2*numTry, ranList, Hub%seed)

    nsite = Hub%n_End - Hub%n_start + 1
    siteList(1:nSite) = Hub%n_Start+(/(i,i=0,nSite-1)/)

    ! generate sites
    do i = 1, numtry
       tmp = int(ranList(i)*nSite) + 1
       site(i) = siteList(tmp)
       ! compress the list
       do j = tmp+1, nSite
          siteList(j-1) = siteList(j) 
       end do
       nSite = nSite - 1
    end do
    
    ! generate slice
    do i = 1, numtry
       tmp = int(ranList(i+numTry)*L) + 1
       slice(i) = tmp
    end do

    call ran0(numTry, ranList, Hub%seed)

    ! Global move
    do i = 1, numTry
       si = site(i)
       !sj = slice(i)

       if (SimType .eq. Hubbard_model) then
          write(*,*) ' No Sweep3 for Hubbard model.'

       else if (SimType .eq. Holstein_model) then
          E_old = 0.0_wp
          E_new = 0.0_wp

          !subroutine PHONON_ACTION(L, norm, x, Sb)
          call Phonon_Action(L, norm_phonon, CHSF(si,:), E_old)
          dx = 2.0d0*Hub%lambda(map(si))/((Hub%omega)**2)/Hub%dtau
          do j = 1, L
             ! update fields at all time slices at site index si
             CHSF(si,j) = -CHSF(si,j) + dx

#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, sj)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, sj)
#endif
             Hub%V_up(si,j) = exp(Hub%lambda(map(si))*CHSF(si,j))
             if (neg_u) then
               Hub%V_dn(si,j) = Hub%V_up(si,j)
             else
               Hub%V_dn(si,j) = exp(-Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             end if
          end do

          call Phonon_Action(L, norm_phonon, CHSF(si,:), E_new)
       end if
       
       ! Store the value of G first
       Hub%G_up%tmp = Hub%G_up%G
       if (compute_dn) then
          G_dn_tmp = Hub%G_dn%G
       end if
       copy_sgn_up = Hub%G_up%sgn
       copy_sgn_dn = Hub%G_dn%sgn

       ! Compute G with new configuration
#if defined(DQMC_ASQRD)
       call cpp_gfun_computeg(Hub%G_up%cpp_data, L, n, &
            Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, Hub%SB_up%B%B, &
            Hub%SB_up%L, Hub%SB_up%nOrth, new_up)
#else
       call DQMC_ComputeG(L, n, Hub%G_up%sgn, Hub%G_up%G, Hub%V_up, &
            Hub%SB_up, Hub%G_up%pvt, .true., new_up, Hub%G_up%sxx)
#endif
       if (compute_dn) then
#if defined(DQMC_ASQRD)
          call cpp_gfun_computeg(Hub%G_up%cpp_data, L, n, &
               Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, Hub%SB_dn%B%B, &
               Hub%SB_dn%L, Hub%SB_dn%nOrth, new_dn)
#else
          call DQMC_ComputeG(L, n, Hub%G_dn%sgn, Hub%G_dn%G, Hub%V_dn, &
               Hub%SB_dn, Hub%G_dn%pvt, .true., new_dn, Hub%G_up%sxx)
#endif
       else
          new_dn = new_up
          Hub%G_dn%sgn =  Hub%G_up%sgn
       end if

       ! Compute the Det ratio
       ! NB: the determinant computed by GetG is log(abs(det(G)))
       !     Here we need log(abs(Z))= -log(abs(det(G)))
       rat = det_up + det_dn - new_up - new_dn

!       if (rat > ZERO) then
!          ratexp = ONE
!       else
          ratexp = exp(rat)
!       end if
       ratexp = ratexp * exp(E_new - E_old)

       ! Compare the ratio to a random number
       ! add random number
       
       if (ratexp >= ranList(i)) then    
          ! accept
          det_up = new_up
          det_dn = new_dn
          accept = accept + 1

          ! update G's counter
          !Hub%G_up%wps = Hub%G_up%nWrap
          !Hub%G_dn%wps = Hub%G_dn%nWrap
       else                  
          ! reject
          ! recover the old values
          Hub%G_up%G = Hub%G_up%tmp
          if (compute_dn) then
             Hub%G_dn%G = G_dn_tmp
          end if
          Hub%G_up%sgn = copy_sgn_up
          Hub%G_dn%sgn = copy_sgn_dn
          reject = reject + 1

          !sj = slice(i)
          do j = 1, L
!             dx = 2.0d0*Hub%lambda(Hub%S%map(si))/(Hub%omega**2)/Hub%dtau
             CHSF(si,j) = -CHSF(si,j) + dx
#if defined(DQMC_ASQRD)
             call cpp_gfun_invalid_cache(Hub%G_up%cpp_data, sj)
             call cpp_gfun_invalid_cache(Hub%G_dn%cpp_data, sj)
#endif
             !Hub%V_up(si,sj) = exp(Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             !Hub%V_dn(si,sj) = exp(-Hub%lambda(Hub%S%map(sj))*CHSF(si,sj))
             Hub%V_up(si,j) = exp(Hub%lambda(map(si))*CHSF(si,j))
             if (neg_u) then
               Hub%V_dn(si,j) = Hub%V_up(si,j)
             else
               Hub%V_dn(si,j) = exp(-Hub%lambda(Hub%S%map(si))*CHSF(si,j))
             end if

          end do
       end if
    end do

    !Update determinant value
    Hub%G_up%det = det_up
    if (Hub%comp_dn) then
       Hub%G_dn%det = det_dn
    endif

    ! update G's counter
    Hub%G_up%wps = Hub%G_up%nWrap
    Hub%G_dn%wps = Hub%G_dn%nWrap

    ! update accept and reject counts
    !Hub%naccept = Hub%naccept + accept
    !Hub%nreject = Hub%nreject + (numTry-accept)
    Hub%nAcceptGlobal2 = Hub%nAcceptGlobal2 + accept
    Hub%nRejectGlobal2 = Hub%nRejectGlobal2 + (numTry-accept)


    do k = 1, L
       write (*,'(A,f6.2,1x,f6.2,1x,f6.2,1x,f6.2)') 'Sweep3', CHSF(1,k), CHSF(2,k), CHSF(3,k), CHSF(4,k)
    end do

    contains
      subroutine PHONON_ACTION(L, norm, x, Sb)
        implicit none
        integer, intent(in)     :: L
        real(wp), intent(in)    :: norm(2)
        real(wp), intent(in)    :: x(L)
        real(wp), intent(inout) :: Sb

        integer  :: iL
        real(wp) :: Ptmp, Ktmp

        Ptmp = 0.0_wp
        Ktmp = 0.0_wp
        do iL = 1, L-1
          Ptmp = Ptmp - x(iL)*x(iL)
          Ktmp = Ktmp - ( x(iL) - x(iL+1) )*( x(iL) - x(iL+1) )
        end do
        Ptmp = Ptmp - x(L)*x(L)
        Ktmp = Ktmp - ( x(L) - x(1) )*( x(L) - x(1) )
        Sb = Ktmp * norm(1) + Ptmp * norm(2)
      return
      end subroutine PHONON_ACTION

  end subroutine DQMC_Hol_Sweep3_Cont

  !---------------------------------------------------------------------!

  subroutine DQMC_Hub_Run(Hub,Info)
    !
    ! Purpose
    ! =======
    !   This subroutine is the main subroutine for DQMC.
    !   There are four major wroks
    !
    !      1. Compute Green function.
    !      2. Perform warmup sweep.
    !      3. Perform actual sweep.
    !      4. Analyze the measurement. (see DQMC_Phy0)
    !
    ! Arguments
    ! =========
    !
    !   Info == 1: Print runtime information
    !   Info == 0: Silent mode
    !
    type(Hubbard), intent(inout) :: Hub    ! Hubbard model

    ! ... local scalar ...
    integer  :: i, j, nIter, nBin, Info

    ! ... Executable ...

    ! Warmup sweep
    do i = 1, Hub%nWarm
       if (Info==1 .and. mod(i,10)==0) write(*,'(A,i6,1x,i3)')' Warmup Sweep, nwrap  : ', i, Hub%G_up%nwrap
       ! The second parameter means no measurement should be made.
       call DQMC_Hub_Sweep(Hub, NO_MEAS0)
       call DQMC_Hub_Sweep2(Hub, Hub%nTry)
    end do
 
    ! We divide all the measurement into nBin,
    ! each having nPass/nBin pass.
    nBin   = Hub%P0%nBin 
    nIter  = Hub%nPass/nBin
    do i = 1, nBin
       do j = 1, nIter
          call DQMC_Hub_Sweep(Hub, Hub%nMeas)
          call DQMC_Hub_Sweep2(Hub, Hub%nTry)
       end do

       ! Accumulate results for each bin
       if (Info==1) write(*,'(a,2i6)') ' Measurement Sweep, bin, iter : ', i, j
       call DQMC_Phy0_Avg(Hub%P0)
       if (Hub%meas2) then
          if(Hub%P2%diagonalize)then
            call DQMC_Phy2_Avg(Hub%P2, Hub%S)
          else
            call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
          endif
       end if
    end do

    ! Get average result
    call DQMC_Phy0_GetErr(Hub%P0)
    if (Hub%meas2) then
       call DQMC_Phy2_GetErr(Hub%P2)
    end if

  end subroutine DQMC_Hub_Run
  
  !--------------------------------------------------------------------!

  subroutine DQMC_Hub_FullMeas(Hub, nnb, A_up, A_dn, sgn_up, sgn_dn)

     !type(Hubbard), target, intent(inout) :: Hub
     type(Hubbard), intent(inout) :: Hub
     integer,  intent(in)         :: nnb
     real(wp), intent(in)         :: A_up(nnb, nnb)
     real(wp), intent(in)         :: A_dn(nnb, nnb)
     real(wp), intent(in)         :: sgn_up
     real(wp), intent(in)         :: sgn_dn

     integer :: it, i, j, nb

     ! Extra space for dtau^2-correct estimate
     type(G_fun), target :: G_up_local, G_dn_local
     real(wp), pointer   :: G_up(:,:) 
     real(wp), pointer   :: G_dn(:,:) 
     
     !Duplicate the Green's function 
     call DQMC_Gfun_Duplicate(G_up_local, Hub%G_up)
     if (.not.Hub%neg_u .or. Hub%comp_dn) then
        call DQMC_Gfun_Duplicate(G_dn_local, Hub%G_dn)
     else
        call DQMC_Gfun_clone(G_dn_local, G_up_local)
     endif
     G_up   => G_up_local%G
     G_dn   => G_dn_local%G

     !write(*,*) "Hub%S%checklist(PHASE)=",Hub%S%checklist(STRUCT_PHASE)

     nb = nnb / Hub%n
     !Perform static measurement on all stored time slices
     do it = 0, nb-1
        !Load the diagonal of Aup/Adn in G_up/G_dn
        i = it*Hub%n + 1 
        j = i + Hub%n - 1
        G_up = A_up(i:j, i:j)

        ! Modified in order to take into account U < 0 model on non-bipartite lattices.
        ! Before the modification, the code calls DQMC_Gfun_CopyUp() which assumes bipartite lattice.
        ! This causes segmentation fault on, for example, the triagular lattice which lacks particle-hole
        ! symmetry. As a result, S%P is not defined. We fix it by separating U < 0 model from the rest
        ! of the if-else statements, and calling DQMC_Gfun_clone().
        if ( Hub%comp_dn ) then
           ! Get G_dn directly when :
           !     1) U > 0, mu_up .neq. 0, or mu_dn .neq. 0. 
           !     2) U > 0, mu_up = mu_dn = 0, but "PHASE" S%P is not defined. E.g. non-bipartite lattice or S%P is simply lacking.
           !     3) U > 0, mu_up = mu_dn = 0, but t_up .neq. t_dn.
           !     4) U < 0, but mu_up .neq. mu_dn, or t_up .neq. t_dn.
           !     5) U = 0, but mu_up .neq. mu_dn, or t_up .neq. t_dn.
           G_dn = A_dn(i:j, i:j)
        else if ( Hub%neg_u ) then
           ! G_dn and G_up are identical when :
           !     1) U < 0, mu_up = mu_dn, and t_up = t_dn
           !     2) U = 0, mu_up = mu_dn, and t_up = t_dn
           call DQMC_Gfun_Clone(G_dn_local, G_up_local)
        else 
           ! Note that here we are left with the last condition: (.not.neg_u) and (.not.comp_dn).
           ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
           ! So we can safely call DQMC_Gfun_CopyUp() which assumes particle-hole symmetry.
           ! Use particle-hole symmetry to get G_dn when :
           !     1) U > 0, at half-filling, t_up = t_dn, and "PHASE" S%P is defined.
           call DQMC_Gfun_CopyUp(G_dn_local, G_up_local, Hub%S%P)
        endif

        call DQMC_GetG_2nd_order(G_up_local, Hub%B_up)
        if (Hub%comp_dn .or. .not.Hub%neg_u) then
           call DQMC_GetG_2nd_order(G_dn_local, Hub%B_dn)
        endif

        call DQMC_Phy0_Meas(Hub%n, Hub%P0, G_up_local%GS, G_dn_local%GS, Hub%U, &
               Hub%mu_up, Hub%mu_dn, Hub%t_up, Hub%t_dn, sgn_up, sgn_dn, Hub%S)
        if (SimType .eq. Holstein_model) then
          call DQMC_Phy0_Meas_Holstein(Hub%n, Hub%L, Hub%P0, G_up_local%GS, G_dn_local%GS, &
                 Hub%U, sgn_up, sgn_dn, Hub%S, Hub%CHSF, norm_phonon, Hub%dtau)
        end if

        if (Hub%meas2) then
           ! Pair measurement
           call DQMC_Phy2_Meas(Hub%n, Hub%P2%M1, Hub%P2%M2, Hub%P2, Hub%S%B, &
              G_up_local%GS, G_dn_local%GS, sgn_up*sgn_dn)
        end if

     enddo

     call DQMC_Gfun_Free(G_up_local)
     call DQMC_Gfun_Free(G_dn_local)

  end subroutine 

  !-------------------------------------------------------------------!

  subroutine DQMC_Hub_Meas(Hub, slice)

     !type(Hubbard), target, intent(inout) :: Hub
     type(Hubbard), intent(inout) :: Hub
     integer, intent(inout)       :: slice

     type(G_fun), target :: G_up_local, G_dn_local
     real(wp), pointer   :: G_up(:,:) 
     real(wp), pointer   :: G_dn(:,:) 
     real(wp), pointer   :: sgn_up    
     real(wp), pointer   :: sgn_dn    
     real(wp)            :: randn(1)
 
     ! Warning: if slice = 0, DQMC_Hub_Meas() would return meaningless results.
     if (slice <= 0 .or. slice > Hub%L) then
       write(*,*) " In subroutine DQMC_Hub_Meas(Hub, slice), the argument 'slice' is out of bound."
       write(*,*) " It will now be reset randomly."
       call ran0(1, randn, Hub%seed)
       slice = ceiling(randn(1)*Hub%L)
       write(*,*) " New time slice index is", slice
     end if
 
     !Duplicate the Green's function
     call DQMC_Gfun_Duplicate(G_up_local, Hub%G_up)
   
     if (.not.Hub%neg_u .or. Hub%comp_dn) then
        call DQMC_Gfun_Duplicate(G_dn_local, Hub%G_dn)
     else
        call DQMC_Gfun_clone(G_dn_local, G_up_local)
     endif
     G_up   => G_up_local%G
     G_dn   => G_dn_local%G
     sgn_up => G_up_local%sgn
     sgn_dn => G_dn_local%sgn

     !Recompute G from scratch
     G_up_local%ilb = -1       
     call DQMC_GetG(slice, G_up_local, Hub%SB_up)
     if ( Hub%comp_dn ) then
        G_dn_local%ilb = -1       
        call DQMC_GetG(slice, G_dn_local, Hub%SB_dn)
     elseif ( Hub%neg_u ) then
        sgn_dn = sgn_up
     else
        ! Note that here we have (.not.neg_u) and (.not.comp_dn).
        ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
        ! So we can safely call DQMC_Gfun_CopyUp() which uses particle-hole symmetry.
        call DQMC_Gfun_CopyUp(G_dn_local, G_up_local, Hub%S%P)
     endif

     !Get G correct to 2nd order
     call DQMC_GetG_2nd_order(G_up_local, Hub%B_up)
     if (Hub%comp_dn .or. .not.Hub%neg_u) then
        call DQMC_GetG_2nd_order(G_dn_local, Hub%B_dn)
     endif
     
     ! Basic measurement
     call DQMC_Phy0_Meas(Hub%n, Hub%P0, G_up_local%GS, G_dn_local%GS, Hub%U, &
            Hub%mu_up, Hub%mu_dn, Hub%t_up, Hub%t_dn, sgn_up, sgn_dn, Hub%S)
     if (SimType .eq. Holstein_model) then
       call DQMC_Phy0_Meas_Holstein(Hub%n, Hub%L, Hub%P0, G_up_local%GS, G_dn_local%GS, &
              Hub%U, sgn_up, sgn_dn, Hub%S, Hub%CHSF, norm_phonon, Hub%dtau)
     end if


     if (Hub%meas2) then
        ! Pair measurement
        call DQMC_Phy2_Meas(Hub%n, Hub%P2%M1, Hub%P2%M2, Hub%P2, Hub%S%B, &
           G_up_local%GS, G_dn_local%GS, sgn_up*sgn_dn)
     end if

     call DQMC_Gfun_Free(G_up_local)
     call DQMC_Gfun_Free(G_dn_local)

  end subroutine DQMC_Hub_Meas

  !-------------------------------------------------------------------!

  subroutine DQMC_Hub_Init_Vmat(Hub)
    !
    ! Purpose
    ! =======
    !    The element of V(i) is either exp(nu) or exp(-nu)
    !    where nu = acosh(exp(U*dtau/2)). (see reference [1].) 
    !    The values of exp(nu) and exp(-nu) are stored in a lookup 
    !    table explook.  The decision of wheather V(i,j) is exp(nu) 
    !    or exp(-nu) is given by the list hub, which is a list 
    !    or random +1 and -1. Matrix V for spin up and down have
    !    opposite selection decision.
    ! 
    !
    ! Arguments
    ! =========
    !
     type(Hubbard), intent(inout) :: Hub

    ! ... Local variables ...
     integer :: i, j
     real(wp) :: temp, temp2

    ! ... Executable ...
     if (.not.associated(Hub%V_up)) then
        allocate(Hub%V_up(Hub%n,Hub%L))
     endif

     ! This fix is required to make test program in /EXAMPLE/test work.
     ! Without the fix, the array size of Hub%V_up and Hub%V_dn would be 1 rather then n * L.
     ! This causes segmentation fault when running the test program. 
     if (size(Hub%V_up) /= Hub%n*Hub%L) then
       allocate(Hub%V_up(Hub%n,Hub%L))
     end if

     ! discrete fields
     if (Hub%HSFtype .eq. HSF_DISC) then
        do i = 1, Hub%L
           do j = 1, Hub%n
              Hub%V_up(j,i) = Hub%explook(Hub%HSF(j,i), Hub%S%map(j))
           end do
        end do
        
        if (Hub%neg_u) then
           if (.not.associated(Hub%V_dn)) then
              Hub%V_dn => Hub%V_up
           endif
        else
           if (.not.associated(Hub%V_dn)) then
              allocate(Hub%V_dn(Hub%n,Hub%L))
           endif
           ! This fix is required to make test program in /EXAMPLE/test work.
           ! Without the fix, the array size of Hub%V_up and Hub%V_dn would be 1 rather then n * L.
           ! This causes segmentation fault when running the test program. 
           if (size(Hub%V_dn) /= Hub%n*Hub%L) then
             allocate(Hub%V_dn(Hub%n,Hub%L))
           end if

           do i = 1, Hub%L
              do j = 1, Hub%n
                 Hub%V_dn(j,i) = Hub%explook(-Hub%HSF(j,i), Hub%S%map(j))
              end do
           end do
        end if

     ! continuous fields
     else if (Hub%HSFtype .eq. HSF_CONT) then

        if (SimType .eq. Hubbard_model) then
           if (Hub%neg_u) then
              do i = 1, Hub%L
                 do j = 1, Hub%n
                    temp  = Hub%lambda(Hub%S%map(j))*Hub%CHSF(j,i)
                    temp2 = Hub%lambda(Hub%S%map(j))*Hub%lambda(Hub%S%map(j))
                    Hub%V_up(j,i) = exp(temp) ! + temp2)
                 end do
              end do

              if (.not.associated(Hub%V_dn)) Hub%V_dn => Hub%V_up
           else
              if (.not.associated(Hub%V_dn)) allocate(Hub%V_dn(Hub%n,Hub%L))
              ! The following is required to make test program in /EXAMPLE/test work.
              ! Without the fix, the array size of Hub%V_up and Hub%V_dn would be 1 rather then n * L.
              ! This causes segmentation fault when running the test program. 
              if (size(Hub%V_dn) /= Hub%n*Hub%L) allocate(Hub%V_dn(Hub%n,Hub%L))

              do i = 1, Hub%L
                 do j = 1, Hub%n
                    temp = Hub%lambda(Hub%S%map(j))*Hub%CHSF(j,i)
                    Hub%V_up(j,i) = exp( temp)
                    Hub%V_dn(j,i) = exp(-temp)
                 end do
              end do
           end if

        else if (SimType .eq. Holstein_model) then
           do i = 1, Hub%L
              do j = 1, Hub%n
                 temp = Hub%lambda(Hub%S%map(j))*Hub%CHSF(j,i)
                 Hub%V_up(j,i) = exp(temp)
              end do
           end do
           if (.not.associated(Hub%V_dn)) Hub%V_dn => Hub%V_up
        end if
 
     end if
  end subroutine DQMC_Hub_Init_Vmat

  !-------------------------------------------------------------------!

end module DQMC_Hubbard
