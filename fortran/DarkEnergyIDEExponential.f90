    ! Coupled CDM and canonical scalar field. The dimensionless field
    ! is phi = sqrt(8*pi*G)*psi, with m_c(phi) proportional to exp(-beta*phi).
    module DarkEnergyIDEExponential
    use precision
    use constants, only: c
    use classes, only: TCAMBdata, TPythonInterfacedClass
    use config, only: GlobalError, error_darkenergy, global_error_flag
    use DarkEnergyInterface, only: TDarkEnergyModel
    use results, only: CAMBdata
    use MassiveNu, only: ThermalNuBackground
    use RungeKuttaDP45Module, only: RungeKuttaDP45Settings, TClassRungeKuttaDP45
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    real(dl), parameter :: grho100 = 3._dl*(100._dl*1000._dl/c)**2
    integer, parameter :: NumEqs = 2

    type, extends(TDarkEnergyModel) :: TIDEExponential
        ! Public inputs. Keep these in the same order as camb/dark_energy.py.
        real(dl) :: V0_over_rho100 = 0.7_dl
        real(dl) :: omc_early = 0.12_dl
        real(dl) :: beta = 0._dl
        real(dl) :: lam = 0._dl
        real(dl) :: a_init = 1e-10_dl
        real(dl) :: integrate_tol = 1e-10_dl
        integer :: history_points = 2049
        ! A uniform log(a) table and cubic Hermite data for the forward solution.
        real(dl), private :: log_a_init, dlog_a
        real(dl), private :: grhog0, grhornomass0, grhob0
        real(dl), private, allocatable :: nu_grho(:), nu_mass(:)
        real(dl), private, allocatable :: phi_table(:), u_table(:), du_table(:)
    contains
    procedure :: RequiresResolution => TIDEExponential_RequiresResolution
    procedure :: ResolveBackground => TIDEExponential_ResolveBackground
    procedure :: CDMGrhoa4 => TIDEExponential_CDMGrhoa4
    procedure :: PhiAtA => TIDEExponential_PhiAtA
    procedure :: InitialScaleFactor => TIDEExponential_InitialScaleFactor
    procedure :: CouplingBeta => TIDEExponential_CouplingBeta
    procedure :: ScalarBackgroundAtA => TIDEExponential_ScalarBackgroundAtA
    procedure :: CoupledFieldSource => TIDEExponential_CoupledFieldSource
    procedure :: ExpansionLogSlope => TIDEExponential_ExpansionLogSlope
    procedure :: RhoCOverRho100 => TIDEExponential_RhoCOverRho100
    procedure :: Init => TIDEExponential_Init
    procedure :: BackgroundDensityAndPressure => TIDEExponential_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TIDEExponential_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TIDEExponential_PerturbationEvolve
    procedure, nopass :: PythonClass => TIDEExponential_PythonClass
    procedure, nopass :: SelfPointer => TIDEExponential_SelfPointer
    procedure, private :: EvolveLogA => TIDEExponential_EvolveLogA
    procedure, private :: TermsAtA => TIDEExponential_TermsAtA
    procedure, private :: FieldAtA => TIDEExponential_FieldAtA
    end type TIDEExponential

    public TIDEExponential
    procedure(TClassRungeKuttaDP45) :: RungeKuttaDP45

    contains

    logical function TIDEExponential_RequiresResolution(this)
    class(TIDEExponential), intent(in) :: this

    TIDEExponential_RequiresResolution = .true.

    end function TIDEExponential_RequiresResolution

    subroutine TIDEExponential_ResolveBackground(this, State, H0, omch2)
    class(TIDEExponential), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl), intent(out) :: H0, omch2
    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) :: y(NumEqs), yprime(NumEqs), w(NumEqs, 9)
    real(dl) :: loga, loga_end, a, grho_base, grho_v, grho_c, grho_plusp, grho_total
    integer :: i, ind

    H0 = 0._dl
    omch2 = 0._dl
    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2
    if (allocated(this%phi_table)) deallocate(this%phi_table, this%u_table, this%du_table)
    if (allocated(this%nu_grho)) deallocate(this%nu_grho, this%nu_mass)

    if (.not. ieee_is_finite(this%V0_over_rho100) .or. this%V0_over_rho100 <= 0._dl .or. &
        .not. ieee_is_finite(this%omc_early) .or. this%omc_early <= 0._dl .or. &
        .not. ieee_is_finite(this%beta) .or. .not. ieee_is_finite(this%lam) .or. &
        .not. ieee_is_finite(this%a_init) .or. this%a_init <= 0._dl .or. this%a_init >= 1._dl .or. &
        .not. ieee_is_finite(this%integrate_tol) .or. this%integrate_tol <= 0._dl .or. &
        this%history_points < 16) then
        call GlobalError('IDE exponential: invalid input parameters', error_darkenergy)
        return
    end if

    select type (State)
    class is (CAMBdata)
        this%grhog0 = State%grhog
        this%grhornomass0 = State%grhornomass
        this%grhob0 = State%grhob
        allocate(this%nu_grho(size(State%grhormass)), this%nu_mass(size(State%nu_masses)))
        this%nu_grho = State%grhormass
        this%nu_mass = State%nu_masses
    class default
        call GlobalError('IDE exponential: expected CAMBdata', error_darkenergy)
        return
    end select

    if (.not. ieee_is_finite(this%grhog0) .or. .not. ieee_is_finite(this%grhornomass0) .or. &
        .not. ieee_is_finite(this%grhob0) .or. any(.not. ieee_is_finite(this%nu_grho)) .or. &
        any(.not. ieee_is_finite(this%nu_mass)) .or. this%grhog0 <= 0._dl .or. this%grhob0 <= 0._dl .or. &
        this%grhornomass0 < 0._dl .or. any(this%nu_grho < 0._dl) .or. any(this%nu_mass < 0._dl)) then
        call GlobalError('IDE exponential: invalid uncoupled background density', error_darkenergy)
        return
    end if

    this%log_a_init = log(this%a_init)
    this%dlog_a = -this%log_a_init/(this%history_points - 1)
    allocate(this%phi_table(this%history_points), this%u_table(this%history_points), &
        this%du_table(this%history_points))

    ! Initial conditions are phi(a_init)=0 and dphi/dlog(a)=0.
    y = 0._dl
    loga = this%log_a_init
    rk_settings%reuse_first_derivative = .true.
    ind = 1
    do i = 1, this%history_points
        a = exp(loga)
        call this%EvolveLogA(NumEqs, loga, y, yprime)
        if (global_error_flag /= 0) return
        this%phi_table(i) = y(1)
        this%u_table(i) = y(2)
        this%du_table(i) = yprime(2)
        if (i == this%history_points) exit
        loga_end = this%log_a_init + this%dlog_a*i
        if (i == this%history_points - 1) loga_end = 0._dl
        call RungeKuttaDP45(this, NumEqs, TIDEExponential_EvolveLogA, loga, y, loga_end, &
            this%integrate_tol, ind, rk_settings, NumEqs, w)
        if (ind < 0 .or. global_error_flag /= 0) then
            if (global_error_flag == 0) &
                call GlobalError('IDE exponential: scalar integration failed', error_darkenergy)
            return
        end if
    end do

    call this%TermsAtA(1._dl, y(1), grho_base, grho_v, grho_c, grho_plusp)
    if (global_error_flag /= 0) return
    if (y(2)**2 >= 6._dl) then
        call GlobalError('IDE exponential: non-positive Friedmann denominator', error_darkenergy)
        return
    end if
    grho_total = grho_base/(1._dl - y(2)**2/6._dl)
    if (.not. ieee_is_finite(grho_total) .or. grho_total <= 0._dl) then
        call GlobalError('IDE exponential: invalid present-day expansion', error_darkenergy)
        return
    end if
    H0 = 100._dl*sqrt(grho_total/grho100)
    omch2 = grho_c/grho100

    end subroutine TIDEExponential_ResolveBackground

    subroutine TIDEExponential_Init(this, State)
    class(TIDEExponential), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: phi, u, grho_base, grho_v, grho_c, grho_plusp, grho_total, grho_phi

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2
    if (.not. allocated(this%phi_table)) &
        call GlobalError('IDE exponential: background was not resolved', error_darkenergy)
    if (global_error_flag /= 0) return

    select type (State)
    class is (CAMBdata)
        call this%FieldAtA(1._dl, phi, u)
        call this%TermsAtA(1._dl, phi, grho_base, grho_v, grho_c, grho_plusp)
        if (global_error_flag /= 0) return
        grho_total = grho_base/(1._dl - u**2/6._dl)
        grho_phi = grho_v + grho_total*u**2/6._dl
        if (abs(grho_phi - State%grhov) > 1e-8_dl*grho_total) &
            call GlobalError('IDE exponential: scalar density does not close at a=1', error_darkenergy)
    class default
        call GlobalError('IDE exponential: expected CAMBdata', error_darkenergy)
    end select

    end subroutine TIDEExponential_Init

    subroutine TIDEExponential_EvolveLogA(this, num, loga, y, yprime)
    class(TIDEExponential), intent(in) :: this
    integer, intent(in) :: num
    real(dl), intent(in) :: loga, y(num)
    real(dl), intent(out) :: yprime(num)
    real(dl) :: a, grho_base, grho_v, grho_c, grho_plusp, grho_total, epsilon

    yprime = 0._dl
    if (global_error_flag /= 0) return
    if (any(.not. ieee_is_finite(y)) .or. y(2)**2 >= 6._dl) then
        call GlobalError('IDE exponential: scalar trajectory left physical domain', error_darkenergy)
        return
    end if

    a = exp(loga)
    call this%TermsAtA(a, y(1), grho_base, grho_v, grho_c, grho_plusp)
    if (global_error_flag /= 0) return
    grho_total = grho_base/(1._dl - y(2)**2/6._dl)
    epsilon = 1.5_dl*grho_plusp/grho_total + 0.5_dl*y(2)**2
    yprime(1) = y(2)
    yprime(2) = -(3._dl - epsilon)*y(2) + 3._dl*(this%lam*grho_v + this%beta*grho_c)/grho_total

    if (any(.not. ieee_is_finite(yprime))) &
        call GlobalError('IDE exponential: non-finite scalar derivative', error_darkenergy)

    end subroutine TIDEExponential_EvolveLogA

    subroutine TIDEExponential_TermsAtA(this, a, phi, grho_base, grho_v, grho_c, grho_plusp)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a, phi
    real(dl), intent(out) :: grho_base, grho_v, grho_c, grho_plusp
    real(dl) :: grho_rad, grho_b, grho_nu, grho_pnu, rhonu, pnu
    integer :: i

    grho_rad = this%grhog0 + this%grhornomass0
    grho_b = this%grhob0*a
    grho_c = grho100*this%omc_early*a*exp(-this%beta*phi)
    grho_v = grho100*this%V0_over_rho100*a**4*exp(-this%lam*phi)
    grho_nu = 0._dl
    grho_pnu = 0._dl
    do i = 1, size(this%nu_grho)
        if (this%nu_grho(i) == 0._dl) cycle
        call ThermalNuBackground%rho_P(a*this%nu_mass(i), rhonu, pnu)
        grho_nu = grho_nu + this%nu_grho(i)*rhonu
        grho_pnu = grho_pnu + this%nu_grho(i)*pnu
    end do
    grho_base = grho_rad + grho_b + grho_c + grho_nu + grho_v
    grho_plusp = (4._dl/3._dl)*grho_rad + grho_b + grho_c + grho_nu + grho_pnu

    ! These a^4-weighted contributions may vanish exactly at a=0, but must be
    ! strictly positive on the forward integration interval.
    if (.not. ieee_is_finite(grho_base) .or. grho_base <= 0._dl .or. &
        .not. ieee_is_finite(grho_c) .or. grho_c < 0._dl .or. (a > 0._dl .and. grho_c == 0._dl) .or. &
        .not. ieee_is_finite(grho_v) .or. grho_v < 0._dl .or. (a > 0._dl .and. grho_v == 0._dl) .or. &
        .not. ieee_is_finite(grho_plusp)) &
        call GlobalError('IDE exponential: non-finite or non-positive density', error_darkenergy)

    end subroutine TIDEExponential_TermsAtA

    subroutine TIDEExponential_FieldAtA(this, a, phi, u)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: phi, u
    real(dl) :: x, s, s2, s3, h, h00, h10, h01, h11
    integer :: i

    phi = 0._dl
    u = 0._dl
    if (.not. allocated(this%phi_table)) then
        call GlobalError('IDE exponential: background history unavailable', error_darkenergy)
        return
    end if
    if (a <= this%a_init) return
    if (a >= 1._dl) then
        phi = this%phi_table(this%history_points)
        u = this%u_table(this%history_points)
        return
    end if

    x = (log(a) - this%log_a_init)/this%dlog_a
    i = min(max(int(x) + 1, 1), this%history_points - 1)
    s = x - (i - 1)
    s2 = s*s
    s3 = s2*s
    h = this%dlog_a
    h00 = 2._dl*s3 - 3._dl*s2 + 1._dl
    h10 = s3 - 2._dl*s2 + s
    h01 = -2._dl*s3 + 3._dl*s2
    h11 = s3 - s2
    phi = h00*this%phi_table(i) + h10*h*this%u_table(i) + &
        h01*this%phi_table(i + 1) + h11*h*this%u_table(i + 1)
    u = h00*this%u_table(i) + h10*h*this%du_table(i) + &
        h01*this%u_table(i + 1) + h11*h*this%du_table(i + 1)

    end subroutine TIDEExponential_FieldAtA

    real(dl) function TIDEExponential_PhiAtA(this, a)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: u

    call this%FieldAtA(a, TIDEExponential_PhiAtA, u)

    end function TIDEExponential_PhiAtA

    real(dl) function TIDEExponential_InitialScaleFactor(this)
    class(TIDEExponential), intent(in) :: this

    TIDEExponential_InitialScaleFactor = this%a_init

    end function TIDEExponential_InitialScaleFactor

    real(dl) function TIDEExponential_CouplingBeta(this)
    class(TIDEExponential), intent(in) :: this

    TIDEExponential_CouplingBeta = this%beta

    end function TIDEExponential_CouplingBeta

    subroutine TIDEExponential_ScalarBackgroundAtA(this, a, phi, phidot, adotoa)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: phi, phidot
    real(dl), intent(in), optional :: adotoa
    real(dl) :: u, grho_base, grho_v, grho_c, grho_plusp, grho_total

    phi = 0._dl
    phidot = 0._dl
    if (a <= 0._dl) return
    call this%FieldAtA(a, phi, u)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(phi) .or. .not. ieee_is_finite(u) .or. u**2 >= 6._dl) then
        call GlobalError('IDE exponential: invalid interpolated scalar trajectory', error_darkenergy)
        return
    end if
    if (u == 0._dl) return
    if (present(adotoa)) then
        if (.not. ieee_is_finite(adotoa) .or. adotoa <= 0._dl) then
            call GlobalError('IDE exponential: invalid conformal Hubble rate', error_darkenergy)
            return
        end if
        phidot = u*adotoa
        return
    end if
    call this%TermsAtA(a, phi, grho_base, grho_v, grho_c, grho_plusp)
    if (global_error_flag /= 0) return
    grho_total = grho_base/(1._dl - u**2/6._dl)
    phidot = u*sqrt(grho_total/3._dl)/a

    end subroutine TIDEExponential_ScalarBackgroundAtA

    real(dl) function TIDEExponential_CoupledFieldSource(this, a, grhoc_t, delta_c)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a, grhoc_t, delta_c

    ! grhoc_t is kappa*a**2*rho_c; delta_c includes the perturbation of the CDM mass.
    TIDEExponential_CoupledFieldSource = this%beta*grhoc_t*delta_c

    end function TIDEExponential_CoupledFieldSource

    real(dl) function TIDEExponential_ExpansionLogSlope(this, a)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: phi, u, grho_base, grho_v, grho_c, grho_plusp, grho_total

    TIDEExponential_ExpansionLogSlope = 0._dl
    if (a <= 0._dl) then
        call GlobalError('IDE exponential: scale factor must be positive', error_darkenergy)
        return
    end if
    call this%FieldAtA(a, phi, u)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(phi) .or. .not. ieee_is_finite(u) .or. u**2 >= 6._dl) then
        call GlobalError('IDE exponential: invalid interpolated scalar trajectory', error_darkenergy)
        return
    end if
    call this%TermsAtA(a, phi, grho_base, grho_v, grho_c, grho_plusp)
    if (global_error_flag /= 0) return
    grho_total = grho_base/(1._dl - u**2/6._dl)
    TIDEExponential_ExpansionLogSlope = 1.5_dl*grho_plusp/grho_total + 0.5_dl*u**2

    end function TIDEExponential_ExpansionLogSlope

    real(dl) function TIDEExponential_RhoCOverRho100(this, a)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: phi

    if (a <= 0._dl) then
        call GlobalError('IDE exponential: scale factor must be positive', error_darkenergy)
        TIDEExponential_RhoCOverRho100 = 0._dl
        return
    end if
    phi = this%PhiAtA(a)
    TIDEExponential_RhoCOverRho100 = this%omc_early*exp(-this%beta*phi)/a**3

    end function TIDEExponential_RhoCOverRho100

    real(dl) function TIDEExponential_CDMGrhoa4(this, a, grhoc0)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(in) :: a, grhoc0
    real(dl) :: phi

    phi = this%PhiAtA(a)
    TIDEExponential_CDMGrhoa4 = grho100*this%omc_early*a*exp(-this%beta*phi)

    end function TIDEExponential_CDMGrhoa4

    subroutine TIDEExponential_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    class(TIDEExponential), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), intent(out), optional :: w
    real(dl) :: phi, u, grho_base, grho_v, grho_c, grho_plusp, grho_total, grho_kin

    if (a <= sqrt(tiny(1._dl))) then
        grhov_t = 0._dl
        if (present(w)) w = -1._dl
        return
    end if
    call this%FieldAtA(a, phi, u)
    call this%TermsAtA(a, phi, grho_base, grho_v, grho_c, grho_plusp)
    if (global_error_flag /= 0 .or. u**2 >= 6._dl) then
        grhov_t = 0._dl
        if (present(w)) w = -1._dl
        return
    end if
    grho_total = grho_base/(1._dl - u**2/6._dl)
    grho_kin = grho_total*u**2/6._dl
    grhov_t = (grho_kin + grho_v)/a**2
    if (present(w)) then
        if (grho_kin + grho_v > 0._dl) then
            w = (grho_kin - grho_v)/(grho_kin + grho_v)
        else
            w = -1._dl
        end if
    end if

    end subroutine TIDEExponential_BackgroundDensityAndPressure

    subroutine TIDEExponential_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TIDEExponential), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: phi, phidot, potential

    call this%ScalarBackgroundAtA(a, phi, phidot, adotoa)
    if (global_error_flag /= 0) then
        dgrhoe = 0._dl
        dgqe = 0._dl
        return
    end if
    potential = grho100*this%V0_over_rho100*exp(-this%lam*phi)
    dgrhoe = phidot*ay(w_ix + 1) - this%lam*a**2*potential*ay(w_ix)
    dgqe = k*phidot*ay(w_ix)

    end subroutine TIDEExponential_PerturbedStressEnergy

    subroutine TIDEExponential_PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TIDEExponential), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:), w
    integer, intent(in) :: w_ix
    real(dl) :: phi, phidot, potential

    ! The two variables are delta_phi and its conformal-time derivative. The CDM
    ! density source is added by the scalar equations after this homogeneous part.
    call this%ScalarBackgroundAtA(a, phi, phidot, adotoa)
    if (global_error_flag /= 0) then
        ayprime(w_ix:w_ix + 1) = 0._dl
        return
    end if
    potential = grho100*this%V0_over_rho100*exp(-this%lam*phi)
    ayprime(w_ix) = y(w_ix + 1)
    ayprime(w_ix + 1) = -2._dl*adotoa*y(w_ix + 1) - k*z*phidot - &
        (k**2 + a**2*this%lam**2*potential)*y(w_ix)

    end subroutine TIDEExponential_PerturbationEvolve

    function TIDEExponential_PythonClass()
    character(len=:), allocatable :: TIDEExponential_PythonClass

    TIDEExponential_PythonClass = 'IDEExponential'

    end function TIDEExponential_PythonClass

    subroutine TIDEExponential_SelfPointer(cptr, P)
    use iso_c_binding, only: c_ptr, c_f_pointer
    type(c_ptr) :: cptr
    type(TIDEExponential), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TIDEExponential_SelfPointer

    end module DarkEnergyIDEExponential
