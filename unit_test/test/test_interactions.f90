module test_interactions
  use interactions
  implicit none

contains

  subroutine unit_test_dipolar
    implicit none
    integer :: i,j,k    
    double precision :: Orient_norm
    integer(kind=8), allocatable :: X_dist(:)
    integer(kind=8), allocatable :: Y_dist(:)
    integer(kind=8), allocatable :: Z_dist(:)

    double precision, allocatable :: dist(:)
    double precision, allocatable :: proj(:)
    double precision, allocatable :: cos_t12(:)
    double precision, allocatable :: C12_2(:)

    logical :: test

    nb_imp = 20
    nb_pairs = nb_imp * (nb_imp - 1) / 2
    
    allocate(Imp_x(nb_imp),Imp_y(nb_imp),Imp_z(nb_imp))
    allocate(X_dist(nb_pairs),Y_dist(nb_pairs),Z_dist(nb_pairs))
    allocate(dist(nb_pairs),proj(nb_pairs),cos_t12(nb_pairs),C12_2(nb_pairs))

    Imp_x = (/(i,i=1,nb_imp)/)
    Imp_y = Imp_x
    Imp_z = Imp_x

    Orientation%x = 1
    Orientation%x = 2
    Orientation%x = 3

    call dipolar

    k = 0
    do i=1,nb_imp - 1
       do j=i + 1,nb_imp
          k = k + 1
          X_dist(k) = Imp_x(j) - Imp_x(i)
          Y_dist(k) = Imp_y(j) - Imp_y(i)
          Z_dist(k) = Imp_z(j) - Imp_z(i)
       end do
    end do

    Orient_norm = dsqrt(dble(Orientation%x**2) + &
                        dble(Orientation%y**2) + &
                        dble(Orientation%z**2))

    dist    = dsqrt(dble(X_dist**2) + dble(Y_dist**2) + dble(Z_dist**2))
    proj    = Orientation%x*X_dist + Orientation%y*Y_dist + Orientation%z*Z_dist
    cos_t12 = proj / (Orient_norm * dist)
    C12_2   = pref * (1.d0 - 3.d0 * cos_t12**2) / ((scaling * dist)**3)

    test = .true.
    do i=1,nb_pairs
       if(abs(C12(i) - C12_2(i)) .gt. 1.d-15) test = .false.
    end do

    if (test) then
       write(*,*)'Test dipolar successful'
    else
       write(*,*)'Test dipolar failed'
    end if

    open(10,file='dipolar_benchmark.dat')
    do i=1,nb_pairs
       write(10,*)scaling * dist(i), C12(i)
    end do

    write(*,*)'Dipolar_benchmark written'

  end subroutine unit_test_dipolar

  subroutine unit_test_hyperfine
    implicit none
    integer :: i
    double precision, allocatable :: X_imp(:), Y_imp(:), Z_imp(:)
    double precision, allocatable :: term1(:), term2(:), term3(:)
    double precision, allocatable :: HYPER2(:)
    logical :: test

    nb_imp = 10
    
    allocate(Imp_x(nb_imp),Imp_y(nb_imp),Imp_z(nb_imp))
    allocate(X_imp(nb_imp),Y_imp(nb_imp),Z_imp(nb_imp))
    allocate(term1(nb_imp),term2(nb_imp),term3(nb_imp))
    allocate(HYPER2(nb_imp))

    Imp_x = (/(i,i=1,nb_imp)/)
    Imp_y = Imp_x
    Imp_z = Imp_x

    call hyperfine

    X_imp = scaling * Imp_x
    Y_imp = scaling * Imp_y
    Z_imp = scaling * Imp_z

    term1 = F(X_imp, Y_imp, Z_imp) * dcos(k0 * X_imp)
    term2 = F(Y_imp, Z_imp, X_imp) * dcos(k0 * Y_imp)
    term3 = F(Z_imp, X_imp, Y_imp) * dcos(k0 * Z_imp)

    HYPER2 = p*(term1+term2+term3)**2

    test = .true.
    do i=1,nb_imp
       if(abs(HYPER(i) - HYPER2(i)) .gt. 1.d-15) test = .false.
    end do

    if (test) then
       write(*,*)'Test hyperfine successful'
    else
       write(*,*)'Test hyperfine failed'
    end if

    deallocate (HYPER)
    deallocate(Imp_x,Imp_y,Imp_z)
    nb_imp = 1000  
    allocate(Imp_x(nb_imp),Imp_y(nb_imp),Imp_z(nb_imp))

    Imp_x = (/(i,i=1,nb_imp)/)
    Imp_y = 0
    Imp_z = 0

    call hyperfine

    open(10,file='hyperfine_benchmark.dat')
    do i=1,nb_imp
       write(10,*)scaling * Imp_x(i),HYPER(i)
    end do

    write(*,*)'Hyperfine_benchmark written'

    deallocate(Imp_x,Imp_y,Imp_z)
    deallocate(X_imp,Y_imp,Z_imp)
    deallocate(term1,term2,term3)
    deallocate(HYPER2)

  end subroutine unit_test_hyperfine

  subroutine unit_test_F
    implicit none
    integer :: i,j,k,nb_imp
    double precision, allocatable :: x(:),y(:),z(:)
    double precision, allocatable :: test_Func1(:)
    double precision, allocatable :: test_Func2(:)
    double precision :: a,b,n

    logical :: test

    a = 25.09d0
    b = 14.43d0
    n = dsqrt(0.029d0/0.044d0)

    nb_imp = 10

    allocate(x(nb_imp),y(nb_imp),z(nb_imp))
    allocate(test_Func1(nb_imp),test_Func2(nb_imp))

    x = (/(dble(i),i=1,nb_imp)/)
    y = x
    z = x

    do i=1,nb_imp
       test_Func1(i) = -dsqrt(x(i)**2/((n*b)**2)+(y(i)**2+z(i)**2)/((n*a)**2))
       test_Func1(i) = dexp(test_Func1(i))/dsqrt(pi*(n*a)**2*(n*b))
    end do

    test = .true.
    test_Func2 = F(x,y,z)
    do i=1,nb_imp
       if(abs(test_Func1(i) - test_Func2(i)) .gt. 1.d-15) test = .false.
    end do

    if (test) then
       write(*,*)'Test F successful'
    else
       write(*,*)'Test F failed'
    end if

    deallocate(x,y,z)
    deallocate(test_Func1,test_Func2)

  end subroutine unit_test_F

end module test_interactions
