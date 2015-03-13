module test_interactions
  use interactions
  implicit none

contains

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
