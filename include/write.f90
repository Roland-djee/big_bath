module write
  use types
  implicit none

contains

subroutine couplings_nuclear_outputs
  use constants, only : pi
  implicit none
  ! local variables
  integer(kind=8) :: k,l,m  
  character(len=100) :: filename

  character (len=*), parameter :: fmt05 = "(a1, t7, a, t22, a, t38, a, t58, a, &
       t79, a, t100, a, t121, a, t142, a)"
  character (len=*), parameter :: fmt4 = "(i5, i5, i5, i5, i5, i5, es20.10e3, es20.10e3, &
       es20.10e3, es20.10e3, es20.10e3, es20.10e3)"

  write(filename, '(a, es10.3e3, a)') "couplings_nuclear_J=",abs(Jnuc),"_rad.dat"
  open(14,file=filename)
    
  ! Headers
  write(14, fmt05)"#","spin 1","spin 2","distance","DC","C12","DJ","DC_n"
 
  m = 0
    do k=2,nb_imp - 1
       do l=k + 1,nb_imp
          m = m + 1
          write(14, fmt4)Imp_x(k), Imp_y(k), Imp_z(k),&
             Imp_x(l), Imp_y(l), Imp_z(l),&
             sqrt(dble(Imp_x(k) - Imp_x(l))**2&
             +dble(Imp_y(k) - Imp_y(l))**2&
             +dble(Imp_z(k) - Imp_z(l))**2),&
             DC(m), C12(m), DJ(m), DC_n(m)

       end do
    end do

end subroutine couplings_nuclear_outputs

subroutine couplings_outputs
  implicit none
  ! local variables
  integer(kind=8) :: i,j,k

  character (len=*), parameter :: fmt01 = "(a1, t5, a, t21, a, t41, a)"
  character (len=*), parameter :: fmt02 = "(a1, t10, a, t25, a, t40, a, t55, a)"

  character (len=*), parameter :: fmt1 = "(i5, i5, i5, es20.10e3, es20.10e3)"
  character (len=*), parameter :: fmt2 = "(i5, i5, i5, i5, i5, i5, es20.10e3, es20.10e3)"
  
  open(14,file='hyperfine.dat')
  open(15,file='dipolar.dat')
  open(16,file='dipolar_nuclear.dat')
    
  ! Headers
  write(14, fmt01)"#","spin 1","distance","J [rad/s]"
  write(15, fmt02)"#","spin 1","spin 2","distance","C12 [rad/s]"
  write(16, fmt01)"#","spin 1","distance","C12 nuclear [rad/s]"

  do i = 1, nb_imp
     write(14,fmt1)Imp_x(i), Imp_y(i), Imp_z(I),&
          sqrt(dble(Imp_x(i))**2&
          +dble(Imp_y(i))**2&
          +dble(Imp_z(I))**2),&
          HYPER(i)
     write(16,fmt1)Imp_x(i), Imp_y(i), Imp_z(I),&
          sqrt(dble(Imp_x(i))**2&
          +dble(Imp_y(i))**2&
          +dble(Imp_z(I))**2),&
          C12_n(i)
  end do
 
  k = 0
  do i = 1, nb_imp - 1 
     do j = i + 1, nb_imp
        k = k + 1 
        write(15,fmt2)Imp_x(i), Imp_y(i), Imp_z(I),&
             Imp_x(j), Imp_y(j), Imp_z(j),&
             sqrt(dble(Imp_x(i) - Imp_x(j))**2&
             +dble(Imp_y(i) - Imp_y(j))**2&
             +dble(Imp_z(I) - Imp_z(j))**2),&
             C12(k)       
     end do
  end do
 
end subroutine couplings_outputs

subroutine impurities_outputs
  implicit none
  integer(kind=8) :: i
  
  character (len=*), parameter :: fmt1 = "(i5, i5, i5)"
  
  open(12, file="Impurities_coord.dat")
  do i=1,nb_imp
     write(12, fmt1)Imp_x(i), Imp_y(i), Imp_z(i)
  end do
  close(12)
end subroutine impurities_outputs

subroutine crystal_outputs
  implicit none
  ! local variables
  integer(kind=8) :: i
  
  character (len=*), parameter :: fmt1 = "(i5, i5, i5)"
  
  open(11, file="Crystal_sites_coord.dat")
  do i=1,nb_sites
     write(11, fmt1)Sites_x(i), Sites_y(i), Sites_z(i)
  end do
  close(11)
end subroutine crystal_outputs

end module write
