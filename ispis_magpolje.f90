! ...
! Ispis rezultata proracuna magnetskog polja 
! po tockama kruznice za graficki prikaz (Gnuplot)
subroutine ispis_magpolje(Nq,theta,Htan,Hrad,H)
	implicit none
	integer,intent(in) :: Nq
	real(8),dimension(:),intent(in) :: theta
	complex(8),dimension(:),intent(in) :: Htan,Hrad,H
	real(8),parameter :: PI = 3.141592653589793d0
	real(8) Modul,Kut
	integer k

	open(unit=3,file='Magpolje_tan.txt',action='write')
	write(3,'(2x,"Fi",7x,"Ht(A/m)",5x,"(o)")')
	do k = 1,Nq
		call zmk(Htan(k),Modul,Kut)
		write(3,'(f6.2,2x,f10.3,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
	end do
	close(3)
	open(unit=3,file='Magpolje_rad.txt',action='write')
	write(3,'(2x,"Fi",7x,"Hr(A/m)",5x,"(o)")')
	do k = 1,Nq
		call zmk(Hrad(k),Modul,Kut)
		write(3,'(f6.2,2x,f10.3,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
	end do
	close(3)
	open(unit=3,file='Magpolje.txt',action='write')
	write(3,'(2x,"Fi",7x,"H(A/m)",5x,"(o)")')
	do k = 1,Nq
		call zmk(H(k),Modul,Kut)
		write(3,'(f6.2,2x,f10.3,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
	end do
	close(3)

end subroutine