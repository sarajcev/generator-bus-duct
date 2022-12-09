! ...
! Ispis rezultata struja faznih vodica za graficki prikaz (Gnuplot)
subroutine ispis_struje(Nc,Ns,Nuk,dFis,Ii)
    implicit none
    integer,intent(in) :: Nc,Ns,Nuk
    real(8),intent(in) :: dFis
    complex(8),dimension(:),intent(in) :: Ii
    real(8),parameter :: PI = 3.141592653589793d0
    real(8) Fis
    real(8) Modul,Kut
    integer i

    ! OKLOP: FAZA R
    open(unit=3,file='Struja_oklop_R.txt',action='write')
    write(3,'(2x,"Fi",7x,"Iznos",5x,"Kut")')
    Fis = 0.d0
    do i = 3*Nc+1,3*Nc+Ns
        call zmk(Ii(i),Modul,Kut)
        write(3,'(f6.2,2x,f8.2,2x,f6.1)') Fis*(180.d0/PI),Modul,Kut
        Fis = Fis + dFis
    end do
    close(3)

    ! OKLOP: FAZA S
    open(unit=3,file='Struja_oklop_S.txt',action='write')
    write(3,'(2x,"Fi",7x,"Iznos",5x,"Kut")')
    Fis = 0.d0
    do i = 3*Nc+Ns+1,3*Nc+2*Ns
        call zmk(Ii(i),Modul,Kut)
        write(3,'(f6.2,2x,f8.2,2x,f6.1)') Fis*(180.d0/PI),Modul,Kut
        Fis = Fis + dFis
    end do
    close(3)

    ! OKLOP: FAZA T
    open(unit=3,file='Struja_oklop_T.txt',action='write')
    write(3,'(2x,"Fi",7x,"Iznos",5x,"Kut")')
    Fis = 0.d0
    do i = 3*Nc+2*Ns+1,Nuk
        call zmk(Ii(i),Modul,Kut)
        write(3,'(f6.2,2x,f8.2,2x,f6.1)') Fis*(180.d0/PI),Modul,Kut
        Fis = Fis + dFis
    end do
    close(3)

end subroutine