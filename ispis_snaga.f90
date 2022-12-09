! ...
! Ispis gubitaka snage oklopa za graficki prikaz (Gnuplot)
subroutine ispis_snaga(Ns,dFis,Pori,Posi,Poti)
    implicit none
    integer,intent(in) :: Ns
    real(8),intent(in) :: dFis
    real(8),dimension(:),intent(in) :: Pori,Posi,Poti
    real(8),parameter :: PI = 3.141592653589793d0
    real(8) Fis
    integer i

    ! OKLOP: FAZA R
    open(unit=3,file='Snaga_oklop_R.txt',action='write')
    write(3,'(2x,"Fi(o)",7x,"Pr(W/m)")')
    Fis = 0.d0
    do i = 1,Ns
        write(3,'(f6.2,2x,f8.4)') Fis*(180.d0/PI),Pori(i)
        Fis = Fis + dFis
    end do
    close(3)

    ! OKLOP: FAZA S
    open(unit=3,file='Snaga_oklop_S.txt',action='write')
    write(3,'(2x,"Fi(o)",7x,"Ps(W/m)")')
    Fis = 0.d0
    do i = 1,Ns
        write(3,'(f6.2,2x,f8.4)') Fis*(180.d0/PI),Posi(i)
        Fis = Fis + dFis
    end do
    close(3)

    ! OKLOP: FAZA T
    open(unit=3,file='Snaga_oklop_T.txt',action='write')
    write(3,'(2x,"Fi(o)",7x,"Pt(W/m)")')
    Fis = 0.d0
    do i = 1,Ns
        write(3,'(f6.2,2x,f8.4)') Fis*(180.d0/PI),Poti(i)
        Fis = Fis + dFis
    end do
    close(3)

end subroutine