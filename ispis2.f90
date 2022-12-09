! ...
! Subroutine za ispis dijela rezultata proracuna (u file i na ekran)
subroutine ispis2(I1,I2,I3,Ir,Is,It,Iro,Iso,Ito,Iz,DV,DR,DS,DT,Pr,Ps,Pt,Puc,Por,Pos,Pot,Pus,Puk)
	implicit none

	complex(8),intent(in) :: Ir,Is,It,Iro,Iso,Ito,Iz,I1,I2,I3,DV,DR,DS,DT
	real(8),intent(in) :: Pr,Ps,Pt,Puc,Por,Pos,Pot,Pus,Puk
	real(8) Modul,Kut

	! -------------------------------------
	! Ispis vrijednosti izracunatih struja
	! -------------------------------------
	write(2,'(/,"Struje u slojevima faznog vodica:")')
	write(2,'("---------------------------------------")')
	call zmk(I1,Modul,Kut)
	write(2,'("1. sloj: I1 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(I2,Modul,Kut)
	write(2,'("2. sloj: I2 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(I3,Modul,Kut)
	write(2,'("3. sloj: I3 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Ekvivalentne fazne struje:")')
	write(2,'("---------------------------------------")')
	call zmk(Ir,Modul,Kut)
	write(2,'("Faza R: Ir =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Is,Modul,Kut)
	write(2,'("Faza S: Is =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(It,Modul,Kut)
	write(2,'("Faza T: It =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Ekvivalentne struje oklopa:")')
	write(2,'("---------------------------------------")')
	call zmk(Iro,Modul,Kut)
	write(2,'("Faza R: Iro =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Iso,Modul,Kut)
	write(2,'("Faza S: Iso =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Ito,Modul,Kut)
	write(2,'("Faza T: Ito =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Ekvivalentna struja u zemlji:")')
	write(2,'("---------------------------------------")')
	call zmk(Iz,Modul,Kut)
	write(2,'("Zemlja: Iz =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Razlika napona zvjezdista (f.v.):")')
	write(2,'("---------------------------------------")')
	call zmk(DV,Modul,Kut)
	write(2,'("DV =",f10.3," (V) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Padovi napona oklopa (EMS):")')
	write(2,'("---------------------------------------")')
	call zmk(DR,Modul,Kut)
	write(2,'("Faza R: DR =",f10.3," (V) /",f6.1," (o)")') Modul,Kut
	call zmk(DS,Modul,Kut)
	write(2,'("Faza S: DS =",f10.3," (V) /",f6.1," (o)")') Modul,Kut
	call zmk(DT,Modul,Kut)
	write(2,'("Faza T: DT =",f10.3," (V) /",f6.1," (o)")') Modul,Kut

	write(2,'(/,"Gubici snage oklopljenog voda:")')
	write(2,'("---------------------------------------")')
	write(2,'("Vodic - Faza R:",f10.2," (W/m)")') Pr
	write(2,'("Vodic - Faza S:",f10.2," (W/m)")') Ps
	write(2,'("Vodic - Faza T:",f10.2," (W/m)")') Pt
	write(2,'("---------------------------------------")')
	write(2,'("Ukupno (f. v.):",f10.2," (W/m)")') Puc
	write(2,'("---------------------------------------")')
	write(2,'("Oklop - Faza R:",f10.2," (W/m)")') Por
	write(2,'("Oklop - Faza S:",f10.2," (W/m)")') Pos
	write(2,'("Oklop - Faza T:",f10.2," (W/m)")') Pot
	write(2,'("---------------------------------------")')
	write(2,'("Ukupno (oklop):",f10.2," (W/m)")') Pus
	write(2,'("---------------------------------------")')
	write(2,'("Ukupni  gubici:",f10.2," (W/m)")') Puk
	close(2)

	write(*,'(/,"Struje u slojevima faznog vodica:")')
	write(*,'("---------------------------------------")')
	call zmk(I1,Modul,Kut)
	write(*,'("1. sloj: I1 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(I2,Modul,Kut)
	write(*,'("2. sloj: I2 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(I3,Modul,Kut)
	write(*,'("3. sloj: I3 =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Ekvivalentne fazne struje:")')
	write(*,'("---------------------------------------")')
	call zmk(Ir,Modul,Kut)
	write(*,'("Faza R: Ir =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Is,Modul,Kut)
	write(*,'("Faza S: Is =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(It,Modul,Kut)
	write(*,'("Faza T: It =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Ekvivalentne struje oklopa:")')
	write(*,'("---------------------------------------")')
	call zmk(Iro,Modul,Kut)
	write(*,'("Faza R: Iro =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Iso,Modul,Kut)
	write(*,'("Faza S: Iso =",f10.3," (A) /",f6.1," (o)")') Modul,Kut
	call zmk(Ito,Modul,Kut)
	write(*,'("Faza T: Ito =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Ekvivalentna struja u zemlji:")')
	write(*,'("---------------------------------------")')
	call zmk(Iz,Modul,Kut)
	write(*,'("Zemlja: Iz =",f10.3," (A) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Razlika napona zvjezdista (f.v.):")')
	write(*,'("---------------------------------------")')
	call zmk(DV,Modul,Kut)
	write(*,'("DV =",f10.3," (V) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Padovi napona oklopa (EMS):")')
	write(*,'("---------------------------------------")')
	call zmk(DR,Modul,Kut)
	write(*,'("Faza R: DR =",f10.3," (V) /",f6.1," (o)")') Modul,Kut
	call zmk(DS,Modul,Kut)
	write(*,'("Faza S: DS =",f10.3," (V) /",f6.1," (o)")') Modul,Kut
	call zmk(DT,Modul,Kut)
	write(*,'("Faza T: DT =",f10.3," (V) /",f6.1," (o)")') Modul,Kut

	write(*,'(/,"Gubici snage oklopljenog voda:")')
	write(*,'("---------------------------------------")')
	write(*,'("Vodic - Faza R:",f10.2," (W/m)")') Pr
	write(*,'("Vodic - Faza S:",f10.2," (W/m)")') Ps
	write(*,'("Vodic - Faza T:",f10.2," (W/m)")') Pt
	write(*,'("---------------------------------------")')
	write(*,'("Ukupno (f. v.):",f10.2," (W/m)")') Puc
	write(*,'("---------------------------------------")')
	write(*,'("Oklop - Faza R:",f10.2," (W/m)")') Por
	write(*,'("Oklop - Faza S:",f10.2," (W/m)")') Pos
	write(*,'("Oklop - Faza T:",f10.2," (W/m)")') Pot
	write(*,'("---------------------------------------")')
	write(*,'("Ukupno (oklop):",f10.2," (W/m)")') Pus
	write(*,'("---------------------------------------")')
	write(*,'("Ukupni  gubici:",f10.2," (W/m)")') Puk

end subroutine