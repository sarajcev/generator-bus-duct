! ########################################################################
! #                                                                      #
! #    RASPODJELA STRUJE, GUBITAKA SNAGE U OKLOPLJENIM GENERATORSKIM     #
! #    VODOVIMA OKRUGLOG POPRECNOG PRESJEKA, TE PRORACUN RASPODJELE      #
! #        ELEKTROMAGNETSKOG POLJA U OKOLISU OKLOPLJENOG VODA            #
! #                                                                      #
! ########################################################################
! Matematicki model temelji se na metodi dionih vodica i SGU iz doktorske
! disertacije Prof. I. Sarajceva. Oklopljeni vod se sastoji od jednofazno
! oklopljenih supljih okruglih vodica u okruglom (cilindricnom) oklopu.
! ========================================================================
! Autor: Dr.sc. Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)

PROGRAM gbdem_main
implicit none

! Ulazne varijable:
integer N1,N2,N3,Ns
real(8) hc,Rc,Kc
real(8) hs,Rs,Ks
real(8) D0
real(8) L,ro,f
complex(8) Zrt,Zst,Ztt
complex(8) Zgu,Ztu
real(8) Vr,Fr,Vs,Fs,Vt,Ft
integer Nq
real(8) xs,ys,Rq
integer Nqk
real(8) xstart,ystart,dx

! Definicija ulaznih varijabli:
!   N1  - broj dionih vodica 1. sloja faznog vodica
!   N2  - broj dionih vodica 2. sloja faznog vodica
!   N3  - broj dionih vodica 3. sloja faznog vodica
!   Ns  - ukupni broj dionih vodica jednog oklopa
!   hc  - debljina faznog vodica u oklopu (mm)
!   Rc  - unutarnji radijus faznog vodica u oklopu (mm)
!   Kc  - specificna vodljivost materijala faznog vodica (S*m/mm2)
!   hs  - debljina oklopa (mm)
!   Rs  - srednji radijus oklopa (mm)
!   Ks  - specificna vodljivost materijala oklopa (S*m/mm2)
!   D0  - razmak medju sredistima faza (mm)
!   L   - ukupna duljina oklopljenih vodova (m)
!   ro  - prosjecna specificna elektricna otpornost tla (ohm*m)
!   f   - frekvencija (Hz)
!   Zrt - impedancija na teretnoj strani oklopljenog voda u fazi R (ohm)
!   Zst - impedancija na teretnoj strani oklopljenog voda u fazi S (ohm)
!   Ztt - impedancija na teretnoj strani oklopljenog voda u fazi T (ohm)
!   Zgu - impedancija uzemljenja (uzemljivaca) na generatorskoj strani
!         oklopljenog voda (ohm)
!   Ztu - impedancija uzemljenja (uzemljivaca) na teretnoj strani
!         oklopljenog voda (ohm)
!   Vr  - iznos faznog napona na generatorskoj strani u fazi R (V)
!   Fr  - kut faznog napona na generatorskoj strani u fazi R (o)
!   Vs  - iznos faznog napona na generatorskoj strani u fazi S (V)
!   Fs  - kut faznog napona na generatorskoj strani u fazi R (o)
!   Vt  - iznos faznog napona na generatorskoj strani u fazi T (V)
!   Ft  - kut faznog napona na generatorskoj strani u fazi R (o)

!   Nq  - broj tocaka na kruznici na kojoj se racuna raspodjela magnetskog polja
!   xs  - x-koordinata sredista kruznice (mm)
!   ys  - y-koordinata sredista kruznice (mm)
!   Rq  - radijus kruznice (mm)

!   Nqk - broj tocaka duz pravca za proracun magnetskog polja
!   xstart - x-koordinata pocetne tocke pravca (mm)
!   ystart - y-koordinata pocetne tocke pravca (mm)
!   dx  - korak proracuna duz promatranog pravca (mm)

! ------------------------------
! Ulazna datoteka: inputdata.txt
! ------------------------------
! N1 N2 N3 Ns
! hc Rc Kc
! hs Rs Ks
! D0 L ro f
! Re(Zrt) Im(Zrt)
! Re(Zst) Im(Zst)
! Re(Ztt) Im(Ztt)
! Re(Zgu) Im(Zgu)
! Re(Ztu) Im(Ztu)
! Vr Fr
! Vs Fs
! Vt Ft
! Nq xs ys Rq
! Nk xstart ystart dx
! ------------------------------

! Napomena: Jedinstveni koordinatni sustav oklopljenih generator. vodova,
! na temelju kojeg se zadaje sva geometrija sustava, polozen je u centru
! (sredistu) lijevog faznog vodica. Koordinatna os x polozena je horizo-
! ntalno i prolazi sredistima svih faznih vodica, dok je koordinatna os
! y polozena okomito i prolazi sredistem lijevog faznog vodica. Naime, 
! faze oklopljenog voda polozene su horizontalno u ravnini (duz osi x).

! Napomena: Pravac za prorasun raspodjele magnetskog polja paralelan je
! s koordinatnom osi x (globalni koord. sustav)!

! Ostale varijable:
integer Nc,Nuk
real(8) dFic1,sc1,dFic2,sc2,dFic3,sc3,Fic
real(8) dFis,ss,Fis
real(8),dimension(:),allocatable :: x,y
real(8),dimension(:,:),allocatable :: d
real(8),dimension(:),allocatable :: R1
complex(8),dimension(:,:),allocatable :: Z
complex(8) Zii,Zik
real(8) De,omega,mio
complex(8) V1,V2,V3
complex(8),dimension(:),allocatable :: V,Ii
complex(8) Ir,Is,It,Iro,Iso,Ito,Iz,I1,I2,I3
complex(8) DV,DR,DS,DT,DV2
real(8) Pr,Ps,Pt,Por,Pos,Pot,Puc,Pus,Puk
real(8),dimension(:),allocatable :: Pori,Posi,Poti
real(8) Re,Im,Modul,Kut
real(8) ac,bc1,bc2,bc3,as,bs
real(8) K1,K2
real(8) Rphr,Rphs,Rpht,Rokr,Roks,Rokt,Reqr,Reqs,Reqt
real(8) dTHq,THq
real(8),dimension(:),allocatable :: xq,yq
real(8),dimension(:),allocatable :: delta,theta
real(8) Hk,He,Hkk,Hek
real(8),dimension(:,:),allocatable :: Hksi,Heta
complex(8),dimension(:),allocatable :: Htan,Hrad,H
complex(8) csum,Ampere,Isum
real(8),dimension(:),allocatable :: Fm
real(8) tstart1,tend1,tstart2,tend2,time1,time2,trun
real(8),dimension(:),allocatable :: xk,yk
real(8),dimension(:),allocatable :: thetak
real(8),dimension(:,:),allocatable :: Hksi_pr,Heta_pr
complex(8),dimension(:),allocatable :: Htan_pr,Hrad_pr,H_pr
complex(8),dimension(:),allocatable :: Bfc,Bfp
integer izbor
integer i,k

! Definicija ostalih vaznijih varijabli:
!   Nc  - ukupni broj dionih vodica jednog faznog vodica
!   Nuk - sveukupni broj dionih vodica oklopljenog generatorskog sustava
!   x,y - vektori koordinata dionih vodica oklopljenog sustava
!   d   - matrica vlastitih i medjusobnih udaljenosti dionih vodica
!   R1  - jedinicni otpori dionih vodica
!   Z   - matrica vlastitih i medjusobnih impedancija sustava
!   De  - dubina prodiranja povratnih struja u tlu
!   Ir  - izracunata ekvivalentna struja faze R
!   Is  - izracunata ekvivalentna struja faze S
!   It  - izracunata ekvivalentna struja faze T
!   Iro - izracunata ekvivalentna struja oklopa faze R
!   Iso - izracunata ekvivalentna struja oklopa faze S
!   Ito - izracunata ekvivalentna struja oklopa faze T
!   Iz  - izracunata ekvivalentna struja u zemlji
!   DV  - razlika napona izmedju dvaju krajeva oklopljenog voda

! Lapack varijable zsysv:
character(1) UPLO
integer N,NRHS,LDA,LDB,LWORK,INFO
integer,dimension(:),allocatable :: IPIV
complex(8),dimension(:),allocatable :: WORK
complex(8),dimension(:,:),allocatable :: B

! Konstante programa:
real(8),parameter :: PI = 3.141592653589793d0
complex(8),parameter :: one = dcmplx(1.d0,0.d0)
complex(8),parameter :: zero = dcmplx(0.d0,0.d0)

interface
    subroutine ispis_struje(Nc,Ns,Nuk,dFis,Ii)
        implicit none
        integer,intent(in) :: Nc,Ns,Nuk
        real(8),intent(in) :: dFis
        complex(8),dimension(:),intent(in) :: Ii
    end subroutine
    subroutine ispis_snaga(Ns,dFis,Pori,Posi,Poti)
        implicit none
        integer,intent(in) :: Ns
        real(8),intent(in) :: dFis
        real(8),dimension(:),intent(in) :: Pori,Posi,Poti
    end subroutine
    subroutine ispis_magpolje(Nq,theta,Htan,Hrad,H)
        implicit none
        integer,intent(in) :: Nq
        real(8),dimension(:),intent(in) :: theta
        complex(8),dimension(:),intent(in) :: Htan,Hrad,H
    end subroutine
end interface


! ======================================================================
!                          POCETAK PRORACUNA
! ======================================================================
call cpu_time(tstart1)
! Citanje ulaznih podataka iz vanjskog filea (inputdata.txt)
open(unit=1,file='inputdata.txt',action='read')
read(1,*) N1,N2,N3,Ns
Nc = N1 + N2 + N3
read(1,*) hc,Rc,Kc
read(1,*) hs,Rs,Ks
read(1,*) D0,L,ro,f
read(1,*) Re,Im
Zrt = dcmplx(Re,Im)
read(1,*) Re,Im
Zst = dcmplx(Re,Im)
read(1,*) Re,Im
Ztt = dcmplx(Re,Im)
read(1,*) Re,Im
Zgu = dcmplx(Re,Im)
read(1,*) Re,Im
Ztu = dcmplx(Re,Im)
read(1,*) Vr,Fr
read(1,*) Vs,Fs
read(1,*) Vt,Ft
read(1,*) Nq,xs,ys,Rq
read(1,*) Nqk,xstart,ystart,dx
close(1)

! ----------------------------------------------------------------------
! Ispis ulaznih podataka u vanjski file (output.txt)
! ----------------------------------------------------------------------
open(unit=2,file='output.txt',action='write')
write(2,'("########################################################################")')
write(2,'("#                                                                      #")')
write(2,'("#    RASPODJELA STRUJE I GUBITAKA SNAGE U OKLOPLJENIM GENERATORSKIM    #")')
write(2,'("#                VODOVIMA OKRUGLOG POPRECNOG PRESJEKA                  #")')
write(2,'("#                                                                      #")')
write(2,'("########################################################################")')
write(2,'("Proracun raspodjele struja  i gubitaka snage u oklopljenim cilindricnim")')
write(2,'("generatorskim vodovima uzemljenima na samo jednom ili pak na oba kraja.")')
write(2,'("------------------------------------------------------------------------")')
write(2,'("Autor: Dr.sc. Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)")')
write(2,'(//,"ULAZNI PODACI:",/)')
write(2,'("Fazni vodic se dijeli na tri sloja:")')
write(2,'("Broj dionih vodica za 1. sloj.......................:",i4)') N1
write(2,'("Broj dionih vodica za 2. sloj.......................:",i4)') N2
write(2,'("Broj dionih vodica za 3. sloj.......................:",i4)') N3
write(2,'("Ukupni broj podjela jednog faznog vodica............:",i4)') Nc
write(2,'("Debljina stjenke faznog vodica......................:",f6.1," (mm)")') hc
write(2,'("Unutarnji radijus faznog vodica.....................:",f6.1," (mm)")') Rc
write(2,'("Specificna vodljivost materijala faznog vodica......:",f8.2," (S*m/mm2)")') Kc
write(2,'(/,"Broj podjela jednog oklopa na dione vodice..........:",i4)') Ns
write(2,'("Debljina stjenke oklopa.............................:",f6.1," (mm)")') hs
write(2,'("Srednji radijus oklopa..............................:",f6.1," (mm)")') Rs
write(2,'("Specificna vodljivost materijala oklopa.............:",f8.2," (S*m/mm2)")') Ks
write(2,'(/,"Razmak izmedju oklopljenih faza.....................:",f8.1," (mm)")') D0
write(2,'("Duljina oklopljenih vodova..........................:",f8.1," (m)")') L
write(2,'("Specificna elektricna otpornost tla.................:",f6.1," (Ohm*m)")') ro
write(2,'("Nazivna frekvencija.................................:",f6.1," (Hz)")') f
write(2,'(/,"Impedancija tereta u fazi R.........................:",f7.3," j",f7.3," (ohm)")') Zrt
write(2,'("Impedancija tereta u fazi S.........................:",f7.3," j",f7.3," (ohm)")') Zst
write(2,'("Impedancija tereta u fazi T.........................:",f7.3," j",f7.3," (ohm)")') Ztt
write(2,'(/,"Impedancija uzemljivaca na generatorskoj strani.....:",f7.3," j",f7.3," (ohm)")') Zgu
write(2,'("Impedancija uzemljivaca na teretnoj strani..........:",f7.3," j",f7.3," (ohm)")') Ztu
write(2,'(/,"Fazni napon faze R..................................:",f8.1," (V)",f6.1," (o)")') Vr,Fr
write(2,'("Fazni napon faze S..................................:",f8.1," (V)",f6.1," (o)")') Vs,Fs
write(2,'("Fazni napon faze T..................................:",f8.1," (V)",f6.1," (o)")') Vt,Ft
write(2,'(/,"Kruznica za proracun raspodjele magnetskog polja:")')
write(2,'("Broj tocaka na kruznici za proracun polja...........:",i5)') Nq
write(2,'("Koordinate sredista kruznice (xs,ys)................:",f8.1," ,",f8.1," (mm)")') xs,ys
write(2,'("Radijus kruznice....................................:",f8.1," (mm)")') Rq
write(2,'(/,"Pravac za proracun raspodjele magnetskog polja:")')
write(2,'("Broj tocaka duz pravca za proracun polja...........:",i5)') Nqk
write(2,'("Koordinate pocetne tocke pravca (xs,ys)............:",f8.1," ,",f8.1," (mm)")') xstart,ystart
write(2,'("Korak proracuna (pravac je paralelan s osi x)......:",f8.1," (mm)")') dx
close(2)
! ----------------------------------------------------------------------
! Ispis ulaznih podataka na ekranu (echo)
! ----------------------------------------------------------------------
write(*,'("########################################################################")')
write(*,'("#                                                                      #")')
write(*,'("#    RASPODJELA STRUJE I GUBITAKA SNAGE U OKLOPLJENIM GENERATORSKIM    #")')
write(*,'("#                VODOVIMA OKRUGLOG POPRECNOG PRESJEKA                  #")')
write(*,'("#                                                                      #")')
write(*,'("########################################################################")')
write(*,'("Proracun raspodjele struja  i gubitaka snage u oklopljenim cilindricnim")')
write(*,'("generatorskim vodovima uzemljenima na samo jednom ili pak na oba kraja.")')
write(*,'("------------------------------------------------------------------------")')
write(*,'("Autor: Dr.sc. Petar Sarajcev, dipl.ing. (petar.sarajcev@fesb.hr)")')
write(*,'(//,"ULAZNI PODACI:",/)')
write(*,'("Fazni vodic se dijeli na tri sloja:")')
write(*,'("Broj dionih vodica za 1. sloj.......................:",i4)') N1
write(*,'("Broj dionih vodica za 2. sloj.......................:",i4)') N2
write(*,'("Broj dionih vodica za 3. sloj.......................:",i4)') N3
write(*,'("Broj podjela jednog faznog vodica na dione vodice...:",i4)') Nc
write(*,'("Debljina stjenke faznog vodica......................:",f6.1," (mm)")') hc
write(*,'("Unutarnji radijus faznog vodica.....................:",f6.1," (mm)")') Rc
write(*,'("Specificna vodljivost materijala faznog vodica......:",f8.2," (S*m/mm2)")') Kc
write(*,'(/,"Broj podjela jednog oklopa na dione vodice..........:",i4)') Ns
write(*,'("Debljina stjenke oklopa.............................:",f6.1," (mm)")') hs
write(*,'("Srednji radijus oklopa..............................:",f6.1," (mm)")') Rs
write(*,'("Specificna vodljivost materijala oklopa.............:",f8.2," (S*m/mm2)")') Ks
write(*,'(/,"Razmak izmedju oklopljenih faza.....................:",f8.1," (mm)")') D0
write(*,'("Duljina oklopljenih vodova..........................:",f8.1," (m)")') L
write(*,'("Specificna elektricna otpornost tla.................:",f6.1," (Ohm*m)")') ro
write(*,'("Nazivna frekvencija.................................:",f6.1," (Hz)")') f
write(*,'(/,"Impedancija tereta u fazi R.........................:",f6.3," j",f6.3," (ohm)")') Zrt
write(*,'("Impedancija tereta u fazi S.........................:",f6.3," j",f6.3," (ohm)")') Zst
write(*,'("Impedancija tereta u fazi T.........................:",f6.3," j",f6.3," (ohm)")') Ztt
write(*,'(/,"Impedancija uzemljivaca na generatorskoj strani.....:",f6.3," j",f6.3," (ohm)")') Zgu
write(*,'("Impedancija uzemljivaca na teretnoj strani..........:",f6.3," j",f6.3," (ohm)")') Ztu
write(*,'(/,"Fazni napon faze R..................................:",f8.1," (V)",f6.1," (o)")') Vr,Fr
write(*,'("Fazni napon faze S..................................:",f8.1," (V)",f6.1," (o)")') Vs,Fs
write(*,'("Fazni napon faze T..................................:",f8.1," (V)",f6.1," (o)")') Vt,Ft
write(*,'(/,"Kruznica za proracun raspodjele magnetskog polja:")')
write(*,'("Broj tocaka na kruznici za proracun polja...........:",i5)') Nq
write(*,'("Koordinate sredista kruznice (xs,ys)................:",f8.1," ,",f8.1," (mm)")') xs,ys
write(*,'("Radijus kruznice....................................:",f8.1," (mm)")') Rq
write(*,'(/,"Pravac za proracun raspodjele magnetskog polja:")')
write(*,'("Broj tocaka duz pravca za proracun polja...........:",i5)') Nqk
write(*,'("Koordinate pocetne tocke pravca (xs,ys)............:",f8.1," ,",f8.1," (mm)")') xstart,ystart
write(*,'("Korak proracuna (pravac je paralelan s osi x)......:",f8.1," (mm)")') dx

write(*,'(/,"Priprema ...")')

! ======================================================================
! Formiranje matrice vlastitih i medjusobnih udaljenosti dionih vodica
! ======================================================================
! Dimenzionalno uskladjenje ulaznih velicina
hc = hc * 1.d-3   ! (mm) => (m)
Rc = Rc * 1.d-3   ! (mm) => (m)
hs = hs * 1.d-3   ! (mm) => (m)
Rs = Rs * 1.d-3   ! (mm) => (m)
D0 = D0 * 1.d-3   ! (mm) => (m)

xs = xs * 1.d-3   ! (mm) => (m)
ys = ys * 1.d-3   ! (mm) => (m)
Rq = Rq * 1.d-3   ! (mm) => (m)

xstart = xstart * 1.d-3 ! (mm) => (m)
ystart = ystart * 1.d-3 ! (mm) => (m)
dx = dx * 1.d-3         ! (mm) => (m)

! ---------------------------------
! Proracun koordinata dionih vodica
! ---------------------------------
! 1. sloj faznog vodica
dFic1 = (2.d0*PI)/real(N1)
sc1 = dFic1 * Rc
! 2. sloj faznog vodica
dFic2 = (2.d0*PI)/real(N2)
sc2 = dFic2 * (Rc + hc/3.d0)
! 3. sloj faznog vodica
dFic3 = (2.d0*PI)/real(N3)
sc3 = dFic3 * (Rc + (2.d0/3.d0)*hc)
! Oklop
dFis = (2.d0*PI)/real(Ns)
ss = dFis * Rs

! Ukupni broj svih dionih vodica sustava
Nuk = 3*Nc + 3*Ns

allocate(x(Nuk))
allocate(y(Nuk))

! ---------------------------
! Fazni vodici (sve tri faze)
! ---------------------------
allocate(delta(Nuk))
! 1. sloj
Fic = 0.d0
do i = 1,N1
    delta(i) = Fic
    x(i) = Rc * dcos(Fic)
    y(i) = Rc * dsin(Fic)
    Fic = Fic + dFic1
end do
! 2. sloj
Fic = 0.d0
do i = N1+1,N1+N2
    delta(i) = Fic
    x(i) = (Rc + hc/3.d0) * dcos(Fic)
    y(i) = (Rc + hc/3.d0) * dsin(Fic)
    Fic = Fic + dFic2
end do
! 3. sloj
Fic = 0.d0
do i = N1+N2+1,Nc
    delta(i) = Fic
    x(i) = (Rc + (2.d0/3.d0)*hc) * dcos(Fic)
    y(i) = (Rc + (2.d0/3.d0)*hc) * dsin(Fic)
    Fic = Fic + dFic3
end do
! Ostale faze (svi slojevi)
do i = Nc+1,2*Nc
    delta(i) = delta(i-Nc)
    x(i) = x(i-Nc) + D0
    y(i) = y(i-Nc)
end do
do i = 2*Nc+1,3*Nc
    delta(i) = delta(i-2*Nc)
    x(i) = x(i-2*Nc) + 2.d0*D0
    y(i) = y(i-2*Nc)
end do
! --------------------
! Oklop (sve tri faze)
! --------------------
! Oklop lijeve faze
Fis = 0.d0
do i = 3*Nc+1,3*Nc+Ns
    delta(i) = Fis
    x(i) = Rs * dcos(Fis)
    y(i) = Rs * dsin(Fis)
    Fis = Fis + dFis
end do
! Oklop srednje
do i = 3*Nc+Ns+1,3*Nc+2*Ns
    delta(i) = delta(i-Ns)
    x(i) = x(i-Ns) + D0
    y(i) = y(i-Ns)
end do
! Oklop desne faze
do i = 3*Nc+2*Ns+1,Nuk
    delta(i) = delta(i-2*Ns)
    x(i) = x(i-2*Ns) + 2.d0*D0
    y(i) = y(i-2*Ns)
end do

! ----------------------------------------------
! Tocke na kruznici za proracun magnetskog polja
! ----------------------------------------------
allocate(xq(Nq))
allocate(yq(Nq))
allocate(theta(Nq))
dTHq = (2.d0*PI)/real(Nq)
THq = 0.d0
do k = 1,Nq
    theta(k) = THq
    xq(k) = xs + Rq * cos(THq)
    yq(k) = ys + Rq * sin(THq)
    THq = THq + dTHq
end do

! ----------------------------------------------
! Tocke na pravcu za proracun magnetskog polja
! ----------------------------------------------
allocate(xk(Nqk))
allocate(yk(Nqk))
allocate(thetak(Nqk))
xk(1) = xstart
yk(1) = ystart
do k = 2,Nqk
    thetak(k) = 0.d0
    xk(k) = xk(k-1) + dx
    yk(k) = yk(k-1)
end do

! ------------------------------------------------
! PRIPREMA ZA PRORACUN MAGNETSKOG POLJA - KRUZNICA
! ------------------------------------------------
allocate(Hksi(Nq,Nuk))
allocate(Heta(Nq,Nuk))
DO k = 1,Nq
    do i = 1,Nuk
        if (i <= 3*Nc) then
            ! Fazni vodici
            if (i <= N1) then
                ! 1. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc1,Hk,He)
            else if ((i > N1).and.(i <= N1+N2)) then
                ! 2. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc2,Hk,He)
            else if ((i > N1+N2).and.(i <= Nc)) then
                ! 3. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc3,Hk,He)
            else if ((i > Nc).and.(i <= Nc+N1)) then
                ! 1. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc1,Hk,He)
            else if ((i > Nc+N1).and.(i <= Nc+N1+N2)) then
                ! 2. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc2,Hk,He)
            else if ((i > Nc+N1+N2).and.(i <= 2*Nc)) then
                ! 3. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc3,Hk,He)
            else if ((i > 2*Nc).and.(i <= 2*Nc+N1)) then
                ! 1. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc1,Hk,He)
            else if ((i > 2*Nc+N1).and.(i <= 2*Nc+N1+N2)) then
                ! 2. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc2,Hk,He)
            else if ((i > 2*Nc+N1+N2).and.(i <= 3*Nc)) then
                ! 3. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hc/3.d0,sc3,Hk,He)
            end if
        else
            ! Oklopi
            call mag_polje(delta(i),x(i),y(i),xq(k),yq(k),hs,ss,Hk,He)
        end if
        ! Priprema za izracun komponenti magnetskog 
        ! polja u svakoj od tocaka (kruznice)
        Hksi(k,i) = Hk
        Heta(k,i) = He
    end do
END DO

deallocate(xq)
deallocate(yq)

! ----------------------------------------------
! PRIPREMA ZA PRORACUN MAGNETSKOG POLJA - PRAVAC
! ----------------------------------------------
allocate(Hksi_pr(Nqk,Nuk))
allocate(Heta_pr(Nqk,Nuk))
DO k = 1,Nqk
    do i = 1,Nuk
        if (i <= 3*Nc) then
            ! Fazni vodici
            if (i <= N1) then
                ! 1. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc1,Hkk,Hek)
            else if ((i > N1).and.(i <= N1+N2)) then
                ! 2. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc2,Hkk,Hek)
            else if ((i > N1+N2).and.(i <= Nc)) then
                ! 3. sloj -> 1. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc3,Hkk,Hek)
            else if ((i > Nc).and.(i <= Nc+N1)) then
                ! 1. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc1,Hkk,Hek)
            else if ((i > Nc+N1).and.(i <= Nc+N1+N2)) then
                ! 2. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc2,Hkk,Hek)
            else if ((i > Nc+N1+N2).and.(i <= 2*Nc)) then
                ! 3. sloj -> 2. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc3,Hkk,Hek)
            else if ((i > 2*Nc).and.(i <= 2*Nc+N1)) then
                ! 1. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc1,Hkk,Hek)
            else if ((i > 2*Nc+N1).and.(i <= 2*Nc+N1+N2)) then
                ! 2. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc2,Hkk,Hek)
            else if ((i > 2*Nc+N1+N2).and.(i <= 3*Nc)) then
                ! 3. sloj -> 3. faza
                call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hc/3.d0,sc3,Hkk,Hek)
            end if
        else
            ! Oklopi
            call mag_polje(delta(i),x(i),y(i),xk(k),yk(k),hs,ss,Hkk,Hek)
        end if
        ! Priprema za izracun komponenti magnetskog 
        ! polja u svakoj od tocaka (kruznice)
        Hksi_pr(k,i) = Hkk
        Heta_pr(k,i) = Hek
    end do
END DO

deallocate(xk)
deallocate(yk)

! --------------------------------------------------------------------
! Formiranje matrice vlastitih i medjusobnih udaljenosti dionih vodica
! --------------------------------------------------------------------
allocate(d(Nuk,Nuk))
DO i = 1,Nuk
    do k = i,Nuk
        IF (i == k) THEN
            ! ------------
            ! Vlastite SGU
            ! ------------
            if (i <= 3*Nc) then
                ! Fazni vodici
                if (i <= N1) then
                    ! 1. sloj -> 1. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc1)
                else if ((i > N1).and.(i <= N1+N2)) then
                    ! 2. sloj -> 1. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc2)
                else if ((i > N1+N2).and.(i <= Nc)) then
                    ! 3. sloj -> 1. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc3)
                else if ((i > Nc).and.(i <= Nc+N1)) then
                    ! 1. sloj -> 2. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc1)
                else if ((i > Nc+N1).and.(i <= Nc+N1+N2)) then
                    ! 2. sloj -> 2. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc2)
                else if ((i > Nc+N1+N2).and.(i <= 2*Nc)) then
                    ! 3. sloj -> 2. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc3)
                else if ((i > 2*Nc).and.(i <= 2*Nc+N1)) then
                    ! 1. sloj -> 3. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc1)
                else if ((i > 2*Nc+N1).and.(i <= 2*Nc+N1+N2)) then
                    ! 2. sloj -> 3. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc2)
                else if ((i > 2*Nc+N1+N2).and.(i <= 3*Nc)) then
                    ! 3. sloj -> 3. faza
                    d(i,i) = 0.2235d0 * (hc/3.d0 + sc3)
                end if
            else
                ! Oklop
                d(i,i) = 0.2235d0 * (hs + ss)
            end if
        ELSE
            ! --------------
            ! Medjusobne SGU
            ! --------------
            d(i,k) = dsqrt((x(i)-x(k))**2+(y(i)-y(k))**2)
        END IF
        d(k,i) = d(i,k)  ! simetricno u odnosu na gl. dijagonalu
    end do
END DO

deallocate(x)
deallocate(y)

! ======================================================================
! Proacun ekvivalentnih otpora dionih vodica faza i oklopa
! ======================================================================
allocate(R1(Nuk))
! ------------------
! Fazni vodici
! ------------------
! 1. sloj -> 1. faza
do i = 1,N1
    R1(i) = L / (Kc * (((hc/3.d0)*sc1)*1.d6))
end do
! 2. sloj -> 1. faza
do i = N1+1,N1+N2
    R1(i) = L / (Kc * (((hc/3.d0)*sc2)*1.d6))
end do
! 3. sloj -> 1. faza
do i = N1+N2+1,Nc
    R1(i) = L / (Kc * (((hc/3.d0)*sc3)*1.d6))
end do
! 1. sloj -> 2. faza
do i = Nc+1,Nc+N1
    R1(i) = L / (Kc * (((hc/3.d0)*sc1)*1.d6))
end do
! 2. sloj -> 2. faza
do i = Nc+N1+1,Nc+N1+N2
    R1(i) = L / (Kc * (((hc/3.d0)*sc2)*1.d6))
end do
! 3. sloj -> 2. faza
do i = Nc+N1+N2+1,2*Nc
    R1(i) = L / (Kc * (((hc/3.d0)*sc3)*1.d6))
end do
! 1. sloj -> 3. faza
do i = 2*Nc+1,2*Nc+N1
    R1(i) = L / (Kc * (((hc/3.d0)*sc1)*1.d6))
end do
! 2. sloj -> 3. faza
do i = 2*Nc+N1+1,2*Nc+N1+N2
    R1(i) = L / (Kc * (((hc/3.d0)*sc2)*1.d6))
end do
! 3. sloj -> 3. faza
do i = 2*Nc+N1+N2+1,3*Nc
    R1(i) = L / (Kc * (((hc/3.d0)*sc3)*1.d6))
end do
! ----------------
! Oklop (sve faze)
! ----------------
do i = 3*Nc+1,Nuk
    R1(i) = L / (Ks * ((hs*ss)*1.d6))
end do

call cpu_time(tend1)
time1 = tend1-tstart1

! ######################################################################
! # DALJNJI PRORACUN OVISI O NACINU UZEMLJENJA OKLOPA GENERATORSKIH    #
! # VODOVA (NA OBA KRAJA ILI SAMO NA JEDNOM KRAJU). ZATO SE NASTAVAK   #
! # PRORACUNA PROVODI ODVOJENO, OVISNO O NACINU UZEMLJENJA OKLOPA.     #
! ######################################################################
write(*,'(/,"Odaerite nacin uzemljenja oklopa generatorskih vodova:")')
write(*,'("---------------------------------------------------------")')
write(*,'("1) Oklop je uzemljen na oba kraja")')
write(*,'("2) Oklop je uzemljen samo na jednom kraju")')
write(*,'("3) Oklop je uzemljen samo na jednom kraju i kratko-spojen")')
write(*,'("---------------------------------------------------------")')
write(*,'("Izbor: ")',advance='no')
read(*,*) izbor

! Objasnjenje varijanti nacina tretiranja oklopa (uzemljenja i kratkog
! spajanja) generatorskih oklopljenih vodova: 
! * U slucaju pod (1) oklop je uzemljen na oba kraja (pri cemu moze biti
!   dodatno medjusobno spojen ukratko na pocetku i/ili kraju trase te na
!   proizvoljnom dodatnom broju mjesta duz trase).
! * U slucaju pod (2) oklop je uzemljen samo na jednom kraju (gdje su mu
!   i oklopi kratko spojeni), dok istovremeno nije kratko spojen na ne-
!   uzemljenom kraju. Postoje pritom tri napona (neuzemljenog kraja) ok-
!   lopa prema zemlji.
! * U slucaju pod (3) oklop je ponovno uzemljen samo na jednom kraju uz
!   medjusobno kratko spojanje oklopa na pocetku i na kraju trase (te na
!   proizvoljnom eventualnom broju mjesta duz trase). Postoji samo jedan
!   napon neuzemljenog kraja oklopa (koji su medjusobno kratko spojeni)
!   prema zemlji / uzemljivacu.

call cpu_time(tstart2)

write(*,'(/,"Racunam ...")')

SELECT CASE(izbor)
    CASE(1)
    ! ******************************************************************
    ! *    OKLOPLJENI GENERATORSKI VOD JE UZEMLJEN NA OBA KRAJA        *
    ! ******************************************************************
    ! Formiranje matrice vlastitih i medjusobnih impedancija sustava
    ! ==================================================================
    omega = 2.d0*PI*f
    mio = 4.d0*PI*1.d-7
    K1 = (omega*mio*L)/8.d0
    K2 = (omega*mio*L)/(2.d0*PI)
    De = 658.d0 * dsqrt(ro/f)

    allocate(Z(Nuk+1,Nuk+1))

    DO i = 1,Nuk
        do k = i,Nuk
            IF (i == k) THEN
                ! ---------------------------
                ! Dijagonalni clanovi matrice
                ! ---------------------------
                if (i <= Nc) then
                    ! Elementi faze R
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zrt
                else if ((i > Nc).and.(i <= 2*Nc)) then
                    ! Elementi faze S
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zst
                else if ((i > 2*Nc).and.(i <= 3*Nc)) then
                    ! Elementi faze T
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Ztt
                else if ((i> 3*Nc).and.(k > 3*Nc)) then
                    ! Elementi oklopa
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zgu + Ztu
                end if
            ELSE
                ! ------------------------------
                ! Vandijagonalni clanovi matrice
                ! ------------------------------
                if ((i <= Nc).and.(k <= Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zrt
                else if (((i > Nc).and.(i <= 2*Nc)).AND.((k > Nc).and.(k <= 2*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zst
                else if (((i > 2*Nc).and.(i <= 3*Nc)).AND.((k > 2*Nc).and.(k <= 3*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Ztt
                else if ((i > 3*Nc).and.(k > 3*Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zgu + Ztu
                else
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik
                end if
                Z(k,i) = Z(i,k)  ! Simetricno u odnosu na gl. dijagonalu
            END IF
        end do
    END DO
    deallocate(d)

    ! ---------------------------------------------------------
    ! Dodavanje nula i jedinica u matrici Z u skladu s teorijom
    ! (generatorske oklop. sabirnice su dio izolirane mreze).
    ! ---------------------------------------------------------
    Z(Nuk+1,1:3*Nc) = one
    Z(Nuk+1,3*Nc+1:Nuk) = zero
    Z(1:3*Nc,Nuk+1) = one
    Z(3*Nc+1:Nuk,Nuk+1) = zero
    Z(Nuk+1,Nuk+1) = zero

    ! ==================================================================
    !      FORMIRANJE VEKTORA ZADANIH NAPONA (GENERATORSKA STRANA)
    ! ==================================================================
    Re = Vr * dcos(Fr*(PI/180.d0))
    Im = Vr * dsin(Fr*(PI/180.d0))
    V1 = dcmplx(Re,Im)
    Re = Vs * dcos(Fs*(PI/180.d0))
    Im = Vs * dsin(Fs*(PI/180.d0))
    V2 = dcmplx(Re,Im)
    Re = Vt * dcos(Ft*(PI/180.d0))
    Im = Vt * dsin(Ft*(PI/180.d0))
    V3 = dcmplx(Re,Im)

    allocate(V(Nuk+1))
    V(1:Nc) = V1
    V(Nc+1:2*Nc) = V2
    V(2*Nc+1:3*Nc) = V3
    V(3*Nc+1:Nuk+1) = zero

    ! ==================================================================
    !      RJESAVANJE SUSTAVA LINEARNIH ALGEBARSKIH JEDNADZBI
    ! ==================================================================
    ! Sustav jednadzbi rjesava se koristenjem LAPACK rutine ZSYSV
    UPLO = 'U'
    N = Nuk + 1
    NRHS = 1
    allocate(B(N,NRHS))
    B(:,NRHS) = V       ! Formiranje desne strane sustava
    LDA = N; LDB = N
    allocate(IPIV(N))
    LWORK = 4*N
    allocate(WORK(LWORK))
    ! -------------------------------------------------------
    CALL ZSYSV(UPLO,N,NRHS,Z,LDA,IPIV,B,LDB,WORK,LWORK,INFO )
    ! -------------------------------------------------------
    deallocate(WORK)
    deallocate(IPIV)

    deallocate(Z)
    deallocate(V)

    if (INFO == 0) then
        ! -------------------------------------
        ! Proracun je uspjesno proveden
        ! -------------------------------------
        write(*,'(/,"Proracun uspjesan! -> ZSYSV: INFO = ",i2)') INFO
        ! Prebacivanje rjesenja u vektor struja
        allocate(Ii(Nuk))
        Ii = B(1:Nuk,NRHS)
        DV = B(Nuk+1,NRHS)
        deallocate(B)
        
        ! Struje u slojevima faznog vodica
        I1 = sum(Ii(1:N1))
        I2 = sum(Ii(N1+1:N1+N2))
        I3 = sum(Ii(N1+N2+1:N1+N2+N3))

        ! -------------------------------------
        ! Proracun faznih struja
        ! -------------------------------------
        Ir = sum(Ii(1:Nc))
        Is = sum(Ii(Nc+1:2*Nc))
        It = sum(Ii(2*Nc+1:3*Nc))
        ! -------------------------------------
        ! Proracun struja oklopa
        ! -------------------------------------
        Iro = sum(Ii(3*Nc+1:3*Nc+Ns))
        Iso = sum(Ii(3*Nc+Ns+1:3*Nc+2*Ns))
        Ito = sum(Ii(3*Nc+2*Ns+1:Nuk))
        ! -------------------------------------
        ! Proracun struje u zemlji
        ! -------------------------------------
        Iz = sum(Ii(1:Nuk))

        ! ==================================================================
        ! Proracun gubitaka snage oklopljenog voda
        ! ==================================================================
        ! Faza R:
        Pr = 0.d0
        do i = 1,Nc
            Pr = Pr + cdabs(Ii(i))**2 * R1(i)
        end do
        Pr = Pr/L    ! (W/m)
        ! Faza S:
        Ps = 0.d0
        do i = Nc+1,2*Nc
            Ps = Ps + cdabs(Ii(i))**2 * R1(i)
        end do
        Ps = Ps/L    ! (W/m)
        ! Faza T:
        Pt = 0.d0
        do i = 2*Nc+1,3*Nc
            Pt = Pt + cdabs(Ii(i))**2 * R1(i)
        end do
        Pt = Pt/L    ! (W/m)

        ! Oklop - faza R:
        Por = 0.d0
        do i = 3*Nc+1,3*Nc+Ns
            Por = Por + (cdabs(Ii(i))**2) * R1(i)
        end do
        Por = Por/L    ! (W/m)
        ! Oklop - faza S:
        Pos = 0.d0
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            Pos = Pos + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pos = Pos/L    ! (W/m)
        ! Oklop - faza T:
        Pot = 0.d0
        do i = 3*Nc+2*Ns+1,Nuk
            Pot = Pot + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pot = Pot/L    ! (W/m)

        ! ------------------------------------
        ! Ukupni gubici snage oklopljenog voda
        ! ------------------------------------
        Puc = Pr + Ps + Pt     ! (W/m)
        Pus = Por + Pos + Pot  ! (W/m)
        Puk = Puc + Pus        ! (W/m)

        ! --------------------------------------
        ! Gubici snage oklopa za graficki prikaz
        ! --------------------------------------
        allocate(Pori(Ns))
        allocate(Posi(Ns))
        allocate(Poti(Ns))
        do i = 3*Nc+1,3*Nc+Ns
            k = i - 3*Nc
            Pori(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            k = i - (3*Nc+Ns)
            Posi(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+2*Ns+1,Nuk
            k = i - (3*Nc+2*Ns)
            Poti(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do

        deallocate(R1)

        ! ==================================================================
        ! Proracun nadomjesnih otpora za gubitke snage oklopljenog voda
        ! ==================================================================
        Rphr = Pr/(cdabs(Ir)**2)
        Rphs = Ps/(cdabs(Is)**2)
        Rpht = Pt/(cdabs(It)**2)

        Rokr = Por/(cdabs(Iro)**2)
        Roks = Pos/(cdabs(Iso)**2)
        Rokt = Pot/(cdabs(Ito)**2)

        Reqr = Rphr + Rokr
        Reqs = Rphs + Roks
        Reqt = Rpht + Rokt

        ! -----------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - KRUZNICA
        ! -----------------------------------------------------
        allocate(Htan(Nq))
        allocate(Hrad(Nq))
        allocate(H(Nq))
        Htan(:) = zero
        Hrad(:) = zero
        do k = 1,Nq
            do i = 1,Nuk
                Htan(k) = Htan(k) + Ii(i) * ( Hksi(k,i)*cos(theta(k)-delta(i)) + &
                Heta(k,i)*sin(theta(k)-delta(i)) )
                Hrad(k) = Hrad(k) + Ii(i) * (-Hksi(k,i)*sin(theta(k)-delta(i)) + &
                Heta(k,i)*cos(theta(k)-delta(i)) )
            end do
            H(k) = sqrt(Htan(k)**2+Hrad(k)**2)
        end do
        
        ! Provjera izracunatih vrijednosti pomocu Amperovog zakona
!       csum = zero
!       do k = 1,Nq
!           csum = csum + Htan(k)
!       end do
!       Ampere = -Rq*dTHq * csum
!       Isum = sum(Ii(1:Nuk))  !Ir + Iro

        deallocate(Hksi)
        deallocate(Heta)

        ! ---------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - PRAVAC
        ! ---------------------------------------------------
        allocate(Htan_pr(Nqk))
        allocate(Hrad_pr(Nqk))
        allocate(H_pr(Nqk))
        Htan_pr(:) = zero
        Hrad_pr(:) = zero
        do k = 1,Nqk
            do i = 1,Nuk
                Htan_pr(k) = Htan_pr(k) + Ii(i) * ( Hksi_pr(k,i)*cos(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*sin(thetak(k)-delta(i)) )
                Hrad_pr(k) = Hrad_pr(k) + Ii(i) * (-Hksi_pr(k,i)*sin(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*cos(thetak(k)-delta(i)) )
            end do
            H_pr(k) = sqrt(Htan_pr(k)**2+Hrad_pr(k)**2)
        end do

        deallocate(Hksi_pr)
        deallocate(Heta_pr)

        deallocate(delta)

        ! ---------------------------------------------------------------------
        ! PRORACUN GUSTOCE MAGNETSKOG TOKA OKLOPLJENOG VODA - KRUZNICA I PRAVAC
        ! ---------------------------------------------------------------------
        mio = 4.d0 * PI * 1.d-7
        allocate(Bfc(Nqk))
        allocate(Bfp(Nqk))
        do k = 1,Nq
            Bfc(k) = mio * H(k) * 1.d3     ! (mT)
        end do
        do k = 1,Nqk
            Bfp(k) = mio * H_pr(k) * 1.d3  ! (mT)
        end do

        ! ------------------------------------------
        ! PRORACUN ELEKTROMAGNETSKIH SILA VODICA
        ! ------------------------------------------
        allocate(Fm(Nq))
        Fm(:) = zero
        mio = 4.d0*PI*1.d-7
        do k = 1,Nq
            do i = 1,Nuk
                Fm(k) = Fm(k) + mio*cdabs(Ii(i))*cdabs(H(k))
            end do
        end do

        ! ==================================================================
        ! Ispis rezultata proracuna u file
        ! ==================================================================
        open(unit=2,file='output.txt',status='old',action='write',position='append')
        write(2,'(//,"REZULTATI PRORACUNA:")')
        write(2,'(/,"Oklopljeni generatorski vod je uzemljen na oba kraja!")')
        write(2,'("Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(2,'("Dimenzije dionog vodica (axb):")')
        ac = (hc/3.d0) * 1.d3
        bc1 = (Rc*1.d3) * dFic1
        bc2 = ((Rc*1.d3)+ac) * dFic2
        bc3 = ((Rc*1.d3)+2.d0*ac) * dFic3
        write(2,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(2,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(2,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        as = hs * 1.d3
        bs = (Rs*1.d3) * dFis
        write(2,'("Oklop................:",f6.2," x",f6.2," (mm)")') as,bs

        ! -------------------------------------
        ! Ispis rezultata proracuna na ekran
        ! -------------------------------------
        write(*,'(//,"REZULTATI PRORACUNA:")')
        write(*,'(/,"Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(*,'(/,"Dimenzije dionog vodica (axb):")')
        write(*,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(*,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(*,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        write(*,'("Oklop................:",f6.2," x",f6.2," (mm)")') as,bs
        
        ! ---------------------------------------------
        ! Poziv subroutine za ispis rezultata proracuna
        ! ---------------------------------------------
        call ispis1(I1,I2,I3,Ir,Is,It,Iro,Iso,Ito,Iz,DV,&
        Pr,Ps,Pt,Puc,Por,Pos,Pot,Pus,Puk)

        ! -----------------------------------------------------
        ! Ispis ekvivalentnih otpora za proracun gubitaka snage
        ! -----------------------------------------------------
        write(*,'(/,"Ekvivalentni otpori za proracun gubitaka:")')
        write(*,'("-----------------------------------------")')
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqr
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqs
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqt
        write(*,'("-----------------------------------------")')
        write(*,'("Ukupni ekv. otpor:",es12.4," (ohm/m)")') Reqr+Reqs+Reqt

        ! Potvrda ispravnog proracuna magnetskog polja
!       call zmk(Ampere,Modul,Kut)
!       write(*,'(/,"Amperov zakon:",f8.2,2x,f6.1)') Modul,Kut
!       call zmk(Isum,Modul,Kut)
!       write(*,'(/,"Suma struja:",f8.2,2x,f6.1)') Modul,Kut

        ! ---------------------------------------------
        ! Ispis rezultata - STRUJE - za graficki prikaz
        ! ---------------------------------------------
        call ispis_struje(Nc,Ns,Nuk,dFis,Ii)

        deallocate(Ii)

        ! ---------------------------------------------------
        ! Ispis rezultata - GUBICI SNAGE - za graficki prikaz
        ! ---------------------------------------------------
        call ispis_snaga(Ns,dFis,Pori,Posi,Poti)

        deallocate(Pori)
        deallocate(Posi)
        deallocate(Poti)

        ! -----------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - KRUZNICA
        ! -----------------------------------------------------------------
        call ispis_magpolje(Nq,theta,Htan,Hrad,H)

        deallocate(H)
        deallocate(Htan)
        deallocate(Hrad)

        open(unit=3,file='Bfield_circ.txt',action='write')
        write(3,'(2x,"Fi",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nq
            call zmk(Bfc(k),Modul,Kut)
            write(3,'(f8.3,2x,f10.4,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
        end do
        close(3)
        deallocate(Bfc)

        ! ----------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKE SILE - za graficki prikaz - KRUZNICA
        ! ----------------------------------------------------------------
        open(unit=3,file='EM_Force.txt',action='write')
        write(3,'(2x,"Fi",7x,"F(N/m)")')
        do k = 1,Nq
            write(3,'(f6.2,2x,f10.3)') theta(k)*(180.d0/PI),Fm(k)
        end do
        close(3)

        deallocate(Fm)
        deallocate(theta)

        ! ---------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - PRAVAC
        ! ---------------------------------------------------------------
        allocate(xk(Nqk))
        xk(1) = xstart*1.d3
        do k = 2,Nqk
            xk(k) = xk(k-1) + dx*1.d3
        end do
        open(unit=3,file='Mag_pr_tan.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Ht(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Htan_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr_rad.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Hr(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Hrad_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"H(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(H_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Bfield_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Bfp(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        deallocate(xk)
        deallocate(Bfp)
        deallocate(H_pr)
        deallocate(Htan_pr)
        deallocate(Hrad_pr)
        deallocate(thetak)

    else
        ! ------------------------------------
        ! Proracun nije uspio -> kraj programa
        ! ------------------------------------
        write(*,'(/,"Pogreska: Proracun neuspjesan! -> ZSYSV: INFO = ",i2)') INFO
        deallocate(B)
        deallocate(R1)
        write(*,'(/,"Kraj programa!")')
        STOP
    end if

    CASE(2)
    ! ******************************************************************
    ! *  OKLOPLJENI GENERATORSKI VOD JE UZEMLJEN NA SAMO JEDNOM KRAJU  *
    ! *         (OKLOPI NISU MEDJUSOBNO KRATKO SPOJENI)                *
    ! ******************************************************************
    ! Formiranje matrice vlastitih i medjusobnih impedancija sustava
    ! ==================================================================
    omega = 2.d0*PI*f
    mio = 4.d0*PI*1.d-7
    K1 = (omega*mio*L)/8.d0
    K2 = (omega*mio*L)/(2.d0*PI)
    De = 658.d0 * dsqrt(ro/f)

    allocate(Z(Nuk+4,Nuk+4))

    ! Inicijalizacija matrice s nulama
    Z(1:Nuk+4,1:Nuk+4) = zero

    DO i = 1,Nuk
        do k = i,Nuk
            IF (i == k) THEN
                ! ---------------------------
                ! Dijagonalni clanovi matrice
                ! ---------------------------
                if (i <= Nc) then
                    ! Elementi faze R
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zrt
                else if ((i > Nc).and.(i <= 2*Nc)) then
                    ! Elementi faze S
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zst
                else if ((i > 2*Nc).and.(i <= 3*Nc)) then
                    ! Elementi faze T
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Ztt
                else if ((i> 3*Nc).and.(k > 3*Nc)) then
                    ! Elementi oklopa
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zgu + Ztu
                end if
            ELSE
                ! ------------------------------
                ! Vandijagonalni clanovi matrice
                ! ------------------------------
                if ((i <= Nc).and.(k <= Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zrt
                else if (((i > Nc).and.(i <= 2*Nc)).AND.((k > Nc).and.(k <= 2*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zst
                else if (((i > 2*Nc).and.(i <= 3*Nc)).AND.((k > 2*Nc).and.(k <= 3*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Ztt
                else if ((i > 3*Nc).and.(k > 3*Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zgu + Ztu
                else
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik
                end if
                Z(k,i) = Z(i,k)  ! Simetricno u odnosu na gl. dijagonalu
            END IF
        end do
    END DO
    deallocate(d)

    ! ---------------------------------------------------------
    ! Dodavanje jedinica u matrici Z u skladu s teorijom
    ! (generatorske oklop. sabirnice su dio izolirane mreze).
    ! ---------------------------------------------------------
    Z(Nuk+1,1:3*Nc) = one
    Z(Nuk+2,3*Nc+1:3*Nc+Ns) = one
    Z(Nuk+3,3*Nc+Ns+1:3*Nc+2*Ns) = one
    Z(Nuk+4,3*Nc+2*Ns+1:Nuk) = one
    Z(1:3*Nc,Nuk+1) = one
    Z(3*Nc+1:3*Nc+Ns,Nuk+2) = one
    Z(3*Nc+Ns+1:3*Nc+2*Ns,Nuk+3) = one
    Z(3*Nc+2*Ns+1:Nuk,Nuk+4) = one

    ! ==================================================================
    !      FORMIRANJE VEKTORA ZADANIH NAPONA (GENERATORSKA STRANA)
    ! ==================================================================
    Re = Vr * dcos(Fr*(PI/180.d0))
    Im = Vr * dsin(Fr*(PI/180.d0))
    V1 = dcmplx(Re,Im)
    Re = Vs * dcos(Fs*(PI/180.d0))
    Im = Vs * dsin(Fs*(PI/180.d0))
    V2 = dcmplx(Re,Im)
    Re = Vt * dcos(Ft*(PI/180.d0))
    Im = Vt * dsin(Ft*(PI/180.d0))
    V3 = dcmplx(Re,Im)

    allocate(V(Nuk+4))
    V(1:Nc) = V1
    V(Nc+1:2*Nc) = V2
    V(2*Nc+1:3*Nc) = V3
    V(3*Nc+1:Nuk+4) = zero

    ! ==================================================================
    !      RJESAVANJE SUSTAVA LINEARNIH ALGEBARSKIH JEDNADZBI
    ! ==================================================================
    ! Sustav jednadzbi rjesava se koristenjem LAPACK rutine ZSYSV
    UPLO = 'U'
    N = Nuk + 4
    NRHS = 1
    allocate(B(N,NRHS))
    B(:,NRHS) = V       ! Formiranje desne strane sustava
    LDA = N; LDB = N
    allocate(IPIV(N))
    LWORK = 4*N
    allocate(WORK(LWORK))
    ! -------------------------------------------------------
    CALL ZSYSV(UPLO,N,NRHS,Z,LDA,IPIV,B,LDB,WORK,LWORK,INFO )
    ! -------------------------------------------------------
    deallocate(WORK)
    deallocate(IPIV)

    deallocate(Z)
    deallocate(V)

    if (INFO == 0) then
        ! -------------------------------------
        ! Proracun je uspjesno proveden
        ! -------------------------------------
        write(*,'(/,"Proracun uspjesan! -> ZSYSV: INFO = ",i2)') INFO
        ! Prebacivanje rjesenja u vektor struja
        allocate(Ii(Nuk))
        Ii = B(1:Nuk,NRHS)
        DV = B(Nuk+1,NRHS)
        DR = B(Nuk+2,NRHS)
        DS = B(Nuk+3,NRHS)
        DT = B(Nuk+4,NRHS)
        deallocate(B)

        ! Struje u slojevima faznog vodica
        I1 = sum(Ii(1:N1))
        I2 = sum(Ii(N1+1:N1+N2))
        I3 = sum(Ii(N1+N2+1:N1+N2+N3))

        ! -------------------------------------
        ! Proracun faznih struja
        ! -------------------------------------
        Ir = sum(Ii(1:Nc))
        Is = sum(Ii(Nc+1:2*Nc))
        It = sum(Ii(2*Nc+1:3*Nc))
        ! -------------------------------------
        ! Proracun struja oklopa
        ! -------------------------------------
        Iro = sum(Ii(3*Nc+1:3*Nc+Ns))
        Iso = sum(Ii(3*Nc+Ns+1:3*Nc+2*Ns))
        Ito = sum(Ii(3*Nc+2*Ns+1:Nuk))
        ! -------------------------------------
        ! Proracun struje u zemlji
        ! -------------------------------------
        Iz = sum(Ii(1:Nuk))

        ! ==================================================================
        ! Proracun gubitaka snage oklopljenog voda
        ! ==================================================================
        ! Faza R:
        Pr = 0.d0
        do i = 1,Nc
            Pr = Pr + cdabs(Ii(i))**2 * R1(i)
        end do
        Pr = Pr/L    ! (W/m)
        ! Faza S:
        Ps = 0.d0
        do i = Nc+1,2*Nc
            Ps = Ps + cdabs(Ii(i))**2 * R1(i)
        end do
        Ps = Ps/L    ! (W/m)
        ! Faza T:
        Pt = 0.d0
        do i = 2*Nc+1,3*Nc
            Pt = Pt + cdabs(Ii(i))**2 * R1(i)
        end do
        Pt = Pt/L    ! (W/m)
        ! Oklop - faza R:
        Por = 0.d0
        do i = 3*Nc+1,3*Nc+Ns
            Por = Por + (cdabs(Ii(i))**2) * R1(i)
        end do
        Por = Por/L    ! (W/m)
        ! Oklop - faza S:
        Pos = 0.d0
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            Pos = Pos + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pos = Pos/L    ! (W/m)
        ! Oklop - faza T:
        Pot = 0.d0
        do i = 3*Nc+2*Ns+1,Nuk
            Pot = Pot + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pot = Pot/L    ! (W/m)

        ! ------------------------------------
        ! Ukupni gubici snage oklopljenog voda
        ! ------------------------------------
        Puc = Pr + Ps + Pt     ! (W/m)
        Pus = Por + Pos + Pot  ! (W/m)
        Puk = Puc + Pus        ! (W/m)

        ! --------------------------------------
        ! Gubici snage oklopa za graficki prikaz
        ! --------------------------------------
        allocate(Pori(Ns))
        allocate(Posi(Ns))
        allocate(Poti(Ns))
        do i = 3*Nc+1,3*Nc+Ns
            k = i - 3*Nc
            Pori(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            k = i - (3*Nc+Ns)
            Posi(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+2*Ns+1,Nuk
            k = i - (3*Nc+2*Ns)
            Poti(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do

        deallocate(R1)

        ! ==================================================================
        ! Proracun nadomjesnih otpora za gubitke snage oklopljenog voda
        ! ==================================================================
        Rphr = Pr/(cdabs(Ir)**2)
        Rphs = Ps/(cdabs(Is)**2)
        Rpht = Pt/(cdabs(It)**2)

        Rokr = Por/(cdabs(Iro)**2)
        Roks = Pos/(cdabs(Iso)**2)
        Rokt = Pot/(cdabs(Ito)**2)

        Reqr = Rphr + Rokr
        Reqs = Rphs + Roks
        Reqt = Rpht + Rokt

        ! -----------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - KRUZNICA
        ! -----------------------------------------------------
        allocate(Htan(Nq))
        allocate(Hrad(Nq))
        allocate(H(Nq))
        Htan(:) = zero
        Hrad(:) = zero
        do k = 1,Nq
            do i = 1,Nuk
                Htan(k) = Htan(k) + Ii(i) * ( Hksi(k,i)*cos(theta(k)-delta(i)) + &
                Heta(k,i)*sin(theta(k)-delta(i)) )
                Hrad(k) = Hrad(k) + Ii(i) * (-Hksi(k,i)*sin(theta(k)-delta(i)) + &
                Heta(k,i)*cos(theta(k)-delta(i)) )
            end do
            H(k) = sqrt(Htan(k)**2+Hrad(k)**2)
        end do
        
        deallocate(Hksi)
        deallocate(Heta)

        ! ---------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - PRAVAC
        ! ---------------------------------------------------
        allocate(Htan_pr(Nqk))
        allocate(Hrad_pr(Nqk))
        allocate(H_pr(Nqk))
        Htan_pr(:) = zero
        Hrad_pr(:) = zero
        do k = 1,Nqk
            do i = 1,Nuk
                Htan_pr(k) = Htan_pr(k) + Ii(i) * ( Hksi_pr(k,i)*cos(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*sin(thetak(k)-delta(i)) )
                Hrad_pr(k) = Hrad_pr(k) + Ii(i) * (-Hksi_pr(k,i)*sin(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*cos(thetak(k)-delta(i)) )
            end do
            H_pr(k) = sqrt(Htan_pr(k)**2+Hrad_pr(k)**2)
        end do

        deallocate(Hksi_pr)
        deallocate(Heta_pr)

        deallocate(delta)

        ! ---------------------------------------------------------------------
        ! PRORACUN GUSTOCE MAGNETSKOG TOKA OKLOPLJENOG VODA - KRUZNICA I PRAVAC
        ! ---------------------------------------------------------------------
        mio = 4.d0 * PI * 1.d-7
        allocate(Bfc(Nqk))
        allocate(Bfp(Nqk))
        do k = 1,Nq
            Bfc(k) = mio * H(k) * 1.d3     ! (mT)
        end do
        do k = 1,Nqk
            Bfp(k) = mio * H_pr(k) * 1.d3  ! (mT)
        end do

        ! ------------------------------------------
        ! PRORACUN ELEKTROMAGNETSKIH SILA VODICA
        ! ------------------------------------------
        allocate(Fm(Nq))
        Fm(:) = zero
        mio = 4.d0*PI*1.d-7
        do k = 1,Nq
            do i = 1,Nuk
                Fm(k) = Fm(k) + mio*cdabs(Ii(i))*cdabs(H(k))
            end do
        end do

        ! -------------------------------------
        ! Ispis rezultata proracuna u file
        ! -------------------------------------
        open(unit=2,file='output.txt',status='old',action='write',position='append')
        write(2,'(//,"REZULTATI PRORACUNA:")')
        write(2,'(/,"Oklopljeni generatorski vod je uzemljen samo na jednom kraju!")')
        write(2,'("Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(2,'("Dimenzije dionog vodica (axb):")')
        ac = (hc/3.d0) * 1.d3
        bc1 = (Rc*1.d3) * dFic1
        bc2 = ((Rc*1.d3)+ac) * dFic2
        bc3 = ((Rc*1.d3)+2.d0*ac) * dFic3
        write(2,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(2,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(2,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        as = hs * 1.d3
        bs = (Rs*1.d3) * dFis
        write(2,'("Oklop......:",f6.2," x",f6.2," (mm)")') as,bs

        ! -------------------------------------
        ! Ispis rezultata proracuna na ekran
        ! -------------------------------------
        write(*,'(//,"REZULTATI PRORACUNA:")')
        write(*,'(/,"Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(*,'(/,"Dimenzije dionog vodica (axb):")')
        write(*,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(*,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(*,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        write(*,'("Oklop................:",f6.2," x",f6.2," (mm)")') as,bs

        ! ---------------------------------------------
        ! Poziv subroutine za ispis rezultata proracuna
        ! ---------------------------------------------
        call ispis2(I1,I2,I3,Ir,Is,It,Iro,Iso,Ito,Iz,DV,DR,DS,DT,&
        Pr,Ps,Pt,Puc,Por,Pos,Pot,Pus,Puk)

        ! -----------------------------------------------------
        ! Ispis ekvivalentnih otpora za proracun gubitaka snage
        ! -----------------------------------------------------
        write(*,'(/,"Ekvivalentni otpori za proracun gubitaka:")')
        write(*,'("-----------------------------------------")')
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqr
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqs
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqt
        write(*,'("-----------------------------------------")')
        write(*,'("Ukupni ekv. otpor:",es12.4," (ohm/m)")') Reqr+Reqs+Reqt

        ! ---------------------------------------------
        ! Ispis rezultata - STRUJA - za graficki prikaz
        ! ---------------------------------------------
        call ispis_struje(Nc,Ns,Nuk,dFis,Ii)

        deallocate(Ii)

        ! ---------------------------------------------------
        ! Ispis rezultata - GUBICI SNAGE - za graficki prikaz
        ! ---------------------------------------------------
        call ispis_snaga(Ns,dFis,Pori,Posi,Poti)

        deallocate(Pori)
        deallocate(Posi)
        deallocate(Poti)

        ! -----------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - KRUZNICA
        ! -----------------------------------------------------------------
        call ispis_magpolje(Nq,theta,Htan,Hrad,H)

        deallocate(H)
        deallocate(Htan)
        deallocate(Hrad)

        open(unit=3,file='Bfield_circ.txt',action='write')
        write(3,'(2x,"Fi",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nq
            call zmk(Bfc(k),Modul,Kut)
            write(3,'(f6.2,2x,f10.3,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
        end do
        close(3)
        deallocate(Bfc)

        ! ------------------------------------------------------
        ! Ispis rezultata - MAGNETSKE SILE - za graficki prikaz
        ! ------------------------------------------------------
        open(unit=3,file='EM_Force.txt',action='write')
        write(3,'(2x,"Fi",7x,"F(N/m)")')
        do k = 1,Nq
            write(3,'(f6.2,2x,f10.3)') theta(k)*(180.d0/PI),Fm(k)
        end do
        close(3)

        deallocate(theta)
        deallocate(Fm)

        ! ---------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - PRAVAC
        ! ---------------------------------------------------------------
        allocate(xk(Nqk))
        xk(1) = xstart*1.d3
        do k = 2,Nqk
            xk(k) = xk(k-1) + dx*1.d3
        end do
        open(unit=3,file='Mag_pr_tan.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Ht(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Htan_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr_rad.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Hr(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Hrad_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"H(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(H_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Bfield_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Bfp(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        deallocate(xk)
        deallocate(Bfp)
        deallocate(H_pr)
        deallocate(Htan_pr)
        deallocate(Hrad_pr)
        deallocate(thetak)

    else
        ! ------------------------------------
        ! Proracun nije uspio -> kraj programa
        ! ------------------------------------
        write(*,'(/,"Pogreska: Proracun neuspjesan! -> ZSYSV: INFO = ",i2)') INFO
        deallocate(B)
        deallocate(R1)
        write(*,'(/,"Kraj programa!")')
        STOP
    end if

    CASE(3)
    ! ******************************************************************
    ! *  OKLOPLJENI GENERATORSKI VOD JE UZEMLJEN NA SAMO JEDNOM KRAJU  *
    ! *           (OKLOPI SU MEDJUSOBNO KRATKO SPOJENI)                *
    ! ******************************************************************
    ! Formiranje matrice vlastitih i medjusobnih impedancija sustava
    ! ==================================================================
    omega = 2.d0*PI*f
    mio = 4.d0*PI*1.d-7
    K1 = (omega*mio*L)/8.d0
    K2 = (omega*mio*L)/(2.d0*PI)
    De = 658.d0 * dsqrt(ro/f)

    allocate(Z(Nuk+2,Nuk+2))

    ! Inicijalizacija matrice s nulama
    Z(:,:) = zero

    DO i = 1,Nuk
        do k = i,Nuk
            IF (i == k) THEN
                ! ---------------------------
                ! Dijagonalni clanovi matrice
                ! ---------------------------
                if (i <= Nc) then
                    ! Elementi faze R
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zrt
                else if ((i > Nc).and.(i <= 2*Nc)) then
                    ! Elementi faze S
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zst
                else if ((i > 2*Nc).and.(i <= 3*Nc)) then
                    ! Elementi faze T
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Ztt
                else if ((i> 3*Nc).and.(k > 3*Nc)) then
                    ! Elementi oklopa
                    Re = R1(i) + K1
                    Im = K2 * dlog(De/d(i,i))
                    Zii = dcmplx(Re,Im)
                    Z(i,i) = Zii + Zgu + Ztu
                end if
            ELSE
                ! ------------------------------
                ! Vandijagonalni clanovi matrice
                ! ------------------------------
                if ((i <= Nc).and.(k <= Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zrt
                else if (((i > Nc).and.(i <= 2*Nc)).AND.((k > Nc).and.(k <= 2*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zst
                else if (((i > 2*Nc).and.(i <= 3*Nc)).AND.((k > 2*Nc).and.(k <= 3*Nc))) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Ztt
                else if ((i > 3*Nc).and.(k > 3*Nc)) then
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik + Zgu + Ztu
                else
                    Re = K1
                    Im = K2 * dlog(De/d(i,k))
                    Zik = dcmplx(Re,Im)
                    Z(i,k) = Zik
                end if
                Z(k,i) = Z(i,k)  ! Simetricno u odnosu na gl. dijagonalu
            END IF
        end do
    END DO
    deallocate(d)

    ! ---------------------------------------------------------
    ! Dodavanje jedinica u matrici Z u skladu s teorijom
    ! (generatorske oklop. sabirnice su dio izolirane mreze).
    ! ---------------------------------------------------------
    Z(Nuk+1,1:3*Nc) = one
    Z(1:3*Nc,Nuk+1) = one

    Z(Nuk+2,3*Nc+1:3*Nc+3*Ns) = one
    Z(3*Nc+1:3*Nc+3*Ns,Nuk+2) = one

    ! ==================================================================
    !      FORMIRANJE VEKTORA ZADANIH NAPONA (GENERATORSKA STRANA)
    ! ==================================================================
    Re = Vr * dcos(Fr*(PI/180.d0))
    Im = Vr * dsin(Fr*(PI/180.d0))
    V1 = dcmplx(Re,Im)
    Re = Vs * dcos(Fs*(PI/180.d0))
    Im = Vs * dsin(Fs*(PI/180.d0))
    V2 = dcmplx(Re,Im)
    Re = Vt * dcos(Ft*(PI/180.d0))
    Im = Vt * dsin(Ft*(PI/180.d0))
    V3 = dcmplx(Re,Im)

    allocate(V(Nuk+2))
    V(:) = zero
    V(1:Nc) = V1
    V(Nc+1:2*Nc) = V2
    V(2*Nc+1:3*Nc) = V3

    ! ==================================================================
    !      RJESAVANJE SUSTAVA LINEARNIH ALGEBARSKIH JEDNADZBI
    ! ==================================================================
    ! Sustav jednadzbi rjesava se koristenjem LAPACK rutine ZSYSV
    UPLO = 'U'
    N = Nuk + 2
    NRHS = 1
    allocate(B(N,NRHS))
    B(:,NRHS) = V       ! Formiranje desne strane sustava
    LDA = N; LDB = N
    allocate(IPIV(N))
    LWORK = 4*N
    allocate(WORK(LWORK))
    ! -------------------------------------------------------
    CALL ZSYSV(UPLO,N,NRHS,Z,LDA,IPIV,B,LDB,WORK,LWORK,INFO )
    ! -------------------------------------------------------
    deallocate(WORK)
    deallocate(IPIV)

    deallocate(Z)
    deallocate(V)

    if (INFO == 0) then
        ! -------------------------------------
        ! Proracun je uspjesno proveden
        ! -------------------------------------
        write(*,'(/,"Proracun uspjesan! -> ZSYSV: INFO = ",i3)') INFO
        ! Prebacivanje rjesenja u vektor struja
        allocate(Ii(Nuk))
        Ii = B(1:Nuk,NRHS)
        DV = B(Nuk+1,NRHS)
        DV2 = B(Nuk+2,NRHS)
        deallocate(B)

        ! Struje u slojevima faznog vodica
        I1 = sum(Ii(1:N1))
        I2 = sum(Ii(N1+1:N1+N2))
        I3 = sum(Ii(N1+N2+1:N1+N2+N3))

        ! -------------------------------------
        ! Proracun faznih struja
        ! -------------------------------------
        Ir = sum(Ii(1:Nc))
        Is = sum(Ii(Nc+1:2*Nc))
        It = sum(Ii(2*Nc+1:3*Nc))
        ! -------------------------------------
        ! Proracun struja oklopa
        ! -------------------------------------
        Iro = sum(Ii(3*Nc+1:3*Nc+Ns))
        Iso = sum(Ii(3*Nc+Ns+1:3*Nc+2*Ns))
        Ito = sum(Ii(3*Nc+2*Ns+1:Nuk))
        ! -------------------------------------
        ! Proracun struje u zemlji
        ! -------------------------------------
        Iz = sum(Ii(1:Nuk))

        ! ==================================================================
        ! Proracun gubitaka snage oklopljenog voda
        ! ==================================================================
        ! Faza R:
        Pr = 0.d0
        do i = 1,Nc
            Pr = Pr + cdabs(Ii(i))**2 * R1(i)
        end do
        Pr = Pr/L    ! (W/m)
        ! Faza S:
        Ps = 0.d0
        do i = Nc+1,2*Nc
            Ps = Ps + cdabs(Ii(i))**2 * R1(i)
        end do
        Ps = Ps/L    ! (W/m)
        ! Faza T:
        Pt = 0.d0
        do i = 2*Nc+1,3*Nc
            Pt = Pt + cdabs(Ii(i))**2 * R1(i)
        end do
        Pt = Pt/L    ! (W/m)
        ! Oklop - faza R:
        Por = 0.d0
        do i = 3*Nc+1,3*Nc+Ns
            Por = Por + (cdabs(Ii(i))**2) * R1(i)
        end do
        Por = Por/L    ! (W/m)
        ! Oklop - faza S:
        Pos = 0.d0
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            Pos = Pos + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pos = Pos/L    ! (W/m)
        ! Oklop - faza T:
        Pot = 0.d0
        do i = 3*Nc+2*Ns+1,Nuk
            Pot = Pot + (cdabs(Ii(i))**2) * R1(i)
        end do
        Pot = Pot/L    ! (W/m)

        ! ------------------------------------
        ! Ukupni gubici snage oklopljenog voda
        ! ------------------------------------
        Puc = Pr + Ps + Pt     ! (W/m)
        Pus = Por + Pos + Pot  ! (W/m)
        Puk = Puc + Pus        ! (W/m)

        ! --------------------------------------
        ! Gubici snage oklopa za graficki prikaz
        ! --------------------------------------
        allocate(Pori(Ns))
        allocate(Posi(Ns))
        allocate(Poti(Ns))
        do i = 3*Nc+1,3*Nc+Ns
            k = i - 3*Nc
            Pori(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+Ns+1,3*Nc+2*Ns
            k = i - (3*Nc+Ns)
            Posi(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do
        do i = 3*Nc+2*Ns+1,Nuk
            k = i - (3*Nc+2*Ns)
            Poti(k) = ((cdabs(Ii(i))**2) * R1(i))/L  ! (W/m)
        end do

        deallocate(R1)

        ! ==================================================================
        ! Proracun nadomjesnih otpora za gubitke snage oklopljenog voda
        ! ==================================================================
        Rphr = Pr/(cdabs(Ir)**2)
        Rphs = Ps/(cdabs(Is)**2)
        Rpht = Pt/(cdabs(It)**2)

        Rokr = Por/(cdabs(Iro)**2)
        Roks = Pos/(cdabs(Iso)**2)
        Rokt = Pot/(cdabs(Ito)**2)

        Reqr = Rphr + Rokr
        Reqs = Rphs + Roks
        Reqt = Rpht + Rokt

        ! -----------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - KRUZNICA
        ! -----------------------------------------------------
        allocate(Htan(Nq))
        allocate(Hrad(Nq))
        allocate(H(Nq))
        Htan(:) = zero
        Hrad(:) = zero
        do k = 1,Nq
            do i = 1,Nuk
                Htan(k) = Htan(k) + Ii(i) * ( Hksi(k,i)*cos(theta(k)-delta(i)) + &
                Heta(k,i)*sin(theta(k)-delta(i)) )
                Hrad(k) = Hrad(k) + Ii(i) * (-Hksi(k,i)*sin(theta(k)-delta(i)) + &
                Heta(k,i)*cos(theta(k)-delta(i)) )
            end do
            H(k) = sqrt(Htan(k)**2+Hrad(k)**2)
        end do
        
        deallocate(Hksi)
        deallocate(Heta)

        ! ---------------------------------------------------
        ! PRORACUN MAGNETSKOG POLJA OKLOPLJENOG VODA - PRAVAC
        ! ---------------------------------------------------
        allocate(Htan_pr(Nqk))
        allocate(Hrad_pr(Nqk))
        allocate(H_pr(Nqk))
        Htan_pr(:) = zero
        Hrad_pr(:) = zero
        do k = 1,Nqk
            do i = 1,Nuk
                Htan_pr(k) = Htan_pr(k) + Ii(i) * ( Hksi_pr(k,i)*cos(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*sin(thetak(k)-delta(i)) )
                Hrad_pr(k) = Hrad_pr(k) + Ii(i) * (-Hksi_pr(k,i)*sin(thetak(k)-delta(i)) + &
                Heta_pr(k,i)*cos(thetak(k)-delta(i)) )
            end do
            H_pr(k) = sqrt(Htan_pr(k)**2+Hrad_pr(k)**2)
        end do

        deallocate(Hksi_pr)
        deallocate(Heta_pr)

        deallocate(delta)

        ! ---------------------------------------------------------------------
        ! PRORACUN GUSTOCE MAGNETSKOG TOKA OKLOPLJENOG VODA - KRUZNICA I PRAVAC
        ! ---------------------------------------------------------------------
        mio = 4.d0 * PI * 1.d-7
        allocate(Bfc(Nqk))
        allocate(Bfp(Nqk))
        do k = 1,Nq
            Bfc(k) = mio * H(k) * 1.d3     ! (mT)
        end do
        do k = 1,Nqk
            Bfp(k) = mio * H_pr(k) * 1.d3  ! (mT)
        end do

        ! ------------------------------------------
        ! PRORACUN ELEKTROMAGNETSKIH SILA VODICA
        ! ------------------------------------------
        allocate(Fm(Nq))
        Fm(:) = zero
        mio = 4.d0*PI*1.d-7
        do k = 1,Nq
            do i = 1,Nuk
                Fm(k) = Fm(k) + mio*cdabs(Ii(i))*cdabs(H(k))
            end do
        end do

        ! -------------------------------------
        ! Ispis rezultata proracuna u file
        ! -------------------------------------
        open(unit=2,file='output.txt',status='old',action='write',position='append')
        write(2,'(//,"REZULTATI PRORACUNA:")')
        write(2,'(/,"Oklopljeni generatorski vod je uzemljen samo na jednom kraju (oklopi kratko-spojeni)!")')
        write(2,'("Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(2,'("Dimenzije dionog vodica (axb):")')
        ac = (hc/3.d0) * 1.d3
        bc1 = (Rc*1.d3) * dFic1
        bc2 = ((Rc*1.d3)+ac) * dFic2
        bc3 = ((Rc*1.d3)+2.d0*ac) * dFic3
        write(2,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(2,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(2,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        as = hs * 1.d3
        bs = (Rs*1.d3) * dFis
        write(2,'("Oklop......:",f6.2," x",f6.2," (mm)")') as,bs

        ! -------------------------------------
        ! Ispis rezultata proracuna na ekran
        ! -------------------------------------
        write(*,'(//,"REZULTATI PRORACUNA:")')
        write(*,'(/,"Ukupni broj dionih vodica: Nuk = ",i5)') Nuk
        write(*,'(/,"Dimenzije dionog vodica (axb):")')
        write(*,'("Fazni vodic - 1. sloj:",f6.2," x",f6.2," (mm)")') ac,bc1
        write(*,'("Fazni vodic - 2. sloj:",f6.2," x",f6.2," (mm)")') ac,bc2
        write(*,'("Fazni vodic - 3. sloj:",f6.2," x",f6.2," (mm)")') ac,bc3
        write(*,'("Oklop................:",f6.2," x",f6.2," (mm)")') as,bs

        ! ---------------------------------------------
        ! Poziv subroutine za ispis rezultata proracuna
        ! ---------------------------------------------
        call ispis3(I1,I2,I3,Ir,Is,It,Iro,Iso,Ito,Iz,DV,DV2,&
        Pr,Ps,Pt,Puc,Por,Pos,Pot,Pus,Puk)

        ! -----------------------------------------------------
        ! Ispis ekvivalentnih otpora za proracun gubitaka snage
        ! -----------------------------------------------------
        write(*,'(/,"Ekvivalentni otpori za proracun gubitaka:")')
        write(*,'("-----------------------------------------")')
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqr
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqs
        write(*,'("Faza R (vodic + oklop):",es12.4," (ohm/m)")') Reqt
        write(*,'("-----------------------------------------")')
        write(*,'("Ukupni ekv. otpor:",es12.4," (ohm/m)")') Reqr+Reqs+Reqt

        ! ---------------------------------------------
        ! Ispis rezultata - STRUJA - za graficki prikaz
        ! ---------------------------------------------
        call ispis_struje(Nc,Ns,Nuk,dFis,Ii)

        deallocate(Ii)

        ! ---------------------------------------------------
        ! Ispis rezultata - GUBICI SNAGE - za graficki prikaz
        ! ---------------------------------------------------
        call ispis_snaga(Ns,dFis,Pori,Posi,Poti)

        deallocate(Pori)
        deallocate(Posi)
        deallocate(Poti)

        ! -----------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - KRUZNICA
        ! -----------------------------------------------------------------
        call ispis_magpolje(Nq,theta,Htan,Hrad,H)

        deallocate(H)
        deallocate(Htan)
        deallocate(Hrad)

        open(unit=3,file='Bfield_circ.txt',action='write')
        write(3,'(2x,"Fi",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nq
            call zmk(Bfc(k),Modul,Kut)
            write(3,'(f6.2,2x,f10.3,2x,f6.1)') theta(k)*(180.d0/PI),Modul,Kut
        end do
        close(3)
        deallocate(Bfc)

        ! ------------------------------------------------------
        ! Ispis rezultata - MAGNETSKE SILE - za graficki prikaz
        ! ------------------------------------------------------
        open(unit=3,file='EM_Force.txt',action='write')
        write(3,'(2x,"Fi",7x,"F(N/m)")')
        do k = 1,Nq
            write(3,'(f6.2,2x,f10.3)') theta(k)*(180.d0/PI),Fm(k)
        end do
        close(3)

        deallocate(theta)
        deallocate(Fm)

        ! ---------------------------------------------------------------
        ! Ispis rezultata - MAGNETSKO POLJE - za graficki prikaz - PRAVAC
        ! ---------------------------------------------------------------
        allocate(xk(Nqk))
        xk(1) = xstart*1.d3
        do k = 2,Nqk
            xk(k) = xk(k-1) + dx*1.d3
        end do
        open(unit=3,file='Mag_pr_tan.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Ht(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Htan_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr_rad.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"Hr(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Hrad_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Mag_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"H(A/m)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(H_pr(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        open(unit=3,file='Bfield_pr.txt',action='write')
        write(3,'(2x,"x(mm)",7x,"B(mT)",5x,"(o)")')
        do k = 1,Nqk
            call zmk(Bfp(k),Modul,Kut)
            write(3,'(f8.2,2x,f10.4,2x,f6.1)') xk(k),Modul,Kut
        end do
        close(3)
        deallocate(xk)
        deallocate(Bfp)
        deallocate(H_pr)
        deallocate(Htan_pr)
        deallocate(Hrad_pr)
        deallocate(thetak)

    else
        ! ------------------------------------
        ! Proracun nije uspio -> kraj programa
        ! ------------------------------------
        write(*,'(/,"Pogreska: Proracun neuspjesan! -> ZSYSV: INFO = ",i2)') INFO
        deallocate(B)
        deallocate(R1)
        write(*,'(/,"Kraj programa!")')
        STOP
    end if

    CASE DEFAULT
        write(*,'("Pogresan unos! -> Kraj programa.")')
        deallocate(B)
        deallocate(R1)
        STOP
END SELECT

call cpu_time(tend2)
time2 = tend2-tstart2

trun = time1 + time2

write(2,'(/,"Vrijeme trajanja proracuna:")')
write(2,'("--------------------------")')
write(2,'("Priprema:",f8.1," (sec)")') time1
write(2,'("Proracun:",f8.1," (sec)")') time2
write(2,'("--------------------------")')
write(2,'("Ukupno vrijeme:",f8.1," (sec)")') trun
close(2)
write(*,'(/,"Vrijeme trajanja proracuna:")')
write(*,'("--------------------------")')
write(*,'("Priprema:",f8.1," (sec)")') time1
write(*,'("Proracun:",f8.1," (sec)")') time2
write(*,'("--------------------------")')
write(*,'("Ukupno vrijeme:",f8.1," (sec)")') trun

write(*,'(/,"Kraj programa!",/)')

END PROGRAM
