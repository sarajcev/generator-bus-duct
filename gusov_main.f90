! ------------------------------------------------------------------------
!        *** GUBICI SNAGE OKLOPLJENIH GENERATORSKIH VODOVA ***
! ------------------------------------------------------------------------
! Proracun gubitaka snage (Jouleovi gubici) trofaznog sustava pravokutnih
! oklopljenih sabirnica za evakuaciju snage iz kaverne (strojarnice) hidro-
! elektrane. Matematicki model temelji se na metodi dionih vodica i SGU
! iz doktorske disertacije Prof. I. Sarajceva. Oklopljeni vod se sastoji
! od tri pravokutne sabirnice na odredjenom razmaku unutar pravokutnog
! oklopa odredjene debljine, koji moze ili ne mora imati pregrade. Koristi
! se originalno razvijena metoda dionih vodica (Prof. I. Sarajcev). Racunaju
! se vlastite i medjusobne udaljenosti izmedju dionih vodica (pravokutnih
! segmenata) vodica pravokutnog poprecnog presjeka. Na temelju njih formira
! se sustav jednadzbi, rjesenjem kojeg se dobivaju struje dionih vodica.
! Ove struje su mjerodavne za proracun Jouleovih gubitaka i topline.
! ------------------------------------------------------------------------

program GuSOV
implicit none

integer tip
real(8) XP,YP
real(8) La,Lb,Dx
integer Na,Nb,Ni,Nu
real(8) a,b
real(8) Xo,Yo
real(8) Du,Ho,Hx,h
integer N1,N2,Nuk
real(8) aa,bb
real(8) kapav,kapao,Sp
real(8) Rr,Xr,Rs,Xs,Rt,Xt
real(8) Rgu,Xgu,Rtu,Xtu
real(8) Vr,Fir,Vs,Fis,Vt,Fit
real(8) ro,f,L

real(8),parameter :: pi = 3.1415926535897d0
complex(8),parameter :: one = dcmplx(1.d0,0.d0)
complex(8),parameter :: zero = dcmplx(0.d0,0.d0)
real(8),dimension(:),allocatable :: x,y
real(8),dimension(:),allocatable :: Rp
real(8),dimension(:,:),allocatable :: d
complex(8),dimension(:,:),allocatable :: Z
complex(8) zii,zik
complex(8) Zr,Zs,Zt,Zgu,Ztu,Zuu
real(8) omega,mio,Konst1,Konst2,De
complex(8) V1,V2,V3
complex(8),dimension(:),allocatable :: V
complex(8) Ifr,Ifs,Ift,Ifo,Iz
real(8) Pr,Ps,Pt,Po
real(8) Rer,Res,Ret,Reo
complex(8) Sgr,Sgs,Sgt,Sg,Str,Sts,Stt,St,DS
real(8) MOD,KUT
real(8) Re,Im
real(8) delta
integer idelta
integer i,j,k,m,n,s,p,q

character(1) UPLO
integer MN
integer NRHS,LDA,LDB,LWORK,INFO
complex(8),dimension(:,:),allocatable :: BA
complex(8),dimension(:),allocatable :: WORK
integer,dimension(:),allocatable :: IPIV
integer izbor
complex(8) I1,I2,I3,I4,I5,I6
integer kk
real(8) tstart,tend,time

! Opis ulaznih varijabli:
! tip - tip oklopljenog voda (1 ili 2)
!       1 - segregated (oklopljeni vod koji ima pregrade izmedju
!           faznih vodica), jednopolno izolirani vodici u zajednickom
!           oklopu,
!       2 - non-segregated (oklopljeni vod koji se sastoji od tri fazna
!           vodica u zajednickom oklopu, bez pregrada).
! La - duljina donje stranice pravokutnog faznog vodica oklopljenih
!      sabirnica (mm),
! Lb - duljina bocne stranice pravokutnog faznog vodica oklopljenih
!      sabirnica (mm),
! Na - broj podjela donje stranice pravokutnika na dione vodice,
! Nb - broj podjela bocne stranice pravokutnika na dione vodice,
!      Ukupni broj dionih vodica jednog faznog vodica stoga iznosi
!      Na*Nb. Postoje tri fazna vodica te ce stoga ukupni broj dionih
!      vodica sviju faza biti 3*(Na*Nb).
! kapav - specificna elektricna vodljivost materijala od kojeg je
!         izradjen fazni vodic na radnoj temperaturi (Sm/mm2).
! Xp - X-koordinata donjeg lijevog ugla prvog faznog vodica (m),
! Yp - Y-koordinata donjeg lijevog ugla prvog faznog vodica (m),
! Dx - razmak medju fazama pravokutnih oklopljenih sabirnica (mm),
! Du - duljina donje stranice pravokutnog oklopa (mm),
! Ho - duljina bocne stranice pravokutnog oklopa (mm),
! Hx - udaljenost od lijeve bocne stranice oklopa do prve pregrade (mm),
! N1 - broj podjela donje stranice oklopa na dione vodice,
! N2 - broj podjela bocne stranice oklopa na dione vodice,
! Xo - X-koordinata donjeg lijevog ugla pravokutnog oklopa (m) -> Xo = 0,
! Yo - Y-koordinata donjeg lijevog ugla pravokutnog oklopa (m) -> Yo = 0,
! h  - debljina oklopa (mm),
! kapao - specificna elektricna vodljivost materijala od kojeg je
!         izradjen oklop na radnoj temperaturi (Sm/mm2).
! L -  duljina oklopljenog voda (m),
! ro - specificni elektricni otpor tla (ohm*m)
! Rr - djelatni otpor faze R (L1) mreze na teretnoj strani (ohm),
! Xr - reaktancija faze R (L1) mreze na strani tereta (ohm),
! Rs - djelatni otpor faze S (L2) mreze na teretnoj strani (ohm),
! Xs - reaktancija faze S (L2) mreze na strani tereta (ohm),
! Rt - djelatni otpor faze T (L3) mreze na teretnoj strani (ohm),
! Xt - reaktancija faze T (L3) mreze na strani tereta (ohm),
! Rgu - djelatni otpor uzemljenja uzemljivaca na generatorskoj strani (ohm),
! Xgu - reaktancija uzemljenja uzemljivaca na generatorskoj strani (ohm),
! Rtu - djelatni otpor uzemljenja uzemljivaca na teretnoj strani (ohm),
! Xtu - reaktancija uzemljenja uzemljivaca na teretnoj strani (ohm),
! Vr  - efektivna vrijednost faznog napona faze R (L1) mreze na
!       generatorskoj strani (V),
! Fir - kut faznog napona faze R (L1), (o),
! Vs  - efektivna vrijednost faznog napona faze S (L2) mreze na
!       generatorskoj strani (V),
! Fis - kut faznog napona faze S (L2), (o),
! Vt  - efektivna vrijednost faznog napona faze T (L3) mreze na
!       generatorskoj strani (V),
! Fit - kut faznog napona faze T (L3), (o),
! f - frekvencija napona (Hz).

! Primjer forme ulazne datoteke:
! ------------------------------
! FILE: sgu_in.txt
! ------------------------------
! tip (1 ili 2)
! La(mm)   Lb(mm)    Na     Nb
! Xp(mm)   Yp(mm)    Dx(mm) kapav
! Du(mm)   Ho(mm)    Hx(mm) N1     N2
! Xo=0(m)  Yo=0(m)   h(mm)  kapao
! L(m)     ro(ohm*m) f(Hz)
! Rr(ohm)  Xr(ohm)
! Rs(ohm)  Xs(ohm)
! Rt(ohm)  Xt(ohm)
! Rgu(ohm) Xgu(ohm)
! Rtu(ohm) Xtu(ohm)
! Vr(V)    Fir(o)
! Vs(V)    Fis(o)
! Vt(V)    Fit(o)
! ------------------------------
! Napomena: Varijable se zapisuju u slobodnom formatu.

! Primjer graficke ilustracije oklopljenog voda
! tip 1 - segregated (oklop ima neku debljinu)

!   --------------------------
!  |        |        |        |
!  |   --   |   --   |   --   |
!  |  |  |  |  |  |  |  |  |  |
!  |  |  |  |  |  |  |  |  |  |
!  |  |  |  |  |  |  |  |  |  |
!  |   --   |   --   |   --   |
!  |        |        |        |
!   --------------------------

! Primjer graficke ilustracije oklopljenog voda
! tip 2 - non-segregated (oklop ima neku debljinu)

!   --------------------------
!  |                          |
!  |   --       --       --   |
!  |  |  |     |  |     |  |  |
!  |  |  |     |  |     |  |  |
!  |  |  |     |  |     |  |  |
!  |   --       --       --   |
!  |                          |
!   --------------------------

! Dimenzije oklopljenog voda zadaju se u milimetrima,
! u skladu s prethodnim opisom ulaznih varijabli.


! ----------------------------------------------------------------------
!              CITANJE I ISPIS ULAZNIH PODATAKA
! ----------------------------------------------------------------------
write(*,'(" PROGRAM: GuSOV-3")')
write(*,'(" Proracun gubitaka snage pravokutnih oklopljenih vodova")')

open(unit=1,file='sgu_in.txt')
read(1,*) tip
read(1,*) La,Lb,Na,Nb
read(1,*) Xp,Yp,Dx,kapav
read(1,*) Du,Ho,Hx,N1,N2
read(1,*) Xo,Yo,h,kapao
read(1,*) L,ro,f
read(1,*) Rr,Xr
read(1,*) Rs,Xs
read(1,*) Rt,Xt
read(1,*) Rgu,Xgu
read(1,*) Rtu,Xtu
read(1,*) Vr,Fir
read(1,*) Vs,Fis
read(1,*) Vt,Fit
close(1)

open(unit=2,file='sgu_out.txt')
write(2,'(" -----------------------------------------------------------------")')
write(2,'(" PROGRAM: GuSOV-3")')
write(2,'(" -----------------------------------------------------------------")')
write(2,'(" Proracun gubitaka snage pravokutnih  oklopljenih vodova. Racunaju")')
write(2,'(" se  gubici  snage u oklopljenim  vodovima koji spajaju generatore")')
write(2,'(" s blok transformatorima. Rijec je o tropolno oklopljenim vodovima")')
write(2,'(" pravokutnog  poprecnog  presjeka,  s  pregradama  izmedju  faznih")')
write(2,'(" vodica, takodjer  pravokutnog  poprecnog  presjeka. Oklop je  pak")')
write(2,'(" uzemljen  na  oba  kraja (generatorska  strana i strana  tereta).")')
write(2,'(" Moze  se  razmatrati  bilo  koje  pogonsko  stanje, ukljucujuci i")')
write(2,'(" nastupe kratkih spojeva.")')
write(2,'(" -----------------------------------------------------------------")')
write(2,'(" Dr. sc. Petar Sarajcev, doc. (petar.sarajcev@fesb.hr)")')
write(2,'(" -----------------------------------------------------------------")')

write(2,'(//," ULAZNI PODACI:")')

write(2,'(/," Tip oklopljenog voda:",i2)') tip
write(2,'(" 1 - oklopljeni vod s pregradama  izmedju faznih vodica (segregated)")')
write(2,'(" 2 - oklopljeni vod s vodicima u zajednickom oklopu (non-segregated)")')

write(2,'(/," Fazni vodici:")')
write(2,'(" Duljina donje stranice pravokutnog faznog vodica oklopljenog voda: ",f8.2," (mm)")') La
write(2,'(" Duljina bocne stranice pravokutnog faznog vodica oklopljenog voda: ",f8.2," (mm)")') Lb
write(2,'(" Broj podjela donje stranice pravokutnika na dione vodice ........: ",i4)') Na
write(2,'(" Broj podjela bocne stranice pravokutnika na dione vodice ........: ",i4)') Nb
write(2,'(" Vodljivost materijala faznog vodica na radnoj temperaturi .......: ",f8.4," (Sm/mm2)")') kapav
write(2,'(" X-koordinata donjeg lijevog ugla prvog faznog vodica ............: ",f8.2," (mm)")') Xp
write(2,'(" Y-koordinata donjeg lijevog ugla prvog faznog vodica ............: ",f8.2," (mm)")') Yp
write(2,'(" Razmak medju fazama pravokutnih oklopljenih sabirnica ...........: ",f8.2," (mm)")') Dx

write(2,'(/," Oklop:")')
write(2,'(" Duljina donje stranice pravokutnog oklopa .......................: ",f8.2," (mm)")') Du
write(2,'(" Duljina bocne stranice pravokutnog oklopa .......................: ",f8.2," (mm)")') Ho
if (tip==1) then
	write(2,'(" Razmak medju stranicama faza pravokutnog oklopa .................: ",f8.2," (mm)")') Hx
end if
write(2,'(" Broj podjela donje stranice oklopa na dione vodice ..............: ",i4)') N1
write(2,'(" Broj podjela bocne stranice oklopa na dione vodice ..............: ",i4)') N2
write(2,'(" Vodljivost materijala oklopa na radnoj temperaturi ..............: ",f8.4," (Sm/mm2)")') kapao
write(2,'(" X-koordinata donjeg lijevog ugla oklopa .........................: ",f8.2," (mm)")') Xo
write(2,'(" Y-koordinata donjeg lijevog ugla oklopa .........................: ",f8.2," (mm)")') Yo
write(2,'(" Debljina oklopa .................................................: ",f8.2," (mm)")') h
write(2,'(" Duljina oklopljenog voda ........................................: ",f8.2," (m)")') L

write(2,'(/," Impedancije na strani tereta:")')
write(2,'(" Faza R (L1): Rr =",f8.2," +j",f8.2," (ohm)")') Rr,Xr
write(2,'(" Faza S (L2): Rs =",f8.2," +j",f8.2," (ohm)")') Rs,Xs
write(2,'(" Faza T (L3): Rt =",f8.2," +j",f8.2," (ohm)")') Rt,Xt
write(2,'(" Relativni specificni elektricni otpor tla .......................: ",f8.2," (ohm*m)")') ro

write(2,'(/," Impedancije uzemljivaca:")')
write(2,'(" Generatorska strana: Zgu =",f8.2," +j",f8.2," (ohm)")') Rgu,Xgu
write(2,'(" Teretna strana ....: Ztu =",f8.2," +j",f8.2," (ohm)")') Rtu,Xtu

write(2,'(/," Fazni naponi:")')
write(2,'(" Faza R (L1):  Vr =",f10.2," (V) /",f6.1," (o)")') Vr,Fir
write(2,'(" Faza S (L2):  Vs =",f10.2," (V) /",f6.1," (o)")') Vs,Fis
write(2,'(" Faza T (L3):  Vt =",f10.2," (V) /",f6.1," (o)")') Vt,Fit
write(2,'(" Frekvencija narinutog napona ....................................: ",f6.1," (Hz)")') f

write(2,'(//," REZULTATI PRORACUNA:",/)')

write(*,'(" Racunam ...")')

call cpu_time(tstart)

! ----------------------------------------------------------------------
!       FORMIRANJE LISTE KOORDINATA DIONIH VODICA FAZA L1,L2,L3
! ----------------------------------------------------------------------
La = La*1.e-3        !(mm) => (m)
Lb = Lb*1.e-3        !(mm) => (m)
Xp = Xp*1.e-3        !(mm) => (m)
YP = Yp*1.e-3        !(mm) => (m)
Dx = Dx*1.e-3        !(mm) => (m)

Du = Du*1.e-3        !(mm) => (m)
Ho = Ho*1.e-3        !(mm) => (m)
Hx = Hx*1.e-3        !(mm) => (m)
h = h*1.e-3          !(mm) => (m)

kapao = kapao*1.e6   !(Sm/mm2) => (Sm/m2)
kapav = kapav*1.e6   !(Sm/mm2) => (Sm/m2)

a = La/real(Na)      !Duljina stranice pravokutnika dionog vodica
b = Lb/real(Nb)      !Duljina stranice pravokutnika dionog vodica

Ni = Na*Nb           !Ukupni broj dionih vodica jedne faze (L1)
Nu = 3*Ni            !Ukupni broj dionih vodica svih faza

! ===================================
! Alociranje memorije za koordinate
! ===================================
  Nuk = Nu + 2*N1 + 2*N2 + 4
  if (tip==1) Nuk = Nuk + 2*N2
! -----------------------------------
  allocate(x(Nuk))
  allocate(y(Nuk))
! -----------------------------------
! Alociranje memorije za povrsine
! dionih vodica
  allocate(Rp(Nuk))
! ===================================
! Alociranje memorije za SGU udalj.
! ===================================
  allocate(d(Nuk,Nuk))
! -----------------------------------

! =================================================
! Formiranje liste koordinata dionih vodica faze L1
! =================================================
x(1) = Xp + a/2.
y(1) = Yp + b/2.
do i = 2,Nb
	x(i) = x(i-1)
	y(i) = y(i-1) + b
end do
do k = 1,Na-1
	m = k*Nb + 1
	n = (k+1)*Nb
	do i = m,n
		x(i) = x(i-Nb) + a
		y(i) = y(i-Nb)
	end do
end do
! =================================================
! Formiranje liste koordinata dionih vodica faze L2
! =================================================
delta = Dx+La
do i = Ni+1,2*Ni
	x(i) = x(i-Ni) + delta
	y(i) = y(i-Ni)
end do
! =================================================
! Formiranje liste koordinata dionih vodica faze L3
! =================================================
do i = 2*Ni+1,3*Ni
	x(i) = x(i-Ni) + delta
	y(i) = y(i-Ni)
end do

! =================================
! Izracun otpora dionih vodica faza
! =================================
do i = 1,Nu
	Sp = a*b
	Rp(i) = L/(kapav*Sp)
end do


! ***************************************
! Ispis broja dionih vodica faza i oklopa
! ***************************************
write(2,'(" Broj dionih vodica faza ...................:",i4)') Nu
if (tip==1) then
	write(2,'(" Broj dionih vodica oklopa .................:",i4)') 2*N1 + 4*N2 + 4
else if (tip==2) then
	write(2,'(" Broj dionih vodica oklopa .................:",i4)') 2*N1 + 2*N2 + 4
end if
write(2,'(" Ukupni broj svih dionih vodica ............:",i4)') Nuk
write(2,'(" Duljina stranice (a) dionog vodica faze ...:",f8.2," (mm)")') a*1.e3
write(2,'(" Duljina stranice (b) dionog vodica faze ...:",f8.2," (mm)")') b*1.e3


! ----------------------------------------------------------------------
!         FORMIRANJE LISTE KOORDINATA DIONIH VODICA OKLOPA
! ----------------------------------------------------------------------
! Duljina stranice pravokutnika dionog vodica oklopa
aa = (Du-2*h)/real(N1)
bb = (Ho-2*h)/real(N2)

! Donji lijevi kut oklopa
k = Nu + 1
x(k) = Xo + h/2.
y(k) = Yo + h/2.

! Donja stranica oklopa - od lijevo prema desno
do i = k+1,N1+k
	x(i) = x(i-1) + aa
	y(i) = y(k)
end do

! Donji desni kut oklopa
m = N1 + k + 1
x(m) = x(m-1) + aa/2. + h/2.
y(m) = y(m-1)

! Lijeva stranica oklopa (prvi clan) - od dolje prema gore
x(m+1) = x(k)
y(m+1) = y(k) + bb/2. + h/2.

! Lijeva stranica oklopa (ostatak) - od dolje prema gore
do i = m+2,N2+m
	x(i) = x(k)
	y(i) = y(i-1) + bb
end do

! Gornji lijevi kut oklopa
n = N2 + m + 1
x(n) = x(k)
y(n) = y(n-1) + bb/2. + h/2.

! Gornja stranica oklopa - od lijevo prema desno
idelta = N1 + N2 + 2
do i = n+1,N1+n
	x(i) = x(i-idelta)
	y(i) = y(k) + Ho - h
end do

! Gornji desni kut oklopa
s = N1 + n + 1
x(s) = x(s-1) + aa/2. + h/2.
y(s) = y(s-1)

! Desna stranica oklopa - od dolje prema gore
idelta = N1 + N2 + 2
do i = s+1,N2+s
	x(i) = x(m)
	y(i) = y(i-idelta)
end do

if (tip==1) then
	! Prva (lijeva) pregrada oklopa - od dolje prema gore
	p = N2 + s
	do i = p+1,N2+p
		x(i) = Xo + Hx + h/2.
		y(i) = y(i-N2)
	end do

	! Druga pregrada oklopa - od dolje prema gore
	q = N2 + p
	do i = q+1,N2+q
		x(i) = Xo + 2*Hx + h/2.
		y(i) = y(i-N2)
	end do
end if

! ==============================================
! Izracun jedinicnog otpora dionih vodica oklopa
! ==============================================
do i = k,Nuk
	if (i==k) then
		Sp = h*h
		Rp(i) = L/(kapao*Sp)

	else if(i>k .and. i<=N1+k) then
		Sp = aa*h
		Rp(i) = L/(kapao*Sp)

	else if(i == m) then
		Sp = h*h
		Rp(i) = L/(kapao*Sp)

	else if(i>m .and. i<=N2+m) then
		Sp = bb*h
		Rp(i) = L/(kapao*Sp)

	else if(i == n) then
		Sp = h*h
		Rp(i) = L/(kapao*Sp)

	else if(i>n .and. i<=N1+n) then
		Sp = aa*h
		Rp(i) = L/(kapao*Sp)

	else if(i == s) then
		Sp = h*h
		Rp(i) = L/(kapao*Sp)

	else if(i>s .and. i<=N2+s) then
		Sp = bb*h
		Rp(i) = L/(kapao*Sp)
	end if

	if (tip==1) then
		if(i>p .and. i<=N2+p) then
			Sp = bb*h
			Rp(i) = L/(kapao*Sp)

		else if(i>q .and. i<=N2+q) then
			Sp = bb*h
			Rp(i) = L/(kapao*Sp)
		end if
	end if
end do

write(2,'(" Duljina stranice (a) dionog vodica oklopa..:",f8.2," (mm)")') aa*1.e3
write(2,'(" Duljina stranice (b) dionog vodica oklopa..:",f8.2," (mm)")') bb*1.e3
write(2,'(" Debljina oklopa ...........................:",f8.2," (mm)")') h*1.e3

! =======================================================
! Formiranje matrice SGU svih dionih vodica faza i oklopa
! =======================================================
do i = 1,Nu
	d(i,i) = 0.2235*(a+b)
end do

do i = k,Nuk
	if (i==k) then
		d(i,i) = 0.2235*(h+h)

	else if(i>k .and. i<=N1+k) then
		d(i,i) = 0.2235*(aa+h)

	else if(i == m) then
		d(i,i) = 0.2235*(h+h)

	else if(i>m .and. i<=N2+m) then
		d(i,i) = 0.2235*(bb+h)

	else if(i == n) then
		d(i,i) = 0.2235*(h+h)

	else if(i>n .and. i<=N1+n) then
		d(i,i) = 0.2235*(aa+h)

	else if(i == s) then
		d(i,i) = 0.2235*(h+h)

	else if(i>s .and. i<=N2+s) then
		d(i,i) = 0.2235*(bb+h)
	end if

	if (tip==1) then
		if(i>p .and. i<=N2+p) then
			d(i,i) = 0.2235*(bb+h)

		else if(i>q .and. i<=N2+q) then
			d(i,i) = 0.2235*(bb+h)
		end if
	end if
end do

do i = 1,Nuk
	do j = i,Nuk
		if (i/=j) then
			d(i,j) = sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
		end if
	end do
end do

! ------------------------------------------------------------------
!          FORMIRANJE MATRICE IMPEDANCIJA SVIH DIONIH VODICA
! ------------------------------------------------------------------
 allocate(Z(Nuk+1,Nuk+1))
! =========================================================
! Definiranje impedancija na strani tereta oklopljenog voda
! =========================================================
Zr = dcmplx(Rr,Xr)
Zs = dcmplx(Rs,Xs)
Zt = dcmplx(Rt,Xt)

! Definiranje impedancija uzemljenje na generatorskoj i
! teretnoj strani oklopljenog voda
Zgu = dcmplx(Rgu,Xgu)
Ztu = dcmplx(Rtu,Xtu)
Zuu = Zgu + Ztu

! Izracun konstanti:
omega = 2.d0*pi*f
mio = 4.d0*pi*1.e-7
Konst1 = (omega*mio*L)/8.d0
Konst2 = (omega*mio*L)/(2.d0*pi)
De = 658.d0*dsqrt(ro/f)

! ---------------------------------
! Popunjavanje matrice impedancija:
! ---------------------------------
do i = 1,Nuk
	Re = Rp(i) + Konst1
	Im = Konst2 * dlog(De/d(i,i))
	zii = dcmplx(Re,Im)

	if (i<=Ni) then
		Z(i,i) = zii + Zr
	else if(i>Ni .and. i<=2*Ni) then
		Z(i,i) = zii + Zs
	else if(i>2*Ni .and. i<=3*Ni) then
		Z(i,i) = zii + Zt
	else
		Z(i,i) = zii + Zuu
	end if
end do

DO i = 1,Nuk
	do j = i,Nuk
		IF (i/=j) THEN
			Re = Konst1
			Im = Konst2 * dlog(De/d(i,j))
			zik = dcmplx(Re,Im)
			if (i<=Ni .AND. j<=Ni) then
				Z(i,j) = zik + Zr
			else if ((i>Ni .and. i<=2*Ni) .AND. (j>Ni .and. j<=2*Ni)) then
				Z(i,j) = zik + Zs
			else if ((i>2*Ni .and. i<=3*Ni) .AND. (j>2*Ni .and. j<=3*Ni)) then
				Z(i,j) = zik + Zt
			else if (i>3*Ni .AND. j>3*Ni) then
				Z(i,j) = zik + Zuu
			else
				Z(i,j) = zik
			end if
		END IF
	end do
END DO

! -----------------------------
! Oslobadjanje zauzete memorije
! -----------------------------
 deallocate(x,y)
 deallocate(d)

! ---------------------------------------------
! Dodavanje nula i jedinica u skladu s teorijom
! (doktorska disertacija Prof. I. Sarajcev)
! ---------------------------------------------
Z(Nuk+1,1:3*Ni) = one
Z(Nuk+1,3*Ni+1:Nuk) = zero
Z(1:3*Ni,Nuk+1) = one
Z(3*Ni+1:Nuk,Nuk+1) = zero
Z(Nuk+1,Nuk+1) = zero

! ----------------------------------------------------------------------
!      FORMIRANJE VEKTORA ZADANIH NAPONA (GENERATORSKA STRANA)
! ----------------------------------------------------------------------
Re = Vr * dcos(Fir*(pi/180.d0))
Im = Vr * dsin(Fir*(pi/180.d0))
V1 = dcmplx(Re,Im)
Re = Vs * dcos(Fis*(pi/180.d0))
Im = Vs * dsin(Fis*(pi/180.d0))
V2 = dcmplx(Re,Im)
Re = Vt * dcos(Fit*(pi/180.d0))
Im = Vt * dsin(Fit*(pi/180.d0))
V3 = dcmplx(Re,Im)

! -----------------------------
! Alociranje memorije za vektor
! -----------------------------
allocate(V(Nuk+1))
!Inicijalizacija
V = zero
! Popunjavanje vektora s vrijednostima napona
V(1:Ni) = V1
V(Ni+1:2*Ni) = V2
V(2*Ni+1:3*Ni) = V3

! ----------------------------------------------------------------------
!   RJESAVANJE SUSTAVA LINEARNIH ALGEBARSKIH JEDNADZBI: {V} = [z]*{I}
! ----------------------------------------------------------------------
! Poziva se LAPACK subroutina ZSYSV koja rjesava simetricni sustav
! kompleksnih linearnih algebarskih jednadzbi u dvostrukoj preciznosti
! ==========================================
! Priprema za poziv LAPACK subroutine ZSYSV:
! ==========================================
UPLO = 'U'
MN = Nuk+1
NRHS = 1
LDA = MN
LDB = MN
LWORK = 4*MN
allocate(WORK(LWORK))
allocate(IPIV(MN))
allocate(BA(LDB,NRHS))
BA(:,NRHS) = V

! ------------------------------
! Poziv LAPACK subroutine ZSYSV:
! ------------------------------
CALL ZSYSV(UPLO,MN,NRHS,Z,LDA,IPIV,BA,LDB,WORK,LWORK,INFO)

if (INFO==0) then
	write(*,'(" ZSYSV: Rjesenje sustava jednadzbi je u redu!")')
	write(*,'(" Rezultati proracuna su u izlaznoj datoteci.")')
else
	write(*,'(" Pogreska pri rjesavanju sustava jednadzbi!")')
	write(*,'(" Kraj programa!")')
	stop
end if

! -----------------------------
! Oslobadjanje zauzete memorije
! -----------------------------
deallocate(Z)
deallocate(V)
deallocate(WORK)
deallocate(IPIV)

! =========================================================
! Definiranje raspodjele struja u faznim vodicima i oklopu
! =========================================================
Ifr = zero
do i = 1,Ni
	Ifr = Ifr + BA(i,NRHS)
end do

Ifs = zero
do i = Ni+1,2*Ni
	Ifs = Ifs + BA(i,NRHS)
end do

Ift = zero
do i = 2*Ni+1,3*Ni
	Ift = Ift + BA(i,NRHS)
end do

Ifo = zero
do i = 3*Ni+1,Nuk
	Ifo = Ifo + BA(i,NRHS)
end do

Iz = zero
do i = 1,Nuk
	Iz = Iz + BA(i,NRHS)
end do

! *********************************************************
! Ispis ukupnih struja u faznim vodicima i struje u oklopu:
! *********************************************************
write(2,'(/," ---------------------------------")')
write(2,'(" Struje u faznim vodicima:")')
write(2,'(" ---------------------------------")')
CALL ZMK(Ifr,MOD,KUT)
write(2,'(" IR =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(Ifs,MOD,KUT)
write(2,'(" IS =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(Ift,MOD,KUT)
write(2,'(" IT =",f10.3," (A) /",f6.1," (o)")') MOD,KUT

write(2,'(/," ---------------------------------")')
write(2,'(" Struja u oklopu:")')
write(2,'(" ---------------------------------")')
CALL ZMK(Ifo,MOD,KUT)
write(2,'(" Io =",f10.3," (A) /",f6.1," (o)")') MOD,KUT

write(2,'(/," ---------------------------------")')
write(2,'(" Struja u zemlji:")')
write(2,'(" ---------------------------------")')
CALL ZMK(Iz,MOD,KUT)
write(2,'(" Iz =",f10.3," (A) /",f6.1," (o)")') MOD,KUT

write(2,'(/," ---------------------------------")')
write(2,'(" Razlika napona medju zvjezdistima")')
write(2,'(" na strani generatora i tereta:")')
write(2,'(" ---------------------------------")')
call ZMK(BA(MN,NRHS),MOD,KUT)
write(2,'(" Vz =",f10.3," (V) /",f6.1," (o)")') MOD,KUT


! ==============================================================
!        Provjera proracuna gubitaka snage
! ==============================================================
Sgr = V1 * dconjg(Ifr)
Sgs = V2 * dconjg(Ifs)
Sgt = V3 * dconjg(Ift)
Sg = Sgr + Sgs + Sgt

Str = (Ifr * Zr) * dconjg(Ifr)
Sts = (Ifs * Zs) * dconjg(Ifs)
Stt = (Ift * Zt) * dconjg(Ift)
St = Str + Sts + Stt

DS = Sg - St

write(2,'(/," ---------------------------------")')
write(2,'(" Provjera proracuna (DS = DP + jDQ):")')
write(2,'(" ---------------------------------")')
write(2,'(" DS =",f12.2," +j",f12.2," (W)")') dreal(DS),dimag(DS)


! ==============================================================
! Definiranje raspodjele struja u pojedinim dijelovima oklopa
! ==============================================================
! Donja stranica oklopa (ukljucuje i kuteve)
I1 = zero
do i = k,N1+k+1
	I1 = I1 + BA(i,NRHS)
end do
! Lijeva stranica oklopa (bez kuteva)
I2 = zero
do i = m+1,N2+m
	I2 = I2 + BA(i,NRHS)
end do
! Gornja stranica oklopa (ukljucuje i kuteve)
I3 = zero
do i = n,N1+n+1
	I3 = I3 + BA(i,NRHS)
end do
! Desna stranica oklopa (bez kuteva)
I4 = zero
do i = s+1,N2+s
	I4 = I4 + BA(i,NRHS)
end do
! Lijeva pregrada (bez kuteva)
I5 = zero
do i = p+1,N2+p
	I5 = I5 + BA(i,NRHS)
end do
! Desna pregrada (bez kuteva)
I6 = zero
do i = q+1,N2+q
	I6 = I6 + BA(i,NRHS)
end do

write(2,'(/," -----------------------------------------------")')
write(2,'(" Rspodjela struje u pojedinim dijelovima oklopa:")')
write(2,'(" -----------------------------------------------")')
CALL ZMK(I1,MOD,KUT)
write(2,'(" Donja  stranica oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(I2,MOD,KUT)
write(2,'(" Lijeva stranica oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(I3,MOD,KUT)
write(2,'(" Gornja stranica oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(I4,MOD,KUT)
write(2,'(" Desna  stranica oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(I5,MOD,KUT)
write(2,'(" Lijeva pregrada oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT
CALL ZMK(I6,MOD,KUT)
write(2,'(" Desna  pregrada oklopa: I =",f10.3," (A) /",f6.1," (o)")') MOD,KUT


! ==============================================================
!        Proracun gubitaka snage oklopljenog voda
! ==============================================================
Pr = 0.d0
do i = 1,Ni
	Pr = Pr + zabs(BA(i,NRHS))**2 * Rp(i)
end do

Ps = 0.d0
do i = Ni+1,2*Ni
	Ps = Ps + zabs(BA(i,NRHS))**2 * Rp(i)
end do

Pt = 0.d0
do i = 2*Ni+1,3*Ni
	Pt = Pt + zabs(BA(i,NRHS))**2 * Rp(i)
end do

Po = 0.d0
do i = 3*Ni+1,Nuk
	Po = Po + zabs(BA(i,NRHS))**2 * Rp(i)
end do

! *********************************************************
! Ispis rezultata proracuna gubitaka snage oklopljenog voda
! *********************************************************
write(2,'(/," ---------------------------------")')
write(2,'(" Gubici snage oklopljenog voda:")')
write(2,'(" ---------------------------------")')
write(2,'(" Faza R:  PR =",f10.2," (W)")') Pr
write(2,'(" Faza S:  PS =",f10.2," (W)")') Ps
write(2,'(" Faza T:  PT =",f10.2," (W)")') Pt
write(2,'(" Oklop :  Po =",f10.2," (W)")') Po
write(2,'(" Ukupni gubici snage oklopljenog voda:",f10.2," (W)")') Pr+Ps+Pt+Po

Pr = Pr/L
Ps = Ps/L
Pt = Pt/L
Po = Po/L

! Ispis rezultata proracuna jedinicnih gubitaka snage oklopljenog voda
write(2,'(/," ----------------------------------------")')
write(2,'(" Jedinicni gubici snage oklopljenog voda:")')
write(2,'(" ----------------------------------------")')
write(2,'(" Faza R:  PR1 =",f10.3," (W/m)")') Pr
write(2,'(" Faza S:  PS1 =",f10.3," (W/m)")') Ps
write(2,'(" Faza T:  PT1 =",f10.3," (W/m)")') Pt
write(2,'(" Oklop :  Po1 =",f10.3," (W/m)")') Po
write(2,'(" Ukupni jedinicni gubici snage oklopljenog voda:",f10.2," (W/m)")') Pr+Ps+Pt+Po

Rer = Pr/(zabs(Ifr)**2)
Res = Ps/(zabs(Ifs)**2)
Ret = Pt/(zabs(Ift)**2)
Reo = Po/(zabs(Ifo)**2)

! Ispis rezultata proracuna jedinicnih ekvivalentnih nadomjesnih
! otpora mjerodavnih za gubitke snage (Jouleove gubitke) prijenosa
! razmatranog oklopljenog generatorskog voda
write(2,'(/," -------------------------------------------------------")')
write(2,'(" Ekvivalentni djelatni otpori gubitaka oklopljenog voda:")')
write(2,'(" -------------------------------------------------------")')
write(2,'(" Faza R:  ReqR =",f10.4," (Ohm/m)")') Rer
write(2,'(" Faza S:  ReqS =",f10.4," (Ohm/m)")') Res
write(2,'(" Faza T:  ReqT =",f10.4," (Ohm/m)")') Ret
write(2,'(" Oklop :  Reqo =",f10.4," (Ohm/m)")') Reo

close(2)


! -----------------------------
! Oslobadjanje zauzete memorije
! -----------------------------
  deallocate(BA)
  deallocate(Rp)

call cpu_time(tend)
time = tend-tstart
write(*,'(" Vrijeme trajanja proracuna iznosi:",f10.3," (sek.)")') time

write(*,'(/," Kraj programa!",/)')

end program
