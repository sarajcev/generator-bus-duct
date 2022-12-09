! ...
! Priprema za proracun magnetskog polja u bilo
! kojoj tocki promatranja (kruznice)
subroutine mag_polje(delta,xi,yi,xq,yq,bi,ai,Hk,He)
	implicit none
	real(8),intent(in) :: delta,xi,yi,xq,yq,ai,bi
	real(8),intent(out) :: Hk,He
	real(8),parameter :: PI = 3.141592653589793d0
	real(8) Ko,P1,P2,P3,P4
	real(8) ksi,eta
	
	
	ksi = (xq-xi)*sin(delta) - (yq-yi)*cos(delta)
	eta = (xq-xi)*cos(delta) + (yq-yi)*sin(delta)

	Ko = 1.d0/(PI*ai*bi)
	P1 = 0.5d0*(ksi+ai/2.d0)*dlog(((eta+bi/2.d0)**2+(ksi+ai/2.d0)**2)/((eta-bi/2.d0)**2+(ksi+ai/2.d0)**2))
	P2 = 0.5d0*(ksi-ai/2.d9)*dlog(((eta+bi/2.d0)**2+(ksi-ai/2.d0)**2)/((eta-bi/2.d0)**2+(ksi-ai/2.d0)**2))
	P3 = (eta+bi/2.d0)*(datan((ksi+ai/2.d0)/(eta+bi/2.d0))-datan((ksi-ai/2.d0)/(eta+bi/2.d0)))
	P4 = (eta-bi/2.d0)*(datan((ksi+ai/2.d0)/(eta-bi/2.d0))-datan((ksi-ai/2.d0)/(eta-bi/2.d0)))

	Hk = -Ko * (P1 - P2 + P3 - P4)

	P1 = 0.5d0*(eta+bi/2.d0)*dlog(((ksi+ai/2.d0)**2+(eta+bi/2.d0)**2)/((ksi-ai/2.d0)**2+(eta+bi/2.d0)**2))
	P2 = 0.5d0*(eta-bi/2.d9)*dlog(((ksi+ai/2.d0)**2+(eta-bi/2.d0)**2)/((ksi-ai/2.d0)**2+(eta-bi/2.d0)**2))
	P3 = (ksi+ai/2.d0)*(datan((eta+bi/2.d0)/(ksi+ai/2.d0))-datan((eta-bi/2.d0)/(ksi+ai/2.d0)))
	P4 = (ksi-ai/2.d0)*(datan((eta+bi/2.d0)/(ksi-ai/2.d0))-datan((eta-bi/2.d0)/(ksi-ai/2.d0)))

	He = Ko * (P1 - P2 + P3 - P4)

end subroutine