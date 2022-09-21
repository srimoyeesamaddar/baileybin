;@/mnt/snoesci/snoe/other/see/level3/get_see
; Newbins.pro

; SMB 8/25/94  ---> 10/94  ----> 1/95
; SMB 8/2002 Added modes 1 and 2, mode 2 hardwired for NRL spectra
;     Mode 1 - Everthing is relative to Hinteregger SC21REFW
;     Mode 2 - New Reference spectra on its own grid is used, for now
;              NRLEUV
; SMB 9/2002 - Added cross sections below 18A from Henke et al 1982,
;              Atomic Data and Nuclear Data Tables, 27, 1, 1982.
; SMB 5/2018 - modified to work with NRL generated bins, romove older stuff like Woods Rocket spectra

mode=2

; This procedure will be used to rewrite solar flux arrays, pure absorption and
; ionization cross sections and ionization state branching ratios. Flux weighting of the cross sections is accomplished with the sc21refw
; reference spectrum.

; In general mod* refers to vectors in the new model bins.

; definitions:

;at_n2 = abosrption threshold for N2
;at_o  = absorption threshold for O
;it_n2 = ionization threshold for N2
;it_o2 = ionization threshold for O2
;it_o  = ionization threshold for O
;modflux = solar flux in model bins (new bins)
;wave1 = short wavelength of bin
;wave2 = long wavelength of bin
;modsc1= Hinteregger model scale factors for each wavelength bin class 1
;modsc2=              "                   "             "        class 2

at_n2= 986.30 ; Angstroms
at_o2= 2000.00 ; max bins we will be having
at_o = 913.00
it_n2= 798.00
it_o2= 1025.72 ; keep this at Ly Beta wvln as long as Fennally compilation used
it_o = at_o

;IF mode eq 1 THEN BEGIN


;getbins,wave1,wave2
lmax=46
restore, '/home/srimoyee/Desktop/nrl_files/sav_files/newspectra.sav';srimoyee adding the save files for the NRL new bins 
wave1=wave_gcm1
wave2=wave_gcm2

n_bins=n_elements(wave1)
lmax=n_elements(wave1)
wave=(wave1+wave2)/2.0
index=sort(wave)
wave1=wave1(index)
wave2=wave2(index)
wave=wave(index)
iflines=(wave1 eq wave2)  ; 1 where bins are lines, 0 for continuum bins


; get the Hinteregger model solar flux
openr,1,'sc21refw.dat'
n_hint=1661
refwvln=fltarr(n_hint)
refflux=fltarr(n_hint)

format='$(f8.2,f8.1,a11,i1,f10.6)'
label='a stinking string'
class=intarr(n_hint)
r_correlation=fltarr(n_hint)

for i=0,n_hint-1 do begin
	readf,1,format,a,b,label,c,d
	refwvln(i)=a
	refflux(i)=b/1000.  ; to be in units of 1e9 photons/cm2/s
	class(i)=c
	r_correlation(i)=d 	
endfor
close,1
h_wvln=refwvln
f_flux=refflux

; or use high solar activity spectrum
ff=euv81(200.,200.,1)
refwvln=reform(ff(0,*))
refflux=reform(ff(1,*))

IF mode EQ 2 THEN BEGIN

; mode 2 is just a more general version of mode 1, any file can be used
; a small amout of special coding may need to be done
; wave1, wave2 are beginning and end of bins
; refflux is in units of 1e9 ph per cm2 per s

file='nrleuv_bastille_flare.dat'

ref=read_dat(file)
ind=where(ref(0,*) LE 18.)

refwvln=[reform(ref(0,ind)),refwvln]
print, refwvln
refflux=[reform(ref(1,ind))/1e9,refflux]
n_ref=n_elements(refwvln)

;wave1=ref(0,*)-.25
;wave2=ref(0,*)+.25
;wave1=refwvln
;wave2=refwvln

n_hint=n_elements(refwvln)
;class=fltarr(n_hint)-1.
;r_correlation=fltarr(n_hint)+1.
;wave=(wave1+wave2)/2.
;n_bins=n_elements(wave)
;lmax=n_bins

;iflines = (wave1 EQ wave2)

ENDIF

lines=wave1*iflines
ifreflines=fltarr(n_hint)
for i=0,n_hint-1 do for j=0,lmax-1 do begin
	if ((lines(j) eq refwvln(i)) and (iflines(j) eq 1)) then $
	ifreflines(i)=1
endfor

;36.098 solar flux
f107 = 168.9
f107a= 132.6
ff=euv81(f107,f107a,1,www,fff)
refflux_098 = ff(1,*)
refwvln_098 = ff(0,*)

;36.107 solar flux
f107 = 121.0
f107a= 93.1
ff=euv81(f107,f107a,1,www,fff)
refflux_107 = ff(1,*)
refwvln_107 = ff(0,*)

;36.124 solar flux
f107 = 85.6
f107a= 84.6
ff=euv81(f107,f107a,1,www,fff)
refflux_124 = ff(1,*)
refwvln_124 = ff(0,*)
;stop

; get cross section files
loc=''
read_table,loc+'photon2.tab',2,808,15,conwayn2_file
read_table,loc+'photoo2.tab',2,808,14,conwayo2_file
read_table,loc+'photoo.tab' ,2,645, 8,conwayo_file
read_table,loc+'phfenn.tab' ,2,1944,10,fenn
read_table,loc+'phfenno.tab',2,1944,8,fenno

; need to extend conway files out to all Hinteregger wavelengths (0s at end)

conwayn2=fltarr(15,n_ref)
	conwayn2(0,*)=refwvln
	n=n_elements(conwayn2_file(0,*))
	for i=1,14 do conwayn2(i,*)=$
          exp(interpol(alog(conwayn2_file(i,*)),h_wvln(0:807),refwvln))
        iin2=where(refwvln GT conwayn2_file(0,n-1))
        conwayn2(*,iin2)=0.0

conwayo2=fltarr(14,n_ref)
        conwayo2(0,*)=refwvln
        n=n_elements(conwayo2_file(0,*))
        for i=1,13 do conwayo2(i,*)=$
          exp(interpol(alog(conwayo2_file(i,*)),h_wvln(0:807),refwvln))
        iio2=where(refwvln GT conwayo2_file(0,n-1))
        conwayo2(*,iio2)=0.0

conwayo =fltarr(8,n_ref)
        conwayo(0,*)=refwvln
        n=n_elements(conwayo_file(0,*))
        for i=1,7 do conwayo (i,*)=$
          exp(interpol(alog(conwayo_file(i,*)),h_wvln(0:644),refwvln))
        iio=where(refwvln GT conwayo_file(0,n-1))
        conwayo(*,iio)=0.0

iin2=where(finite(conwayn2) eq 0) & if (iin2(0) ne -1) then conwayn2(iin2)=0.0
iio2=where(finite(conwayo2) eq 0) & if (iio2(0) ne -1) then conwayo2(iio2)=0.0
iio =where(finite(conwayo ) eq 0) & if (iio(0)  ne -1) then conwayo(iio)  =0.0

; need to have fennally et al cross sections on same grid as Hinteregger
; model spectrum, Conway stuff already is
; will make sure long wavelengths are zeroed properly

; Add cross sections for wavelengths below 18 Angstroms
; Stan's old values 
;extra_wvln=[.75,1.5,3.,6.,14.]
;extra_o=[2e-23,2e-22,2e-21,1.2e-20,7e-19]*1e18
;extra_o2=[4e-23,4e-22,4e-21,2.4e-20,1.4e-19]*1e18
;extra_n2=[3e-23,3e-22,3e-21,1.5e-20,9e-20]*1e18
; Henke values

extra_ener=[1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,22.,24.,26.,28.]*1000.

extra_n2=[78956.195,26105.123,11492.653,3467.374,1446.114,726.487,412.348,255.309,168.8,117.527,85.348,49.768,32.239,$

  22.663,17.007,13.457,11.115,9.505,8.360,7.519]*1e-28*1e4/1e-18



extra_o=[122161.625,41708.105,18762.01,5823.478,2473.609,1258.528,720.725,448.952,297.923,207.758,150.824,87.428,55.945,$

  38.643,28.379,21.914,17.641,14.7,12.606,11.071]*1e-28*1e4/1e-18



extra_o2=extra_o*2.

extra_wvln=12397./extra_ener

;extra_wvln=[1.25,1.44,1.54,1.66,1.79,2.10,2.29,2.50,2.75, 2.78,3.36,4.16,$
;            4.73,5.41,5.73,6.07,7.13,8.34,9.89,10.44,11.91,12.25,13.33,$
;            14.56,16.0]
;extra_o=[5.75,8.81,11.,13.9,17.6,28.9,37.6,49.4,65.7,67.7,120.,224.,325.,$
;         476.,560.,661.,1030.,1600.,2530.,2920.,4150.,4480.,5600.,6960.,$
;         8880.]*26.56*1e-24*1e18
;extra_o2=extra_o*2.
;extra_n2=[3.55,5.47,6.85,8.65,11.,18.2,23.7,31.3,41.7,43.,76.9,146.,214.,$
;          318.,376.,446.,708.,1110.,1790.,2070.,2970.,3210.,4020.,5070.,$
;          6500.]*2.*23.26*1e-24*1e18
            
fennsigin2=exp(interpol(alog([extra_n2,reform(fenn(4,*) )]),[extra_wvln,reform(fenn(0,*) )],refwvln))
fennsigio2=exp(interpol(alog([extra_o2,reform(fenn(8,*) )]),[extra_wvln,reform(fenn(0,*) )],refwvln))
fennsigio =exp(interpol(alog([extra_o ,reform(fenno(7,*))]),[extra_wvln,reform(fenno(0,*))],refwvln))
fennsigan2=exp(interpol(alog([extra_n2,reform(fenn(1,*) )]),[extra_wvln,reform(fenn(0,*) )],refwvln))
fennsigao2=exp(interpol(alog([extra_o2,reform(fenn(5,*) )]),[extra_wvln,reform(fenn(0,*) )],refwvln))
fennsigao =exp(interpol(alog([extra_o, reform(fenno(7,*))]),[extra_wvln,reform(fenno(0,*))],refwvln))

iin2=where(finite(fennsigin2) eq 0) & if (iin2(0) ne -1) then fennsigin2(iin2)=0.0
iin2=where(finite(fennsigio2) eq 0) & if (iin2(0) ne -1) then fennsigio2(iin2)=0.0
iin2=where(finite(fennsigio)  eq 0) & if (iin2(0) ne -1) then fennsigio(iin2) =0.0
iin2=where(finite(fennsigan2) eq 0) & if (iin2(0) ne -1) then fennsigan2(iin2)=0.0
iin2=where(finite(fennsigao2) eq 0) & if (iin2(0) ne -1) then fennsigao2(iin2)=0.0
iin2=where(finite(fennsigao)  eq 0) & if (iin2(0) ne -1) then fennsigao(iin2) =0.0

indx=where(refwvln gt at_n2) & fennsigan2(indx)=0.0
indx=where(refwvln gt 1027.) & fennsigao2(indx)=0.0 ; Fennally goes 1027A only
                               fennsigao2(indx)=conwayo2(1,indx)
indx=where(refwvln gt at_o ) & fennsigao (indx)=0.0
indx=where(refwvln gt it_n2) & fennsigin2(indx)=0.0 
indx=where(refwvln gt it_o2) & fennsigio2(indx)=0.0
indx=where(refwvln ge it_o ) & fennsigio (indx)=0.0

; Fennallly file also has Ly beta cross section = 0 for ionization of O2,
; this is clearly wrong so will be fixed.
	indx=where(refwvln eq it_o2)
	IF indx(0) NE -1 THEN fennsigio2(indx(0))=1.0
	
; write out fennally cross sections in Hinteregger bins

openw,1,'fennhbins.tab'
printf,1,'Fennally et al cross sections at Hinteregger wavelengths.'
printf,1,'  ','   WVLN   ','    ','  O Ion   ','   ','  O2 Ion  ','   ',$
         '  O2 Abs  ','   ','  N2 Ion  ','   ','  N2 Abs  '   
format='$(1x,6(f10.4,3x))'

indx=where(refwvln le 1027.0)
n=n_elements(indx)
for i=0,n-1 do begin
	printf,1,format,refwvln(i),fennsigio(i),fennsigio2(i),$
                 fennsigao2(i),fennsigin2(i),fennsigan2(i)
endfor
close,1

; Calculate flux weighted cross sections from Fennally compilation.
; also in this loop get solar flux into model bins

modsigin2=fltarr(n_bins)
modsigio2=fltarr(n_bins)
modsigio =fltarr(n_bins)
modsigan2=fltarr(n_bins)
modsigao2=fltarr(n_bins)
modsigao =fltarr(n_bins)
modflux  =fltarr(n_bins)
modflux_098 = fltarr(n_bins)
modflux_107 = fltarr(n_bins)
modflux_124 = fltarr(n_bins)
modsc1=fltarr(n_bins)
modsc2=fltarr(n_bins)

for ibin=0,n_bins-1 do begin

        ; the following if statement should handle lines or continua, it does
        ; assume though that the line is one of the Hinteregger lines

	if iflines(ibin) eq 1 then $
	ind=where((refwvln ge wave1(ibin)) and (refwvln le wave2(ibin))) else $
	begin
	 ind=where((refwvln ge wave1(ibin)) and (refwvln lt wave2(ibin)))


	 if (ind(0) ne -1) then begin
	  indcont=(where(ifreflines(ind) eq 0))
	  oldind=ind
	  if indcont(0) ne -1 then begin 
		indcont=ind(indcont)
		ind=indcont
	  endif
;	  print,wave1(ibin),wave2(ibin),n_elements(ind),n_elements(oldind)
	 endif
	endelse

	if (ind(0) ne -1) then begin

	totflux=total(refflux(ind))
	modflux(ibin)=totflux
	modflux_098(ibin)=total(refflux_098(ind))
	modflux_107(ibin)=total(refflux_107(ind))
	modflux_124(ibin)=total(refflux_124(ind))
	
	if totflux eq 0.0 then totflux=1.0
		
	modsigin2(ibin)=total(fennsigin2(ind)*refflux(ind))/totflux
	modsigio2(ibin)=total(fennsigio2(ind)*refflux(ind))/totflux
	modsigio (ibin)=total(fennsigio (ind)*refflux(ind))/totflux
	modsigan2(ibin)=total(fennsigan2(ind)*refflux(ind))/totflux
	modsigao2(ibin)=total(fennsigao2(ind)*refflux(ind))/totflux
	modsigao (ibin)=total(fennsigao (ind)*refflux(ind))/totflux

	indc1=where(class(ind) eq 1)
	indc2=where(class(ind) eq 2)
	if indc1(0) ne  -1 then begin
              modsc1(ibin)=total(r_correlation(ind(indc1))*refflux(ind(indc1)))$
                          / totflux
	endif
        if indc2(0) ne  -1 then begin
             modsc2(ibin)=total(r_correlation(ind(indc2))*refflux(ind(indc2)))$
                          / totflux
        endif

	endif

endfor

; this is how /glow uses the scale factors
modsc1=modsc1*modflux*1000.
modsc2=modsc2*modflux*1000.

; now need flux weighted branching ratio

probn2=fltarr(6,n_bins)
probo2=fltarr(6,n_bins)
probo =fltarr(6,n_bins)

for ibin=0,n_bins-1 do begin

        ; the following if statement should handle lines or continua, it does
        ; assume though that the line is one of the Hinteregger lines

        if iflines(ibin) eq 1 then $
        ind=where((refwvln ge wave1(ibin)) and (refwvln le wave2(ibin))) else $
        begin
         ind=where((refwvln ge wave1(ibin)) and (refwvln lt wave2(ibin)))
         if (ind(0) ne -1) then begin
          indcont=(where(ifreflines(ind) eq 0))
          oldind=ind
          if indcont(0) ne -1 then begin
                indcont=ind(indcont)
                ind=indcont
          endif
;          print,wave1(ibin),wave2(ibin),n_elements(ind),n_elements(oldind)
         endif
        endelse

	num=n_elements(ind)
	if (ind(0) ne -1) then begin

        totflux=total(refflux(ind))
        if totflux eq 0.0 then totflux=1.0

; N2:  X, A, B, C, F, Dissoc, conway arrays 5,6,7,8,9,3
; Conway says C,F dissociate but Kirby does not agree, for now (to be
; consistant with Stan) I will go with Kirby so I will correct the Conway
; dissociation array (3) by subtracting off values for C,F.  This should
; be a small correction anyway.

; Im using fragmentation cross section to get dissociation, so I need to 
; remove N2++ 
; note that K shell branching ratios are thrown in to dissociation but
; this is done manually later in this program.

	Conwayn2(3,ind)=(conwayn2(3,ind)-conwayn2(8,ind)-conwayn2(9,ind)-$
                                         conwayn2(13,ind))>0.0

	tots=fltarr(num)
	
	ars=[5,6,7,8,9,3]
	for i=0,num-1 do begin
		for j=0,n_elements(ars)-1 do $
		tots(i)=tots(i)+conwayn2(ars(j),ind(i))
		if tots(i) eq 0 then tots(i)=1e32
	endfor

	ind_tots=where(tots ne 1e32)
	if ind_tots(0) ne -1 then totflux2=total(refflux(ind(ind_tots))) else $
           totflux2=totflux

	probn2(0,ibin)=$
        total(conwayn2(5,ind)/tots*refflux(ind)/totflux2)
        probn2(1,ibin)=$
        total(conwayn2(6,ind)/tots*refflux(ind)/totflux2)
        probn2(2,ibin)=$
        total(conwayn2(7,ind)/tots*refflux(ind)/totflux2)
        probn2(3,ibin)=$
        total(conwayn2(8,ind)/tots*refflux(ind)/totflux2)
        probn2(4,ibin)=$
        total(conwayn2(9,ind)/tots*refflux(ind)/totflux2)
        probn2(5,ibin)=$
        total(conwayn2(3,ind)/tots*refflux(ind)/totflux2)

; O2:  X, a+A, b, Dissoc, Conway arrays 5,6,7 3

        tots=fltarr(num)

        ars=[5,6,7,3]
        for i=0,num-1 do begin
		for j=0,n_elements(ars)-1 do $
                tots(i)=tots(i)+conwayo2(ars(j),ind(i))
		if tots(i) eq 0 then tots(i)=1e32
	endfor

        ind_tots=where(tots ne 1e32)
        if ind_tots(0) ne -1 then totflux2=total(refflux(ind(ind_tots))) else $
           totflux2=totflux
 
        probo2(0,ibin)=$
        total(conwayo2(5,ind)/tots*refflux(ind)/totflux2)
        probo2(1,ibin)=$
        total(conwayo2(6,ind)/tots*refflux(ind)/totflux2)
        probo2(2,ibin)=$
        total(conwayo2(7,ind)/tots*refflux(ind)/totflux2)
        probo2(3,ibin)=$
        total(conwayo2(3,ind)/tots*refflux(ind)/totflux2)


; O: 4S, 2Do, 2Po, 4Pe, 2Pe

        tots=fltarr(num)

        ars=[2,3,4,5,6]

        for i=0,num-1 do begin
		for j=0,n_elements(ars)-1 do $
                tots(i)=tots(i)+conwayo(ars(j),ind(i))
		if tots(i) eq 0 then tots(i)=1e32
	endfor

        ind_tots=where(tots ne 1e32)
        if ind_tots(0) ne -1 then totflux2=total(refflux(ind(ind_tots))) else $
           totflux2=totflux
 
        probo(0,ibin)=$
        total(conwayo(2,ind)/tots*refflux(ind)/totflux2)
        probo(1,ibin)=$
        total(conwayo(3,ind)/tots*refflux(ind)/totflux2)
        probo(2,ibin)=$
        total(conwayo(4,ind)/tots*refflux(ind)/totflux2)
        probo(3,ibin)=$
        total(conwayo(5,ind)/tots*refflux(ind)/totflux2)
        probo(4,ibin)=$
        total(conwayo(6,ind)/tots*refflux(ind)/totflux2)
	endif

endfor

; Correct for places where Hinteregger spectrum has no data.

ind1=where((modsigin2 ne 0.0) and (wave2 ge 50.0))
 
indin2=where((modsigin2 eq 0) and (wave2 ge  50.0) and ( wave2 le it_n2))
indio2=where((modsigio2 eq 0) and (wave2 ge  50.0) and ( wave2 le it_o2))
indio =where((modsigio  eq 0) and (wave2 ge  50.0) and ( wave2 le it_o ))
indan2=where((modsigan2 eq 0) and (wave2 ge  50.0) and ( wave2 le at_n2))
indao2=where((modsigao2 eq 0) and (wave2 ge  50.0) and ( wave2 le at_o2))
indao =where((modsigao  eq 0) and (wave2 ge  50.0) and ( wave2 le at_o ))

if indio2(0) ne -1 then begin
print,'Correcting for regions of no solar irradiance:'

wv1=wave(ind1)

wvin2=wave(indin2)
wvio2=wave(indio2)
wvio =wave(indio )
wvan2=wave(indan2)
wvao2=wave(indao2)
wvao =wave(indao )


modsigin2(indin2)=exp(interpol(alog(modsigin2(ind1)),wv1,wvin2))
modsigio2(indio2)=exp(interpol(alog(modsigio2(ind1)),wv1,wvio2))
modsigio (indio )=exp(interpol(alog(modsigio (ind1)),wv1,wvio ))
modsigan2(indan2)=exp(interpol(alog(modsigan2(ind1)),wv1,wvan2))
modsigao2(indao2)=exp(interpol(alog(modsigao2(ind1)),wv1,wvao2))
modsigao (indao )=exp(interpol(alog(modsigao (ind1)),wv1,wvao ))
for i=0,5 do begin
	probn2(i,indin2)=exp(interpol(alog(probn2(i,ind1)),wv1,wvin2))
	probo2(i,indio2)=exp(interpol(alog(probo2(i,ind1)),wv1,wvio2))
	probo (i,indio )=exp(interpol(alog(probo (i,ind1)),wv1,wvio))
endfor
ENDIF

ii=where(finite(modsigin2) eq 0) & if (ii(0) ne -1) then modsigin2(ii)=0.0
ii=where(finite(modsigio2) eq 0) & if (ii(0) ne -1) then modsigio2(ii)=0.0
ii=where(finite(modsigio ) eq 0) & if (ii(0) ne -1) then modsigio(ii) =0.0
ii=where(finite(modsigan2) eq 0) & if (ii(0) ne -1) then modsigan2(ii)=0.0
ii=where(finite(modsigao2) eq 0) & if (ii(0) ne -1) then modsigao2(ii)=0.0
ii=where(finite(modsigao ) eq 0) & if (ii(0) ne -1) then modsigao(ii) =0.0
ii=where(finite(probn2) eq 0) & if (ii(0) ne -1) then probn2(ii)=0.0
ii=where(finite(probo2) eq 0) & if (ii(0) ne -1) then probo2(ii)=0.0
ii=where(finite(probo) eq 0) & if (ii(0) ne -1) then probo(ii)=0.0
;stop
; Correct longwavelength endpoints ( a double check).

indx=where(wave1 gt at_n2) & modsigan2(indx)=0.0
;indx=where(wave1 gt at_o2) & modsigao2(indx)=0.0
indx=where(wave1 gt at_o ) & modsigao (indx)=0.0
indx=where(wave1 gt it_n2) & modsigin2(indx)=0.0 & probn2(*,indx)=0.0
indx=where(wave1 gt it_o2) & modsigio2(indx)=0.0 & probo2(*,indx)=0.0
indx=where(wave1 gt it_o ) & modsigio (indx)=0.0 & probo (*,indx)=0.0

IF 2 EQ 1 THEN BEGIN ; should not need this anymore, but will keep in case ;)

; We don't expect the shortest wavelength bins to ever change so they will
; be hardwired in, they are not in the Conway, Kirby or Fennally compilations.
; *****the .5 to 1A bin needs to be checked, Im making assumptions

modsigin2(0:4)=[.3e-22,.3e-21,.3e-20,.15e-19,.90e-19]*1e18
modsigio2(0:4)=[.4e-22,.4e-21,.4e-20,.24e-19,.14e-18]*1e18
modsigio (0:4)=[.2e-22,.2e-21,.2e-20,.12e-19,.70e-18]*1e18
modsigan2(0:4)=modsigin2(0:4)
modsigao2(0:4)=modsigio2(0:4)
modsigao (0:4)=modsigio (0:4)
probn2(0,0:4)  =probn2(0,5)
probn2(1,0:4)  =probn2(1,5)
probn2(5,0:4)  =probn2(5,5)
probo2(3,0:4)  =1.0
probo (0,0:4)  =probo(0,5)
probo (1,0:4)  =probo(1,5)
probo (2,0:4)  =probo(2,5)
probo (3,0:4)  =probo(3,5)
probo (4,0:4)  =probo(4,5)

ENDIF

; Correct solar flux by extracting lines from bins.
; this should not be needed since main loop ind variable take care of this
; SMB 1-14-95
;or ibin=0,n_bins-1 do begin
	
;if wave1(ibin) eq wave2(ibin) then begin
;	for jbin=0,n_bins-1 do begin
;		if ((wave1(ibin)  ge wave1(jbin)) and $
;                           (wave2(ibin)  le wave2(jbin)) and $
;                           (jbin ne ibin)) then begin        $
;                         modflux(jbin)=modflux(jbin)-modflux(ibin)
;		  modflux_098(jbin)=modflux_098(jbin)-modflux_098(ibin)
;		  modflux_107(jbin)=modflux_107(jbin)-modflux_107(ibin)
;		  modflux_124(jbin)=modflux_124(jbin)-modflux_124(ibin)
;		endif
;	endfor
;endif
;ndfor

; Solar flux from rocket flights is obtained if we are using 131 (1nm) bins.
; turned off for now... SMB 20190517
;if lmax eq 131 and 2 eq 1 then begin;
;
;indxuv1=where((wave1 ge 18.)  and (wave2 le 60. ))
;indxuv2=where((wave1 ge 60.)  and (wave2 le 170.))
;indxuv3=where((wave1 ge 170.) and (wave2 le 300.))
;
;	modflux_107(indxuv1)=modflux_107(indxuv1)*1.35
;	modflux_107(indxuv3)=modflux_107(indxuv3)*2.06
;
;	modflux_124(indxuv1)=modflux_124(indxuv1)*2.23
;	modflux_124(indxuv2)=modflux_124(indxuv2)*2.25
;	modflux_124(indxuv3)=modflux_124(indxuv3)*1.57
;
;	rocket_sflux,wv1r,wv2r,sflux_098,sflux_107,sflux_124,error_098,$
;                               error_107,error_124
;	if total(wave1-wv1r) ne 0.0 then begin  ; make sure wvln intervals are=
;		print,'Bin mismatch with rocket solar flux.'
;		print,'Stopping!'
;		stop
;	endif
;
;	moderror_098=fltarr(n_elements(modflux_098))
;	moderror_107=fltarr(n_elements(modflux_098))
;	moderror_124=fltarr(n_elements(modflux_098))
;
;	; Tom uses -1 to indicate regions of no measurement.
;	
;	ind098=where(sflux_098 ne -1.0)
;	modflux_098(ind098)=sflux_098(ind098)/1e9 ; get in same units
;	moderror_098(ind098)=error_098(ind098)
;	ind107=where(sflux_107 ne -1.0)
;	modflux_107(ind107)=sflux_107(ind107)/1e9  
;	moderror_107(ind107)=error_107(ind107)
;	ind124=where(sflux_124 ne -1.0)
;	modflux_124(ind124)=sflux_124(ind124)/1e9  
;	moderror_124(ind124)=error_124(ind124)
;
;
;print,' '
;indxuv=where((wave2 le 300.) and (wave1 ge 170.))
;print,'Ratio of model solar flux to XUV diode measurement:'
;print,'36.098:',total(sflux_098(indxuv)/1e9)/total(modflux_098(indxuv))	
;print,'36.107:',total(sflux_107(indxuv)/1e9)/total(modflux_107(indxuv))	
;print,'36.124:',total(sflux_124(indxuv)/1e9)/total(modflux_124(indxuv))	
;
;print,' '
;indeuv=where((wave2 le 310.) and (wave1 ge 300.))
;print,'Ratio of model solar flux to EUV measurement:'
;print,'36.098:',total(sflux_098(indeuv)/1e9)/total(modflux_098(indeuv))	
;print,'36.107:',total(sflux_107(indeuv)/1e9)/total(modflux_107(indeuv))	
;print,'36.124:',total(sflux_124(indeuv)/1e9)/total(modflux_124(indeuv))	
;print,'Comparison of model solar flux to EUV measurement (meas,model)'
;print,'36.098:',total(sflux_098(indeuv)/1e9),total(modflux_098(indeuv))	
;print,'36.107:',total(sflux_107(indeuv)/1e9),total(modflux_107(indeuv))	
;print,'36.124:',total(sflux_124(indeuv)/1e9),total(modflux_124(indeuv))	
;
;print,'36.098:',total(sflux_098(indeuv))/1e9,total(modflux_098(indxuv))
;print,'36.107:',total(sflux_107(indeuv))/1e9,total(sflux_107(indxuv))/1e9
;print,'36.124:',total(sflux_124(indeuv))/1e9,total(sflux_124(indxuv))/1e9
;print,' '
;
;temp_bin=where((wave1 ge 60.) and (wave2 le 170.))
;modflux_107(temp_bin)=modflux_107(temp_bin) * 1.7   ; best estimate, not meased
;ojbin1=where((wave1 ge 20.) and (wave2 le 100.))
;ojbin2=where((wave1 ge 50.) and (wave2 le 575.))
;print,'unscaled 36.098 2 - 10 nm:',total(modflux_098(ojbin1)),$
;' 5 - 57.5 nm:',total(modflux_098(ojbin2))
;print,'36.107 2 - 10 nm:',total(modflux_107(ojbin1)),$
;' 5 - 57.5 nm:',total(modflux_107(ojbin2))
;print,'36.124 2 - 10 nm:',total(modflux_124(ojbin1)),$
;' 5 - 57.5 nm:',total(modflux_124(ojbin2))
;modflux_107(temp_bin)=modflux_107(temp_bin) / 1.7   ; change back since not msd
;
;indint=where(wave2 le 310)
;indint2=where(wave2(indint) le 300.)
;tot_098=total(modflux_098(indint))+total(modflux_098(indint(indint2)))
;tot_107=total(modflux_107(indint))
;tot_124=total(modflux_124(indint))
;print,''
;print,'36.098 total soft x-rays = ',tot_098
;print,'36.107 total soft x-rays = ',tot_107
;print,'36.124 total soft x-rays = ',tot_124
;print,' '
;
;; write rocket sflux files:
;openw,1,'refsolspec.36098'
;for i=0,n_elements(modflux_098)-1 do begin
;	printf,1,wave1(i)/10.,wave2(i)/10.,modflux_098(i)*1e9,moderror_098(i)
;endfor
;close,1
;
;openw,1,'refsolspec.36107'
;for i=0,n_elements(modflux_107)-1 do begin
;        printf,1,wave1(i)/10.,wave2(i)/10.,modflux_107(i)*1e9,moderror_107(i)
;endfor
;close,1
;
;openw,1,'refsolspec.36124'
;for i=0,n_elements(modflux_124)-1 do begin
;        printf,1,wave1(i)/10.,wave2(i)/10.,modflux_124(i)*1e9,moderror_124(i)
;endfor
;close,1
;
;	
;endif

IF mode EQ 1 THEN BEGIN

; Load some typical solar min values for <18A wavelengths, see Bailey thesis
; for references
lt18_minflux=[1e1,1e2,2e2,1e4,7e7]/1e9
modflux(0:4)=lt18_minflux
modflux_098(0:4)=lt18_minflux
modflux_107(0:4)=lt18_minflux
modflux_124(0:4)=lt18_minflux

ENDIF

; Now write out the cross sections.

	format="$((1x,f7.2,3x,f7.2,6(3x,f4.2),1x,f9.5,1x,f9.5))"

	openw,1,'newbins_n2.dat'
	printf,1,'Output of newbins.pro'
        printf,1, 'N2 branching ratios and cross sections.'
        printf,1, '  '
        printf,1,' Wavelength Bins (A) ' $
        ,' X  ','   ',' A  ','   ',' B  ','   ',' C  ','   ',' F  ', $
        '   ','Diss','   ','TotIon','   ','TotAbs'
        for j=0,lmax-1 do begin
        printf,1,format, $
        wave1(j),wave2(j),probn2(0:5,j),modsigin2(j),modsigan2(j)
        endfor
        close,1

        openw,1,'newbins_o2.dat'
	printf,1,'Ouput of newbins.pro
        printf,1, 'O2 branching ratios and cross sections.'
        printf,1, '  '
        printf,1,' Wavelength Bins (A) ' $
        ,' X  ','   ','a+A ','   ',' b  ','   ','diss','   ','    ', $
        '   ','    ','   ','TotIon','  ',' TotAbs'
        for j=0,lmax-1 do begin
            IF modsigao2(j) GT 999.999 THEN modsigao2(j) = 0.
        printf,1,format, $
        wave1(j),wave2(j),probo2(0:5,j),modsigio2(j),modsigao2(j)
        endfor
        close,1

        openw,1,'newbins_o.dat'
	printf,1,'Output of newbins.pro'
        printf,1, 'O branching ratios and cross sections.'
        printf,1, '  '
        printf,1,' Wavelength Bins (A) ' $
        ,' 4s ','   ','2Do ','   ','2Po ','   ','4Pe ','   ','2Pe ', $
        '   ','    ','   ','TotIon','  ',' TotAbs'
        for j=0,lmax-1 do begin
        printf,1,format, $
        wave1(j),wave2(j),probo(0:5,j),modsigio(j),modsigao(j)
        endfor
        close,1

; write out solar flux file

	openw,1,'newbins_sflux.dat'
	printf,1,'Output of newbins.pro'
	printf,1,'Solar flux in model bins.
	printf,1,'   '
        printf,1,' Wavelength Bins (A) ','    Solar Flux '
	format='$(1x,f7.2,3x,f7.2,3x,6f17.9)'
	for i=0,lmax-1 do begin
		printf,1,format,wave1(i),wave2(i),modflux(i),$
                                modsc1(i),modsc2(i),$
				modflux_098(i),modflux_107(i),modflux_124(i)
	endfor
	close,1
end












