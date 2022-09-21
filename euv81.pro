;
; Procedure EUV81 is the Hinteregger et al., 1981
; F10.7-association proxy model using the
; SC#21REFW reference spectrum
;
; Input:
;     F107=10.7 cm flux,
;     F81=10.7 cm flux 81-day average.
;     OPTION = 0 (ask user), 1 (model res.), 2 (1 Ang. bin),
;      or 3 (1 nm bin) for resulting flux array
;      or 4 (1 nm bin smoothed for w resolution)
;
; Output:   W=wavelength array (Angstroms), (always model resolution !!!)
;     FR=solar min reference flux; (gigaphotons cm-2 s-1),
;     F=scaled flux (gigaphotons cm-2 s-1),
;     C=contrast ratio,
;     K=chromospheric(1)/coronal(2) flag,
;     ID=line identification.
;
; Returns:  Flux Array being FLTARR(2, N)
;      There is an option to get model resolution,
;      1 Angstrom bins or 1 nm bins
;
; Stan Solomon, 3/92
;
; modified: TNW  3/92   permit binning and simplify interface
;

function euv81, f107, f81, option, w_data, f_data

common euv81_data, iread, fr, c, k, id, w, f
if n_elements(fr) le 1 then iread = 1

if n_params(0) lt 2 then begin
  print, 'Usage: flux = EUV81( f107, f107_81_avg, [opt, w, f, ird, fr,'$
    +'c, kk, id] )'
  return, -1
endif

if n_params(0) lt 3 then option = 0

if (option eq 4) then begin
  if (n_params(0) lt 4) then wsm = 1.0 else wsm = w(0)
endif

; if (n_params(0) lt 10) then iread = 1

; These coefficients reduce to the reference spectrum at solar minimum, defined
; for cycle 21 as F10.7 = 67.6, F10.7 81 day average = 71.5:

; x1=1.0 & y1=0.0138 & z1=0.005 & x2=1.0 & y2=0.59425 & z2=0.3811

; These coefficients are the "best fit" values but give some negative values
; at solar minimum:

 x1=1.31 & y1=0.01106 & z1=0.00492 & x2=-6.618 & y2=0.66159 & z2=0.38319

;   calculate R1 and R2 indices from the F10.7 and <F10.7>
r=fltarr(3)
r(0) = 0.
r(1) = x1 + y1*(f81-71.5) + z1*(f107-f81+3.9)
r(2) = x2 + y2*(f81-71.5) + z2*(f107-f81+3.9)

if iread ne 0 then begin
  ww = 0.0 & ffr = 0.0 & kk = 0 & cc = 0.0
  w=fltarr(1661) & fr=w & f=w & c=w & k=intarr(1661) & id=strarr(1661) & iid=''
  print, 'Reading the EUV81 data file...'
  openr,1,'sc21refw.dat'
  for i=0,1660 do begin
    readf,1,format='(f8.2,f8.1,a11,i1,f9.5)',ww,ffr,iid,kk,cc
    w(i)=ww
    fr(i)=ffr/1000.
    id(i)=iid
    k(i)=kk
    c(i)=cc
  endfor
  close,1
  iread = 0 ; file has been read into common variables
endif
;
;   calculate the flux at the model resolution
;
f = fr + fr * c * (r(k) - 1.)

if option eq 0 then $
  read, 'Binning options: 1 = Model Res., 2 = 1 Ang., 3 = 1 nm ? ', option

w1 = min(w)
w2 = max(w)
nf = n_elements(f)

if (option eq 2) or (option eq 4) then begin
          ; 1 Angstrom bins on 1 Angstrom centers
    wmin = fix(w1)
    wmax = fix(w2)
    if wmax ne w2 then wmax = wmax + 1
    na = wmax - wmin + 1
    fa = fltarr(2, na)
    fa(0,*) = findgen(na) + wmin
    for j = 0, nf-1 do begin
      iw = fix( w(j) + 0.5 ) - wmin
      fa(1,iw) = fa(1,iw) + f(j)
    endfor

    if (option eq 4) then begin
              ; smooth and make into 1 nm bins
    wstart = fix( fa(0,0)/10 + 0.9 )
    wend = fix( fa(0,na-1)/10 )
    na = wend - wstart - 1
    new = fltarr(2, na )
    new(0,*) = findgen(na) + wstart + 0.5
    ns = fix(wsm * 10) + 1
    if ns lt 3 then ns = 3
    fas = smooth( smooth( fa(1,*), ns), ns )   ; triangle smooth
    ks = where( fa(0,*) eq wstart * 10 )
    ks = ks(0)
    for j = 0, na-1 do new(1,j) = total( fas(ks+j*10:ks+j*10+9) )
    fa = new
    endif

endif else if (option eq 3) then begin
          ; 1 nm bins on 0.5 nm centers
    wmin = fix(w1/10.) + 0.5
    wmax = fix(w2/10.) - 0.5
    if wmax ne (w2/10.-0.5) then wmax = wmax + 1
    na = fix(wmax - wmin + 1)
    fa = fltarr(2, na)
    fa(0,*) = findgen(na) + wmin
    for j = 0, nf-1 do begin
      iw = fix( fix( w(j)/10.) + 0.5 - wmin )
      if (iw + wmin - 0.5) eq (w(j)/10.) then begin    ; split between bins
        fa(1,iw-1) = fa(1,iw-1) + f(j)/2
        if iw lt na then fa(1,iw) = fa(1,iw) + f(j)/2
      endif else begin
        fa(1,iw) = fa(1,iw) + f(j)
      endelse
    endfor

endif else begin
          ; no binning - use model resolution
    fa = fltarr(2, nf)
    fa(0,*) = w
    fa(1,*) = f
endelse

w_data = w
f_data = f
return, fa
end
