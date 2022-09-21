; procedure to get the solar flux wavelength bins used by glow
; SMB 8/94
; SMB 8.2002 - made work for UAF/GI setup, added some comments

pro getbins,wave1,wave2

; input file for wavelength bins  is glow_sfbins.txt
; file should have 1 title line
; line 2 should start with: number of bins, scale factor to conver to
; A

loc='F:\OneDrive\Documents\ACE Software\Solar\bins\' 
loc='D:\bailey_onedrive\OneDrive\Documents\ACE Software\Solar\bins\'
;openr,1,loc+'glow_sfbins.txt'
openr,1,'glow_sfbins.TXT_OLD'
s='a string'
readf,1,s
readf,1,nbins,scale

wave1=fltarr(nbins)
wave2=fltarr(nbins)

for i=0,nbins-1 do begin
	readf,1,b,c
	wave1(i)=b
	wave2(i)=c
endfor
close,1

wave1=wave1*scale
wave2=wave2*scale

return
end
