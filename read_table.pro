pro read_table,file,n_string,n_lines,n_colums,data

; read a text file and return the table of data

data=fltarr(n_colums,n_lines)

openr,1,file

s='a string'
for i=1,n_string do readf,1,s
readf,1,data
close,1

return
end