#El siguiente script parte del script carga una lista "sod_cla_Serial" a partir de numeros
#de iones almacenados en un archivo de texto "serial.txt". Podría ser util en otro momento...

#set ch [open "serial.txt" r] #serial numbers of ions
#  if {[gets $ch line] != -1} { 
#    set sod_cla_Serials [split $line] 
#	}

#El siguiente script define los límines de zmin y zmax a aplicar para el borde.
set upperZlimit 89.0 ; #  
set lowerZlimit -10.0 ; #  
set padding 3.0 ;# umbral de distancia para calcular la fuerza sobre un átomo cerca del borde.
set springConstant -10 ;#constante de fuerza para aplicar a las moléculas en el padding. 

# El siguente comando que define cómo se obtienen las coordenadas para "proc calcforces"
#wrapmode input ; # the position corresponding to the input files of the simulation (tal como son...)
wrapmode cell ; # the equivalent position in the unit cell centered on the cellOrigin
#wrapmode patch ; #(default) the position in NAMD's internal patch data structure, requires no extra calculation and is almost the same as cell
#wrapmode nearest ; #the equivalent position nearest to the cellOrigin

#El siguiente comando realiza los cálculos.
proc calcforces {step unique} { 
    #atomos ingorados se actualizan cada 100 ciclos. (cambiar el 100 por otro valor si es necesario)
    if {$step % 10 == 0 } {cleardrops;} 
    #inicia los cálculos por cada atomo...
    while {[nextatom]} { 
      set atomid [getid]
      #if { [lsearch $::sod_cla_Serials $atomid] < 0 } {dropatom; continue;} ;#solo activar si se usa una lista de átomos específicos.
      #obtiene las coordenadas
      foreach {x y z} [getcoord] { break; } 
      #acá se revisa si el átomo entra al cálculo (si está suficientemente cerca) ... no parece necesario 
      if {$z < $::upperZlimit - $::padding && $z > $::lowerZlimit + $::padding} { dropatom; continue;} 
      #crea una variable de fuerzas y una variable fz, todas cero.
      set force {0. 0. 0.} 
      set fz 0.0 
      set enerBC 0.0
      # aplica la fuerza si una molécula se encuentra próximo al borde (según valor de padding)
      if {$z > $::upperZlimit - $::padding} { set fz [expr {$::springConstant * ($z - $::upperZlimit + $::padding) }]; } 
      if {$z > $::upperZlimit - $::padding} { set enerBC [expr {0.5 * $::springConstant * ($z - $::upperZlimit + $::padding) * ($z - $::upperZlimit + $::padding)}]; } 
      if {$z < $::lowerZlimit + $::padding} { set fz [expr {$::springConstant * ($z - $::lowerZlimit - $::padding) }]; } 
      if {$z < $::lowerZlimit + $::padding} { set enerBC [expr {0.5 * $::springConstant * ($z - $::lowerZlimit - $::padding) * ($z - $::lowerZlimit - $::padding) }]; } 
      #este comando asigna el valor de fz a la lista "force" (la lista es fx, fy, fz en posiciones 0, 1, 2)
      lset force 2 $fz 
      #añade las fuerzas calculadas en el borde a las fuerzas del sistema
      addforce $force 
      addenergy $enerBC 
    } 
} 
