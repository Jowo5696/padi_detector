set terminal epslatex color
set notitle
set xlabel 'position/mm'
set ylabel 'events/\#'
set grid
set key box top left width -4 # 'samplen x' sets how much space the symbol takes

set style line 1 lt rgb "#1f78b4" pt 13 ps .5 # blue
set style line 2 lt rgb "#33a02c" pt 13 ps .5 # green
set style line 3 lt rgb "#e31a1c" pt 13 ps .5 # red
set style line 4 lt rgb "#fdbf6f" pt 13 ps .5 # yellow
set style line 5 lt rgb "#ff7f00" pt 13 ps .5 # orange
set style line 6 lt rgb "#6a3d9a" pt 13 ps .5 # purple
set style line 7 lt rgb "#a6cee3" pt 13 ps .5 # pastel blue
set style line 8 lt rgb "#b2df8a" pt 13 ps .5 # pastel green
set style line 9 lt rgb "#fb9a99" pt 13 ps .5 # pastel red
set style line 10 lt rgb "#cab2d6" pt 13 ps .5 # pastel purple

# in mm
array d[12]
d[1] = 400
d[2] = 400
d[3] = 400
d[4] = 150 #
d[5] = 200 #
d[6] = 400
d[7] = 447 #
d[8] = 400
d[9] = 400
d[10] = 400
d[11] = 400
d[12] = 400

# in cm
array th[12]
th[1] = 18.65
th[2] = 21
th[3] = 2.96
th[4] = 5
th[5] = 5
th[6] = 5
th[7] = 5
th[8] = .7
th[9] = 9.93
th[10] = 1.48
th[11] = 27.39
th[12] = 4.36

# in cm
x0_cu = 1.436
x0_al = 8.897
array x0[12]
x0[1] = x0_al
x0[2] = x0_al
x0[3] = x0_cu 
x0[4] = x0_al 
x0[5] = x0_al 
x0[6] = x0_al 
x0[7] = x0_al 
x0[8] = x0_cu 
x0[9] = x0_al 
x0[10] = x0_cu 
x0[11] = x0_al 
x0[12] = x0_cu 

# in rad
array theta[12]
# in mm
array c[12]

do for [i = 1:12] {
  set terminal epslatex color
  set output 'output_'.i.'.tex'
  set notitle
  set xlabel 'position/mm'
  set ylabel 'events/\#'
  set grid
  set key box top left width -4 # 'samplen x' sets how much space the symbol takes

  files = '../../src/elsa/run_'.(i-1).'.txt'

  f(x) = a * exp( -(x - b)**n / S**n ) + g
  n = 2 
  a = 8000 # amplitude
  b = 50 # mean
  S = 9 # sigma
  g = 1000 # offset

  # might need to adjust the taken data
  #every ::30::70
  fit f(x) files u 1:2:3 every ::1::97 yerrors via a,b,S,g
  c[i] = S
  chi2 = (FIT_STDFIT*FIT_STDFIT)
  set label sprintf("$\\chi^2 = %.5f$", chi2) at 10,700
  set label sprintf("$\\sigma = %.5f$", c[i]) at 70,700

  theta[i] = atan(c[i]/d[i])

  plot files u 1:2 with lines,\
    f(x)
  reset
}

set terminal epslatex color
set output 'output_correlation.tex'
set notitle
set xlabel '$\sqrt{x/x_0}(1+0.038\ln(x/x_0))$'
set ylabel '$\theta_0=\arctan(\sigma/d)$'
set grid
set key box top left width -4 # 'samplen x' sets how much space the symbol takes

# val needs int not dynamic values
#$val << EOD
#th[1] theta[1]
#th[2] theta[2]
#th[6] theta[6]
#th[9] theta[9]
#th[11] theta[11]
#EOD

set print "correlation.dat"
print sprintf("%.2f %.5f", sqrt(th[1]/x0[1]) * (1 + .038 * log(th[1] / x0[1])), theta[1])
print sprintf("%.2f %.5f", sqrt(th[2]/x0[2]) * (1 + .038 * log(th[2] / x0[2])), theta[2])
print sprintf("%.2f %.5f", sqrt(th[6]/x0[6]) * (1 + .038 * log(th[6] / x0[6])), theta[6])
print sprintf("%.2f %.5f", sqrt(th[9]/x0[9]) * (1 + .038 * log(th[9] / x0[9])), theta[9])
print sprintf("%.2f %.5f", sqrt(th[11]/x0[11]) * (1 + .038 * log(th[11] / x0[11])), theta[11])
set print
print sprintf("%.2f %.5f", sqrt(th[1]/x0[1]) * (1 + .038 * log(th[1] / x0[1])), theta[1])
print sprintf("%.2f %.5f", sqrt(th[2]/x0[2]) * (1 + .038 * log(th[2] / x0[2])), theta[2])
print sprintf("%.2f %.5f", sqrt(th[6]/x0[6]) * (1 + .038 * log(th[6] / x0[6])), theta[6])
print sprintf("%.2f %.5f", sqrt(th[9]/x0[9]) * (1 + .038 * log(th[9] / x0[9])), theta[9])
print sprintf("%.2f %.5f", sqrt(th[11]/x0[11]) * (1 + .038 * log(th[11] / x0[11])), theta[11])

#h(x) = m * x + b
#fit h(x) "correlation.dat" u 1:2:(1) yerrors via m,b
plot "correlation.dat" ls 1
#  h(x)

# NaN with points / lines title '' ls 1 # fake legend

# pt 0 pixel
# pt 1 plus
# pt 13 dot
# lw linewidth
# ps pointsize
# yerrorbars uses last column
# lc rgb #d95f02 (orange) #1b9e77 (green) #7570b3 (blue)
