
import BenchmarkExample: BenchmarkExample
n = 18
# BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = n+1)
# BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = (n+1,2*n+1), mode = 2)
BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = (2*n+1,n+1), coef = 0.96)
# 2,3,4 0.8
# 5 0.85
# 6,7 0.9
# 8 0.91
# 9 0.92
# 10 0.93
# 11,12,13 0.94
# 16 0.955