
import BenchmarkExample: BenchmarkExample
n = 8
# BenchmarkExample.PatchTest.generateMsh("./msh/patchtest_"*string(4278)*".msh", transfinite = n+1)
# BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = (n+1,2*n+1), mode = 2)
BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = (2*n+1,n+1), coef = 0.9)
