
import BenchmarkExample: BenchmarkExample
n = 3
# BenchmarkExample.PatchTest.generateMsh("./msh/patchtest_"*string(4278)*".msh", transfinite = n+1)
BenchmarkExample.PlateWithHole.generateMsh("./msh/PlateWithHole_"*string(n)*".msh", transfinite = n)
