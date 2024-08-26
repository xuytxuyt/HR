
import BenchmarkExample: BenchmarkExample
n = 8
BenchmarkExample.PatchTest.generateMsh("./msh/patchtest_"*string(n)*".msh", transfinite = n+1)