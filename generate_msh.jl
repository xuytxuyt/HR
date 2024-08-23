
import BenchmarkExample: BenchmarkExample
n = 2
BenchmarkExample.PatchTest.generateMsh("./msh/PatchTest_"*string(n-1)*".msh", transfinite = n)