
import BenchmarkExample: BenchmarkExample
n = 31
BenchmarkExample.PatchTest.generateMsh("./msh/patchtest_"*string(4278)*".msh", transfinite = n+1)