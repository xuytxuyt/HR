using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_patch_test.jl")

ndiv = 1
elements, nodes, nodes_p = import_patchtest_mix("msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
nᵖ = 3*length(elements["Ωˢ"])
nᵘ = length(nodes)

E = 1.0;
ν = 0.3;
n = 1
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)
u(x,y) = (1.0+2.0*x+3.0*y)^n
v(x,y) = (4.0+5.0*x+6.0*y)^n
∂u∂x(x,y) = 2.0*n*(1.0+2.0*x+3.0*y)^(n-1)
∂u∂y(x,y) = 3.0*n*(1.0+2.0*x+3.0*y)^(n-1)
∂v∂x(x,y) = 5.0*n*(4.0+5.0*x+6.0*y)^(n-1)
∂v∂y(x,y) = 6.0*n*(4.0+5.0*x+6.0*y)^(n-1)
∂²u∂x²(x,y) = 4.0*n*(n-1)*(1.0+2.0*x+3.0*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6.0*n*(n-1)*(1.0+2.0*x+3.0*y)^abs(n-2)
∂²u∂y²(x,y) = 9.0*n*(n-1)*(1.0+2.0*x+3.0*y)^abs(n-2)
∂²v∂x²(x,y) = 25.0*n*(n-1)*(4.0+5.0*x+6.0*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30.0*n*(n-1)*(4.0+5.0*x+6.0*y)^abs(n-2)
∂²v∂y²(x,y) = 36.0*n*(n-1)*(4.0+5.0*x+6.0*y)^abs(n-2)
σ₁₁(x,y) = Cᵢᵢᵢᵢ*∂u∂x(x,y) + Cᵢᵢⱼⱼ*∂v∂y(x,y)
σ₂₂(x,y) = Cᵢᵢⱼⱼ*∂u∂x(x,y) + Cᵢᵢᵢᵢ*∂v∂y(x,y)
σ₁₂(x,y) = Cᵢⱼᵢⱼ*(∂u∂y(x,y) + ∂v∂x(x,y))
∂σ₁₁∂x(x,y) = Cᵢᵢᵢᵢ*∂²u∂x²(x,y) + Cᵢᵢⱼⱼ*∂²v∂x∂y(x,y)
∂σ₂₂∂y(x,y) = Cᵢᵢⱼⱼ*∂²u∂x∂y(x,y) + Cᵢᵢᵢᵢ*∂²v∂y²(x,y)
∂σ₁₂∂x(x,y) = Cᵢⱼᵢⱼ*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂σ₁₂∂y(x,y) = Cᵢⱼᵢⱼ*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
eval(prescribeForFem)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Ωᵖ"])
set𝝭!(elements["∂Ωˢ"])
set𝝭!(elements["∂Ωᵖ"])
set𝝭!(elements["Γ¹"])
set𝝭!(elements["Γ²"])
set𝝭!(elements["Γ³"])
set𝝭!(elements["Γ⁴"])
set𝝭!(elements["Γᵖ"])

ops = [
    Operator{:∫σᵢⱼσₖₗdΩ}(:E=>E,:ν=>ν),
    Operator{:∫∇σᵢⱼuᵢdΩ}(),
    Operator{:∫σᵢⱼnⱼuᵢdΓ}(),
    Operator{:∫∫vᵢbᵢdxdy}(),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫σᵢⱼnⱼgᵢdΓ}(),
    Operator{:L₂}(:E=>E,:ν=>ν),
    Operator{:H₁}(:E=>E,:ν=>ν),
]

k1 = zeros(3*nᵖ,2*nᵘ)
k2 = zeros(3*nᵖ,2*nᵘ)
k3 = zeros(3*nᵖ,2*nᵘ)
kᵖᵖ = zeros(3*nᵖ,3*nᵖ)
fᵖ = zeros(3*nᵖ)
fᵘ = zeros(2*nᵘ)
ops[1](elements["Ωˢ"],kᵖᵖ)
ops[2](elements["Ωˢ"],elements["Ωᵖ"],k1)
ops[3](elements["∂Ωˢ"],elements["∂Ωᵖ"],k2)
ops[3](elements["Γˢ"],elements["Γᵖ"],k3)
ops[4](elements["Ω"],fᵘ)
# ops[5](elements["Γ¹"],fᵘ)
# ops[5](elements["Γ²"],fᵘ)
# ops[5](elements["Γ³"],fᵘ)
# ops[5](elements["Γ⁴"],fᵘ)
ops[6](elements["Γ¹ˢ"],fᵖ)
ops[6](elements["Γ²ˢ"],fᵖ)
ops[6](elements["Γ³ˢ"],fᵖ)
ops[6](elements["Γ⁴ˢ"],fᵖ)
kᵘᵖ = -k2+k1+k3
k = [ zeros(2*nᵘ,2*nᵘ) kᵘᵖ';kᵘᵖ kᵖᵖ]
f = [fᵘ;fᵖ]
d = k\f
d₁ = d[1:2:2*nᵘ]
d₂ = d[2:2:2*nᵘ]
push!(nodes,:d=>d₁)
# push!(nodes,:d=>d₂)
set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
L₂ = ops[7](elements["Ωᵍ"])
a = log10(L₂)
