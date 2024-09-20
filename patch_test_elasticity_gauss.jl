
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress

include("import_patch_test.jl")

ndiv = 16
elements, nodes = import_patchtest_gauss("msh/patchtest_"*string(ndiv)*".msh")

nₚ = length(nodes)
nₑ = length(elements["Ω"])

E = 1.0
ν = 0.3

set∇𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Γ"])

n = 1
u(x,y) = (1+2*x+3*y)^n
v(x,y) = (4+5*x+6*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1-ν^2)*(ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1-ν^2)*(ν*ε₁₁(x,y) + ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

∂σ₁₁∂x(x,y) = E/(1-ν^2)*(∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1-ν^2)*(∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂x(x,y) + ∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂y(x,y) + ∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)

prescribe!(elements["Ω"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ¹"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Γ¹"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ²"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Γ²"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ³"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Γ³"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ⁴"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Γ⁴"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γ¹"],:α=>(x,y,z)->E*1e3, index=:𝑔)
prescribe!(elements["Γ²"],:α=>(x,y,z)->E*1e3, index=:𝑔)
prescribe!(elements["Γ³"],:α=>(x,y,z)->E*1e3, index=:𝑔)
prescribe!(elements["Γ⁴"],:α=>(x,y,z)->E*1e3, index=:𝑔)
prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ²"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ³"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ³"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ³"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ³"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

𝑎 = ∫∫εᵢⱼσᵢⱼdxdy=>elements["Ω"]
𝑎ᵅ = ∫σᵢⱼnⱼgᵢds=>elements["Γ"]
𝑓 = ∫∫vᵢbᵢdxdy=>elements["Ω"]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

𝑎(k)
𝑎ᵅ(k,f)
𝑓(f)

d = k\f

d₁ = d[1:2:end]
d₂ = d[2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

# 𝐿₂ = L₂(elements["Ωᵍ"])
 𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍ"])