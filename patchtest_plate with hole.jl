
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫σᵢⱼσₖₗdxdy, ∫∫∇σᵢⱼuᵢdxdy, ∫σᵢⱼnⱼuᵢds, ∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress

include("import_plate_with_hole.jl")

ndivs = 2
ndiv = 2
# elements, nodes = import_patchtest_mix("msh/patchtest_u_"*string(nₚ)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
elements, nodes = import_plate_with_hole_mix("msh/PlateWithHole_"*string(ndivs)*".msh","./msh/PlateWithHole_"*string(ndiv)*".msh");

nₛ = 3
nₚ = length(nodes)
nₑ = length(elements["Ω"])

E = 1.0
ν = 0.3
T = 1000.0
a = 1.0

set𝝭!(elements["Ω"])
set𝝭!(elements["∂Ω"])
set∇𝝭!(elements["Ωᵍ"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γᵍ"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["∂Ωˢ"])

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

prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
# prescribe!(elements["Γᵗ"],:g₁=>(x,y,z)->u(x,y))
# prescribe!(elements["Γᵗ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

𝑎 = ∫∫σᵢⱼσₖₗdxdy=>elements["Ωˢ"]
𝑏 = [
    ∫σᵢⱼnⱼuᵢds=>(elements["∂Ωˢ"],elements["∂Ω"]),
    ∫∫∇σᵢⱼuᵢdxdy=>(elements["Ωˢ"],elements["Ω"]),
]
𝑏ᵅ = ∫σᵢⱼnⱼgᵢds=>(elements["Γˢ"],elements["Γ"])
𝑓 = [
    ∫∫vᵢbᵢdxdy=>elements["Ω"],
    ∫vᵢtᵢds=>elements["Γᵗ"],
]

kᵖᵖ = zeros(3*nₛ*nₑ,3*nₛ*nₑ)
fᵖ = zeros(3*nₛ*nₑ)
kᵖᵘ = zeros(3*nₛ*nₑ,2*nₚ)
fᵘ = zeros(2*nₚ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)

d = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(2*nₚ,2*nₚ)]\[fᵖ;-fᵘ]
d₁ = d[3*nₛ*nₑ+1:2:end]
d₂ = d[3*nₛ*nₑ+2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)

# # 𝐿₂ = L₂(elements["Ωᵍ"])
𝐿₂, 𝐻ₑ = Hₑ_PlaneStress(elements["Ωᵍ"])
