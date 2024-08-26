
using ApproxOperator
using ApproxOperator.Heat: ∫∫qᵢpᵢdxdy, ∫pᵢnᵢuds, ∫∫∇𝒑udxdy, ∫pᵢnᵢgⱼds, ∫vbdΩ, ∫vgdΓ, L₂, L₂𝒑

include("import_patch_test.jl")

# nₚ = 49
ndiv = 2
# elements, nodes = import_patchtest_mix("msh/patchtest_u_"*string(nₚ)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
elements, nodes = import_patchtest_mix("msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh");

nₛ = 3
nₚ = length(nodes)
nₑ = length(elements["Ω"])

n = 2
u(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γ¹"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))

𝑎 = ∫∫qᵢpᵢdxdy=>elements["Ωˢ"]
𝑏 = [
    ∫pᵢnᵢuds=>(elements["∂Ωˢ"],elements["∂Ω"]),
    ∫∫∇𝒑udxdy=>(elements["Ωˢ"],elements["Ω"]),
]
𝑏ᵅ = ∫pᵢnᵢgⱼds=>(elements["Γˢ"],elements["Γ"])
𝑓 = ∫vbdΩ=>elements["Ω"]

kᵖᵖ = zeros(2*nₛ*nₑ,2*nₛ*nₑ)
fᵖ = zeros(2*nₛ*nₑ)
kᵖᵘ = zeros(2*nₛ*nₑ,nₚ)
fᵘ = zeros(nₚ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)

d = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(nₚ,nₚ)]\[fᵖ;-fᵘ]

𝑢 = d[2*nₛ*nₑ+1:end]
push!(nodes,:d=>𝑢)

L₂_𝑢 = L₂(elements["Ωᵍ"])
