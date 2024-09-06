
using ApproxOperator, XLSX, TimerOutputs
using ApproxOperator.Heat: ∫∫qᵢpᵢdxdy, ∫pᵢnᵢuds, ∫∫∇𝒑udxdy, ∫pᵢnᵢgⱼds, ∫vbdΩ, ∫vgdΓ, L₂, L₂𝒑, H₁

include("import_patch_test.jl")

# nₚ = 49
ndivu = 64
ndiv = 64
# elements, nodes = import_patchtest_mix("msh/patchtest_u_"*string(nₚ)*".msh","./msh/patchtest_"*string(ndiv)*".msh");
elements, nodes = import_patchtest_mix("msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndivu)*".msh");

const to = TimerOutput()

nₛ = 3
nₚ = length(nodes)
nₑ = length(elements["Ω"])

n = 5
# u(x,y) = (x+y)^n
# ∂u∂x(x,y) = n*(x+y)^abs(n-1)
# ∂u∂y(x,y) = n*(x+y)^abs(n-1)
# ∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
u(x,y) = (1+2*x+3*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)
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

@timeit to "assembly matrix" begin

kᵖᵖ = zeros(2*nₛ*nₑ,2*nₛ*nₑ)
fᵖ = zeros(2*nₛ*nₑ)
kᵖᵘ = zeros(2*nₛ*nₑ,nₚ)
fᵘ = zeros(nₚ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)
end
@timeit to "solve" begin

d = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(nₚ,nₚ)]\[fᵖ;-fᵘ]
end

𝑢 = d[2*nₛ*nₑ+1:end]
push!(nodes,:d=>𝑢)

# L₂_𝑢 = L₂(elements["Ωᵍ"])
𝐻ₑ, 𝐿₂  = H₁(elements["Ωᵍ"])

# XLSX.openxlsx("./xlsx/patchtest.xlsx", mode="rw") do xf
# index = 64
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndivu,index)+1
#     Sheet["O"*string(ind)] = 3*nₑ
#     Sheet["P"*string(ind)] = log10(𝐿₂)
#     Sheet["Q"*string(ind)] = log10(𝐻ₑ)
# end

show(to)