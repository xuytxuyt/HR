
using ApproxOperator, XLSX, TimerOutputs
# using SparseArrays, Pardiso
using ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫∇𝑛vgds, ∫vbdΩ, L₂, H₁

include("import_patch_test.jl")

ndiv = 16
elements, nodes = import_patchtest_gauss("msh/patchtest_"*string(ndiv)*".msh");
# ps = MKLPardisoSolver()
const to = TimerOutput()

nₚ = length(nodes)
nₑ = length(elements["Ω"])
@timeit to "shape function" begin
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set∇𝝭!(elements["Γ"])
end
n = 1
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
prescribe!(elements["Ωᵍ"],:∂u∂z=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ²"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ³"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γ¹"],:α=>(x,y,z)->1E3)
prescribe!(elements["Γ²"],:α=>(x,y,z)->1E3)
prescribe!(elements["Γ³"],:α=>(x,y,z)->1E3)
prescribe!(elements["Γ⁴"],:α=>(x,y,z)->1E3)

𝑎 = ∫∫∇v∇udxdy=>elements["Ω"]
𝑎ᵅ = ∫∇𝑛vgds=>elements["Γ"]
𝑓 = ∫vbdΩ=>elements["Ω"]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "assembly matrix" begin

𝑎(k)
𝑎ᵅ(k,f)
𝑓(f)
end


# k = sparse([kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(nₚ,nₚ)])
# set_matrixtype!(ps,-2)
# k = get_matrix(ps,k,:N)
d = k\f
# f = [fᵖ;-fᵘ]
# @timeit to "solve" pardiso(ps,d,k,f)

push!(nodes,:d=>d)

# L₂_𝑢 = L₂(elements["Ωᵍ"])
𝐻ₑ, 𝐿₂  = H₁(elements["Ωᵍ"])

println(𝐿₂)
println(𝐻ₑ)
# XLSX.openxlsx("./xlsx/patchtest.xlsx", mode="rw") do xf
# index = 64
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndivu,index)+1
#     Sheet["O"*string(ind)] = 3*nₑ
#     Sheet["P"*string(ind)] = log10(𝐿₂)
#     Sheet["Q"*string(ind)] = log10(𝐻ₑ)
# end

show(to)