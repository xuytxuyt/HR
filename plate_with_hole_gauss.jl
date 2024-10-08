
using ApproxOperator, XLSX, TimerOutputs
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫σᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, Hₑ_PlaneStress

include("import_plate_with_hole.jl")

ndiv = 32
elements, nodes = import_plate_with_hole_gauss("msh/PlateWithHole_"*string(ndiv)*".msh",2*ndiv,0.98);
# ps = MKLPardisoSolver()
const to = TimerOutput()

nₚ = length(nodes)
nₑ = length(elements["Ω"])
@timeit to "shape function" begin
    set∇𝝭!(elements["Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set∇𝝭!(elements["Γᵍ"])
    set𝝭!(elements["Γᵗ"])
end
T = 1000.0
E = 3e6
ν = 0.3
a = 1.0
r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
u(x,y) = T*a*(1+ν)/2/E*( r(x,y)/a*2/(1+ν)*cos(θ(x,y)) + a/r(x,y)*(4/(1+ν)*cos(θ(x,y))+cos(3*θ(x,y))) - a^3/r(x,y)^3*cos(3*θ(x,y)) )
v(x,y) = T*a*(1+ν)/2/E*( -r(x,y)/a*2*ν/(1+ν)*sin(θ(x,y)) - a/r(x,y)*(2*(1-ν)/(1+ν)*sin(θ(x,y))-sin(3*θ(x,y))) - a^3/r(x,y)^3*sin(3*θ(x,y)) )
∂u∂x(x,y) = T/E*(1 + a^2/2/r(x,y)^2*((ν-3)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y)))
∂u∂y(x,y) = T/E*(-a^2/r(x,y)^2*((ν+5)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y)))
∂v∂x(x,y) = T/E*(-a^2/r(x,y)^2*((ν-3)/2*sin(2*θ(x,y))+(1+ν)*sin(4*θ(x,y))) + 3*a^4/2/r(x,y)^4*(1+ν)*sin(4*θ(x,y)))
∂v∂y(x,y) = T/E*(-ν - a^2/2/r(x,y)^2*((1-3*ν)*cos(2*θ(x,y))-2*(1+ν)*cos(4*θ(x,y))) - 3*a^4/2/r(x,y)^4*(1+ν)*cos(4*θ(x,y)))
σ₁₁(x,y) = T - T*a^2/r(x,y)^2*(3/2*cos(2*θ(x,y))+cos(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₂₂(x,y) = -T*a^2/r(x,y)^2*(1/2*cos(2*θ(x,y))-cos(4*θ(x,y))) - T*3*a^4/2/r(x,y)^4*cos(4*θ(x,y))
σ₁₂(x,y) = -T*a^2/r(x,y)^2*(1/2*sin(2*θ(x,y))+sin(4*θ(x,y))) + T*3*a^4/2/r(x,y)^4*sin(4*θ(x,y))

prescribe!(elements["Ω"],:E=>(x,y,z)->E)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Γᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z,n₁,n₂)->(1-abs(n₂))*abs(n₁))
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z,n₁,n₂)->(1-abs(n₁))*abs(n₂))
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

𝑎 = ∫∫εᵢⱼσᵢⱼdxdy=>elements["Ω"]
𝑎ᵅ = ∫σᵢⱼnⱼgᵢds=>elements["Γᵍ"]
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]


k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

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

d₁ = d[1:2:end]
d₂ = d[2:2:end]
push!(nodes,:d₁=>d₁,:d₂=>d₂)


𝐻ₑ, 𝐿₂ = Hₑ_PlaneStress(elements["Ωᵍ"])

println(𝐿₂)
println(𝐻ₑ)
XLSX.openxlsx("./xlsx/platewithhole.xlsx", mode="rw") do xf
index = 4,8,16,32
    Sheet = xf[3]
    ind = findfirst(n->n==ndiv,index)+1
    Sheet["A"*string(ind)] = 3*nₑ
    Sheet["B"*string(ind)] = log10(𝐿₂)
    Sheet["C"*string(ind)] = log10(𝐻ₑ)
end

show(to)