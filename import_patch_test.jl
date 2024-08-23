
using Tensors, BenchmarkExample
import Gmsh: gmsh

function import_patchtest_mix(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    integrationOrder_Ω = 2
    integrationOrder_Γ = 2
    integrationOrder_Ωᵍ = 10

    gmsh.initialize()

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_p = get𝑿ᵢ()
    xᵖ = nodes_p.x
    yᵖ = nodes_p.y
    zᵖ = nodes_p.z
    Ω = getElements(nodes_p, entities["Ω"])
    # s, var𝐴 = cal_area_support(Ω)
    # s = 1.5*s*ones(length(nodes_p))
    s = 2.5*ones(length(nodes_p))
    push!(nodes_p,:s₁=>s,:s₂=>s,:s₃=>s)

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder_Ω)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"], integrationOrder_Γ, normal = true)
    elements["Γ²"] = getElements(nodes, entities["Γ²"], integrationOrder_Γ, normal = true)
    elements["Γ³"] = getElements(nodes, entities["Γ³"], integrationOrder_Γ, normal = true)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], integrationOrder_Γ, normal = true)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ¹"], :𝝭=>:𝑠)
    push!(elements["Γ²"], :𝝭=>:𝑠)
    push!(elements["Γ³"], :𝝭=>:𝑠)
    push!(elements["Γ⁴"], :𝝭=>:𝑠)

    # type = PiecewisePolynomial{:Constant2D}
    type = PiecewisePolynomial{:Linear2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integrationOrder_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integrationOrder_Γ)
    elements["Γ¹ˢ"] = getElements(entities["Γ¹"], entities["Γ"], elements["∂Ωˢ"])
    elements["Γ²ˢ"] = getElements(entities["Γ²"], entities["Γ"], elements["∂Ωˢ"])
    elements["Γ³ˢ"] = getElements(entities["Γ³"], entities["Γ"], elements["∂Ωˢ"])
    elements["Γ⁴ˢ"] = getElements(entities["Γ⁴"], entities["Γ"], elements["∂Ωˢ"])
    elements["Γˢ"] = elements["Γ¹ˢ"]∪elements["Γ²ˢ"]∪elements["Γ³ˢ"]∪elements["Γ⁴ˢ"]
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["∂Ωˢ"], :𝝭=>:𝑠)

    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    sp = RegularGrid(xᵖ,yᵖ,zᵖ,n = 3,γ = 5)
    elements["Ωᵖ"] = getElements(nodes_p, entities["Ω"], type, integrationOrder_Ω, sp)
    elements["∂Ωᵖ"] = getElements(nodes_p, entities["Γ"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ¹ᵖ"] = getElements(nodes_p, entities["Γ¹"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ²ᵖ"] = getElements(nodes_p, entities["Γ²"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ³ᵖ"] = getElements(nodes_p, entities["Γ³"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γ⁴ᵖ"] = getElements(nodes_p, entities["Γ⁴"], type, integrationOrder_Γ, sp, normal = true)
    elements["Γᵖ"] = elements["Γ¹ᵖ"]∪elements["Γ²ᵖ"]∪elements["Γ³ᵖ"]∪elements["Γ⁴ᵖ"]

    nₘ = 6
    # nₘ = 21
    𝗠 = (0,zeros(nₘ))
    push!(elements["Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["Ωᵖ"], :𝗠=>𝗠)
    push!(elements["∂Ωᵖ"], :𝝭=>:𝑠)
    push!(elements["∂Ωᵖ"], :𝗠=>𝗠)
    push!(elements["Γ¹ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ¹ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ²ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ²ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ³ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ³ᵖ"], :𝗠=>𝗠)
    push!(elements["Γ⁴ᵖ"], :𝝭=>:𝑠)
    push!(elements["Γ⁴ᵖ"], :𝗠=>𝗠)


    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    return elements, nodes, nodes_p
end

prescribeForFem = quote
    prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
    prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
    prescribe!(elements["Γ¹"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γ²"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γ³"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γ⁴"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γ¹"],:t₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ²"],:t₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ³"],:t₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ⁴"],:t₂=>(x,y,z)->0.0)
    prescribe!(elements["Γ¹ˢ"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ²ˢ"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ³ˢ"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ⁴ˢ"],:g₁=>(x,y,z)->u(x,y))
    prescribe!(elements["Γ¹ˢ"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ²ˢ"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ³ˢ"],:g₂=>(x,y,z)->v(x,y))
    prescribe!(elements["Γ⁴ˢ"],:g₂=>(x,y,z)->v(x,y))
end
