
using Gmsh, Statistics

const lobatto2 = ([-1.0,0.0,0.0,
                    1.0,0.0,0.0],[1.0,1.0])

const lobatto3 = ([-1.0,0.0,0.0,
                    0.0,0.0,0.0,
                    1.0,0.0,0.0],[1/3,4/3,1/3])

const trilobatto3 = ([0.0000000000000000,0.5000000000000000,0.0,
                      0.5000000000000000,0.0000000000000000,0.0,
                      0.5000000000000000,0.5000000000000000,0.0],
                   0.5*[1/3,1/3,1/3])
function import_patchtest_mix(filename1::String,filename2::String)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    Ω = getElements(nodes, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes))
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

    integration_Ω = 2
    integrationOrder_Ωᵍ = 8
    integration_Γ = 2

    gmsh.open(filename2)
    entities = getPhysicalGroups()

    # type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integration_Ω, sp)
    elements["∂Ω"] = getElements(nodes, entities["Γ"], type, integration_Γ, sp, normal = true)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], type, integrationOrder_Ωᵍ, sp)
    elements["Γ¹"] = getElements(nodes, entities["Γ¹"],type, integration_Γ, sp, normal = true)
    elements["Γ²"] = getElements(nodes, entities["Γ²"],type, integration_Γ, sp, normal = true)
    elements["Γ³"] = getElements(nodes, entities["Γ³"],type, integration_Γ, sp, normal = true)
    elements["Γ⁴"] = getElements(nodes, entities["Γ⁴"], type, integration_Γ, sp, normal = true)
    elements["Γ"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ω"], :𝝭)
    push!(elements["∂Ω"], :𝝭)
    push!(elements["Γ¹"], :𝝭)
    push!(elements["Γ²"], :𝝭)
    push!(elements["Γ³"], :𝝭)
    push!(elements["Γ⁴"], :𝝭)
    push!(elements["Ω"],  :𝗠=>𝗠)
    push!(elements["∂Ω"], :𝗠=>𝗠)
    push!(elements["Γ¹"], :𝗠=>𝗠)
    push!(elements["Γ²"], :𝗠=>𝗠)
    push!(elements["Γ³"], :𝗠=>𝗠)
    push!(elements["Γ⁴"], :𝗠=>𝗠)
    push!(elements["Ωᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    set𝝭!(elements["Ω"])
    set𝝭!(elements["∂Ω"])
    set∇𝝭!(elements["Ωᵍ"])
    set𝝭!(elements["Γ"])

    type = PiecewisePolynomial{:Linear2D}
    # type = PiecewisePolynomial{:Quadratic2D}
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integration_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integration_Γ)
    elements["Γ¹ˢ"] = getElements(entities["Γ¹"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ²ˢ"] = getElements(entities["Γ²"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ³ˢ"] = getElements(entities["Γ³"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γ⁴ˢ"] = getElements(entities["Γ⁴"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γˢ"] = elements["Γ¹ˢ"]∪elements["Γ²ˢ"]∪elements["Γ³ˢ"]∪elements["Γ⁴ˢ"]
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    set∇𝝭!(elements["Ωˢ"])
    set𝝭!(elements["∂Ωˢ"])

    # gmsh.finalize()

    return elements, nodes
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4/3^0.5*avg𝐴)^0.5
    return s, var𝐴
end
