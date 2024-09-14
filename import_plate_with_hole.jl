
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

function import_plate_with_hole_mix(filename1::String,filename2::String,n,c)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    gmsh.initialize()

    gmsh.open(filename1)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    w = 0.0
    for i in 0:n-1
        w += c^i
    end
    ds₂ = 4*2^0.5/w
    ds₁ = ds₂*c^(n-1)
    s = zeros(length(nodes))
    for (i,node) in enumerate(nodes) 
        xᵢ = node.x
        yᵢ = node.y
        r = (xᵢ^2+yᵢ^2)^0.5
        s[i] = ds₁ + (r-1)/4/2^0.5*(ds₂-ds₁)
    end
    s .*= 2.5
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
    elements["Γᵍ"] = getElements(nodes, entities["Γᵍ"],type, integration_Γ, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"],type, integration_Γ, sp, normal = true)

    nₘ = 21
    𝗠 = zeros(nₘ)
    ∂𝗠∂x = zeros(nₘ)
    ∂𝗠∂y = zeros(nₘ)
    push!(elements["Ω"], :𝝭)
    push!(elements["∂Ω"], :𝝭)
    push!(elements["Γᵍ"], :𝝭)
    push!(elements["Γᵗ"], :𝝭)
    push!(elements["Ω"],  :𝗠=>𝗠)
    push!(elements["∂Ω"], :𝗠=>𝗠)
    push!(elements["Γᵍ"], :𝗠=>𝗠)
    push!(elements["Γᵗ"], :𝗠=>𝗠)
    push!(elements["Ωᵍ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["Ωᵍ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)



    type = PiecewisePolynomial{:Linear2D}
    # type = PiecewisePolynomial{:Quadratic2D}
    println(entities)
    elements["Ωˢ"] = getPiecewiseElements(entities["Ω"], type, integration_Ω)
    elements["∂Ωˢ"] = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], type, integration_Γ)
    elements["Γᵍˢ"] = getElements(entities["Γᵍ"],entities["Γ"], elements["∂Ωˢ"])
    elements["Γᵗˢ"] = getElements(entities["Γᵗ"],entities["Γ"], elements["∂Ωˢ"])
    push!(elements["Ωˢ"], :𝝭, :∂𝝭∂x, :∂𝝭∂y)
    push!(elements["∂Ωˢ"], :𝝭)

    # gmsh.finalize()

    return elements, nodes
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    s  = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
        s[i]  = (4/3^0.5*𝐴s[i])^0.5
    end
    return s
end
