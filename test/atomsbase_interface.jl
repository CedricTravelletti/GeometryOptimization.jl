#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "AtomsBase interface: check position updating" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using LinearAlgebra
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = collect.(orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions])
    
    strain =  zeros(6) 
    # Compute deformed cell manually for checks.
    deformation_tensor = I + GeometryOptimization.voigt_to_full(strain)
    new_lattice = eachcol(deformation_tensor
                          * GeometryOptimization.bbox_to_matrix(bounding_box(silicon_supercell)))
    
    new_general_pos = ComponentVector(atoms = new_positions, strain = strain)
    new_system = update_positions(silicon_supercell, new_general_pos)
    
    @test position(new_system) == new_positions
    @test bounding_box(new_system) == bounding_box(silicon_supercell)
end

@testitem "AtomsBase interface: check positions+cell updating" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using LinearAlgebra
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = collect.(orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions])
    
    strain = 0.1 * rand(Float64, 6) # 10% random strain.
    # Compute deformed cell manually for checks.
    deformation_tensor = I + GeometryOptimization.voigt_to_full(strain)
    new_lattice = eachcol(deformation_tensor
                          * GeometryOptimization.bbox_to_matrix(bounding_box(silicon_supercell)))
    deformed_new_positions = [deformation_tensor * position for position in new_positions]
    
    new_general_pos = ComponentVector(atoms = new_positions, strain = strain)
    new_system = update_positions(silicon_supercell, new_general_pos)
    
    @test position(new_system) == deformed_new_positions
    @test bounding_box(new_system) == new_lattice
end

@testitem "AtomsBase interface: check positions+cell updating (with mask)" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using LinearAlgebra
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    orig_lattice = bounding_box(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = collect.(orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions])
    new_not_clamped_positions = reduce(vcat, new_positions)

    strain = 0.1 * rand(Float64, 6) # 10% random strain.
    # Compute deformed cell manually for checks.
    deformation_tensor = I + GeometryOptimization.voigt_to_full(strain)
    new_lattice = eachcol(deformation_tensor
                          * GeometryOptimization.bbox_to_matrix(bounding_box(silicon_supercell)))
    deformed_new_positions = [deformation_tensor * position for position in new_positions]
    
    new_general_pos = ComponentVector(atoms = new_not_clamped_positions,
                                      strain = strain)
    new_system = update_not_clamped_positions(silicon_supercell, new_general_pos)
    
    @test position(new_system) == deformed_new_positions
    @test bounding_box(new_system) == new_lattice
end
