#= The differentiable DFTK interface.
=#
@testitem "DFTK differentiable interface: check position updating" setup=[TestCases] begin
    using DFTK
    using GeometryOptimization

    silicon = TestCases.silicon
    model_kwargs = (; temperature=1e-6, functionals=[:lda_x, :lda_c_pw])
    model = model_DFT(silicon; symmetries=false, model_kwargs...)
    
    strain =  [0.1, 0, 0, 0, 0, 0]
    new_model = update_positions(model, model.positions, strain)
    
    @test isapprox(new_model.lattice[1, :], 1.1 * model.lattice[1, :]; atol=1e-5)
    @test isapprox(new_model.positions[1][1], 1.1 * model.positions[1][1]; atol=1e-5)
end
