#=
# Note that by default all particles in the system are assumed optimizable.
# IMPORTANT: Note that we always work in cartesian coordinates.
=#
using DFTK
export minimize_energy!


function update_state(original_system, new_system, state)
    if bounding_box(original_system) != bounding_box(new_system)
        return DFTK.DFTKState()
    else
        return state
    end
end
    
""" 
By default we work in cartesian coordinates.
Note that internally, when optimizing the cartesian positions, atomic units 
are used.
"""
function Optimization.OptimizationFunction(system, calculator; kwargs...)
    mask = not_clamped_mask(system)  # mask is assumed not to change during optim.

    # TODO: Note that this function will dispatch appropriately when called with 
    # a Component vector.
    f = function(x, p)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        state = update_state(system, new_system, calculator.state)
        energy = AtomsCalculators.potential_energy(new_system, calculator; state, kwargs...)
        austrip(energy)
    end

    function g!(G, x, p)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        # TODO: Determine if here we need a call to update_state.
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)

        forces = AtomsCalculators.forces(new_system, calculator; kwargs...)
        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = collect(Iterators.flatten(forces[mask]))

        # NOTE: minus sign since forces are opposite to gradient.
        G .= - austrip.(forces_concat)
    end
    function g!(G::ComponentVector, x::ComponentVector, p)
        deformation_tensor = I + voigt_to_full(austrip.(x.strain))
        new_system = update_not_clamped_positions(system, x * u"bohr")
        
        state = update_state(system, new_system, calculator.state)
        forces = AtomsCalculators.forces(new_system, calculator; state, kwargs...)
        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = collect(Iterators.flatten([deformation_tensor * f for f in forces[mask]]))

        # NOTE: minus sign since forces are opposite to gradient.
        G.atoms .= - austrip.(forces_concat)
        virial = AtomsCalculators.virial(new_system, calculator)
        G.strain .= - full_to_voigt(virial / deformation_tensor)
    end
    OptimizationFunction(f; grad=g!)
end

@doc raw"""
    minimize_energy!(system, calculator; [relaxation, optimizer])

Minimize the system's energy with respect to it's internal degrees of freedom 
(geometry optimization). Degrees of freedom can be either the atomic coordinates 
only, or the atomic coordinates together with the unit cell vectors.

Keyword arguments:
- `relaxation`: Defines if only atoms or atoms and unit cell vectors are optimized.
                Must be one of `["atoms", "unit_cell"]`.
- `optimizer`: Optimizer to use. Use one from Optimization.jl. Defaults to LBFGS.

Returns a named tuple (; optimized_system, optimization_data).

"""
function minimize_energy!(system, calculator; relaxation="atoms", optimizer=Optim.LBFGS(), kwargs...)
    # Use current system parameters as starting positions.
    if relaxation == "atoms"
        x0 = austrip.(not_clamped_positions(system))
    elseif relaxation == "unit_cell"
        x0 = ComponentVector(atoms = austrip.(reduce(vcat, position(system))), strain = zeros(6))
    else
        print("Error: unknown relaxation type. Please use one of [`atoms`, `unit_cell`].")
    end
    f_opt = OptimizationFunction(system, calculator)
    problem = OptimizationProblem(f_opt, x0, nothing)  # Last argument needed in Optimization.jl.
    optimization_data = solve(problem, optimizer; kwargs...)

    # Return final structure.
    optimized_system = update_not_clamped_positions(system, optimization_data.u * u"bohr")
    (; optimized_system, optimization_data)
end
