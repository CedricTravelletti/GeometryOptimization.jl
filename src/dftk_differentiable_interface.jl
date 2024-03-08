# 
# Interface for updating positions of lattice in a DFTK system. 
# The goal of this interface is to be differentiable, 
# so that derivatives of DFTK quantities with respect to positions or strain 
# can be taken.

using DFTK
using LinearAlgebra


function update_positions(model::DFTK.Model, positions, strain)
    """ Update positions and strain of DFTK model in a differentiable way.

    """
    deformation_tensor = I + voigt_to_full(strain)
    positions = [deformation_tensor * x for x in positions]
    lattice = deformation_tensor * model.lattice
    Model(model; lattice, positions)
end
