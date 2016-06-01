package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Abstract multigrid variable for the boundary layer implementation
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
abstract public class BoundaryLayer extends MultigridVariable {
	/**
	 * New boundary layer object
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public BoundaryLayer() throws MultigridSystemNotSetException {
		super();
	}

	/**
	 * Compute the value of the boundary layer based on biomass
	 * composition and the boundary conditions
	 * 
	 * @param b
	 * @param bc
	 */
	public abstract void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc);

	public void setThickness(float h) {
	};
}