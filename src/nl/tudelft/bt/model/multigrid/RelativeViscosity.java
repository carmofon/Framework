package nl.tudelft.bt.model.multigrid;

import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Implements multigrid variable for the relative viscosity, for computing
 * velocity profiles such as in the rube reactor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class RelativeViscosity extends MultigridVariable {
	/**
	 * @throws MultigridSystemNotSetException
	 */
	public RelativeViscosity() throws MultigridSystemNotSetException {
		super();
		_name = "relativeViscosity";
	}

	/**
	 * Initializes the values of the relative viscosity, which will be 1 for any
	 * location in the liquid and 1.0e20f for any point in biofilm or carrier (very
	 * high, to simulate no velocity flux)
	 * 
	 * @param bac
	 * @param bc
	 */
	public void computeValues(Collection bac, BoundaryConditions bc) {
		//pre-convert collection into array of Particulate species for faster
		// assessment
		ParticulateSpecies[] particulates = new ParticulateSpecies[bac.size()];
		int p = 0;
		for (Iterator iter = bac.iterator(); iter.hasNext();) {
			particulates[p++] = (ParticulateSpecies) iter.next();
		}
		//set the value of the finner grid
		float[][][] d = _mg[_order - 1];
		MultigridUtils.setValues(d, 1.0f);
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d[i].length; j++)
				for (int k = 0; k < d[i][j].length; k++)
					// if inside carrier or biofilm, relative diffusivity should
					// be very high (0 velocity flux)
					if (bc.isCarrier(i, j, k)
							| hasBiomass(particulates, i, j, k)) {
						d[i][j][k] = 1.0e20f;
					}
		//update corser grid copies
		updateMultigridCopies();
	}

	/**
	 * Check if there is any biomass in location i, j, k
	 * 
	 * @param particulates
	 * @param i
	 * @param j
	 * @param k
	 * @return true if there is biomass at this location
	 */
	private boolean hasBiomass(ParticulateSpecies[] particulates, int i, int j,
			int k) {
		for (int index = 0; index < particulates.length; index++) {
			if (particulates[index]._mg[_order - 1][i][j][k] > 0)
				return true;
		}
		return false;
	}
}