package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Solute species with static linear gradient in the Y direction
 * 
 * @author wkc
 *
 */
public class StaticGradientYSoluteSpecies extends SoluteSpecies {

	// slope of the static gradient
	private float _slope;
	
	/**
	 * 
	 * @param n	name of species
	 * @param d	diffusivity
	 * @throws MultigridSystemNotSetException
	 */
	public StaticGradientYSoluteSpecies(String n, float d)
			throws MultigridSystemNotSetException {
		super(n, d);
	}
	
	/*
	 *When returning the value calculate based on the current position in the system
	 */
	@Override
	public float getValue() {
		ContinuousCoordinate currentLocation = getCurrentLocation();
		return getBulkConcentration() - (_slope * currentLocation.y);
	}
	
	public void setSlope(float m) {
		_slope = m;
	}
	

}
