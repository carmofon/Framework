package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Calculates the process factor by Hill (non-Michaelis-Menten) kinetics
 * 
 * P = (species ^ n) / (k ^ n + species ^ n)
 * 
 * @author wkc
 *
 */
public class Hill extends ProcessFactor {
	
	protected MultigridVariable _species;
	
	private float _n;
	
	private float _k;
	
	/**
	 * 
	 * @param species	the reactant species
	 * @param constant	the reaction constant
	 * @param coefficient	the Hill coefficient
	 */
	public Hill(MultigridVariable species, float constant, float coefficient) {
		_species = species;
		_k = constant;
		_n = coefficient;
	}

	@Override
	public float getValue() {
		float conc = _species.getValue();
		conc = (conc < 0 ? 0 : conc);
		double p = Math.pow(conc, _n) / (Math.pow(_k, _n) + Math.pow(conc, _n));
		return (float) p;
	}

	@Override
	public float getDerivative(SoluteSpecies c) {
		if (c == _species) {
			float conc = _species.getValue();
			conc = (conc < 0 ? 0 : conc);
			double p = (_n * Math.pow(_k, _n) * Math.pow(conc, _n - 1)) 
				/ (Math.pow(Math.pow(_k, _n) + Math.pow(conc, _n), 2));
			return (float) p;
		}
		else {
			return 0f;
		}
	}

	@Override
	public float getMaximumValue() {
		float conc = _species.getMaximumValue();
		conc = (conc < 0 ? 0 : conc);
		float p = (float) Math.pow(conc, _n) / (float) (Math.pow(_k, _n) + Math.pow(conc, _n));
		return p;
	}

}
