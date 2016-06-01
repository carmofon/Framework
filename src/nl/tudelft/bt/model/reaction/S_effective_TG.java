/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class S_effective_TG extends ProcessFactor {
	protected SoluteSpecies _solute1;
	protected SoluteSpecies _gf1;
	protected SoluteSpecies _gf2;
	private float _ks1;
	private float _kg;
	private float _fraction;

	/**
	 * @param c
	 * @param k
	 */
	public S_effective_TG(
		SoluteSpecies c1,
		SoluteSpecies g1,
		SoluteSpecies g2,
		float ks1,
		float kg,
		float fraction) {
		_solute1 = c1;
		_gf1=g1;
		_gf2=g2;
		_ks1 = ks1;
		_kg=kg;
		_fraction = fraction;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc1 = _solute1.getValue();
		float concg1 = _gf1.getValue();
		float concg2 = _gf2.getValue();
		return (_fraction*conc1 + (1-_fraction)*conc1*(concg1*concg2)/(concg1*concg2+_kg)) / (_ks1+_fraction*conc1 + (1-_fraction)*conc1*(concg1*concg2)/(concg1*concg2+_kg));
	}

	public float getMaximumValue() {
		float conc1 = _solute1.getMaximumValue();
		float concg1 = _gf1.getMaximumValue();
		float concg2 = _gf2.getMaximumValue();
		return (_fraction*conc1 + (1-_fraction)*conc1*(concg1*concg2)/(concg1*concg2+_kg)) / (_ks1+_fraction*conc1 + (1-_fraction)*conc1*(concg1*concg2)/(concg1*concg2+_kg));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _solute1) {
			float conc1 = _solute1.getValue();
			float concg1 = _gf1.getValue();
			float concg2 = _gf2.getValue();
			return (_ks1*(concg1*concg2 +_fraction*_kg)*(_kg + concg1*concg2))/((_kg*_ks1 + conc1*concg1*concg2 + concg1*concg2*_ks1 + conc1*_fraction*_kg)*(_kg*_ks1 + conc1*concg1*concg2 + concg1*concg2*_ks1 + conc1*_fraction*_kg));						
		} else {
			return 0f;
		}
	}
}
