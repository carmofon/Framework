/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.vanni;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements a zeroth order dependency of a species S. The difference between
 * zeroth order dependency and no dependency at all is that if concentration of
 * the species S is zero, the rate value is also 0.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Monod_plus_Good extends ProcessFactor {
	protected MultigridVariable _substrateSpecies;
	protected MultigridVariable _goodSpecies;
	protected float _halfsaturation;
	protected float _inversGpar;

	/**
	 * @param c
	 *            species that defines the factor
	 */
	public Monod_plus_Good(MultigridVariable S, MultigridVariable G, float KN, float B) {
		_substrateSpecies = S;
		_goodSpecies = G;
		_halfsaturation = KN;
		_inversGpar=B;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concS = _substrateSpecies.getValue();
		float concG = _goodSpecies.getValue();
		concS = (concS < 0 ? 0 : concS);
		concG = (concG < 0 ? 0 : concG);
		//return (concS*concG == 0 ? 0 : concS*concG);
		//return concS/(concS+_halfsaturation)+1*_inversGpar*concG;
		return concS/(concS+_halfsaturation)+((_inversGpar*concG*concS)/(concS+_halfsaturation));
	}

	public float getMaximumValue() {
		float concS = _substrateSpecies.getMaximumValue();
		float concG = _goodSpecies.getMaximumValue();
		//return (concS*concG == 0 ? 0 : concS*concG);
		concS = (concS < 0 ? 0 : concS);
		concG = (concG < 0 ? 0 : concG);
		//return concS/(concS+_halfsaturation)+1*_inversGpar*concG;
		return concS/(concS+_halfsaturation)+((_inversGpar*concG*concS)/(concS+_halfsaturation));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _substrateSpecies){
			float concS = _substrateSpecies.getValue();
			float concG = _goodSpecies.getValue();
			//return concG;
			concS = (concS < 0 ? 0 : concS);
			concG = (concG < 0 ? 0 : concG);
			//return _halfsaturation / ((_halfsaturation + concS) * (_halfsaturation + concS));
			return _halfsaturation / ((_halfsaturation + concS) * (_halfsaturation + concS))+(1*_inversGpar*concG* _halfsaturation)/ ((_halfsaturation + concS) * (_halfsaturation + concS));
		} else if (c == _goodSpecies){
			float concS = _substrateSpecies.getValue();
			float concG = _goodSpecies.getValue();
			//return concS;
			//return _inversGpar;
			return _inversGpar*concS/(concS+_halfsaturation);
		}
		return 0f;
	}
}