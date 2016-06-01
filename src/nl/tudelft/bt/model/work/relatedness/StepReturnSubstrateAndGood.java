/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.relatedness;

import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class StepReturnSubstrateAndGood extends ProcessFactor {
	protected MultigridVariable _substrate;

	protected MultigridVariable _publicGood;

	private float _Ks;

	private float _goodThreshold;

	/**
	 * @param c
	 * @param k
	 */
	public StepReturnSubstrateAndGood(MultigridVariable substrate,
			MultigridVariable publicGood, float Ks, float goodThreshold) {
		_substrate = substrate;
		_publicGood = publicGood;
		_Ks = Ks;
		_goodThreshold = goodThreshold;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concSubstrate = _substrate.getValue();
		float concPublicGood = _publicGood.getValue();
		if (concPublicGood < _goodThreshold) {
			return 0;
		} else {
			return concSubstrate / (_Ks + concSubstrate);
		}
	}

	public float getMaximumValue() {
		float concSubstrate = _substrate.getMaximumValue();
		float concPublicGood = _publicGood.getMaximumValue();
		concPublicGood = (concPublicGood < 0 ? 0 : concPublicGood);
		return concSubstrate / (_Ks + concSubstrate);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		float concPublicGood = _publicGood.getValue();
		if (concPublicGood < _goodThreshold)
			return 0;
		if (c == _substrate) {
			float conc1 = _substrate.getValue();
			return _Ks / (_Ks + conc1) / (_Ks + conc1);
		} else if (c == _publicGood) {
			return 0;
		}
		return 0f; // derivative for other chemicals
	}
}
