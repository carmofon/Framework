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
public class StepReturnGoodOnly extends ProcessFactor {
	protected MultigridVariable _publicGood;
	private float _goodThreshold;

	/**
	 * @param c
	 * @param k
	 */
	public StepReturnGoodOnly(MultigridVariable publicGood, float goodThreshold){
		_publicGood = publicGood;
		_goodThreshold = goodThreshold;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concPublicGood = _publicGood.getValue();
		if (concPublicGood > _goodThreshold){
			return 1;
		} else {
			return 0;
		}
	}

	public float getMaximumValue() {
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		return 0;
	}
}
