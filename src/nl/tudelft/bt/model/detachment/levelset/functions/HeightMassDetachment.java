/* 
 * Created on 10-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset.functions;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;

/**
 * Implements mass detachment that is proportional to height
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class HeightMassDetachment
	extends DetachmentSpeedFunction {
	float _detachmentRate;
	
	/**
	 * @param rate
	 */
	public HeightMassDetachment(float rate) {
		_detachmentRate = rate;
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate, float)
	 */
	public float getValue(ContinuousCoordinate c) {
		float density = Model.model().biomassContainer.getElementDensity(c);
		return _detachmentRate*c.x/density;
	}
	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.detachment.DetachmentSpeedFunction#detachmentIsOff()
	 */
	public boolean detachmentIsOff() {
		return _detachmentRate == 0;
	}
}
