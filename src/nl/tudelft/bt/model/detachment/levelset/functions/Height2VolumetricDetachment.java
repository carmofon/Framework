/* 
 * Created on 10-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset.functions;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements volumetric detachment that is proportional to height2
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Height2VolumetricDetachment
	extends DetachmentSpeedFunction {
	float _detachmentRate;
	
	/**
	 * @param rate
	 */
	public Height2VolumetricDetachment(float rate) {
		_detachmentRate = rate;
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate, float)
	 */
	public float getValue(ContinuousCoordinate c) {
		return _detachmentRate*ExtraMath.sq(c.x);
	}
	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.detachment.DetachmentSpeedFunction#detachmentIsOff()
	 */
	public boolean detachmentIsOff() {
		return _detachmentRate == 0;
	}
}
