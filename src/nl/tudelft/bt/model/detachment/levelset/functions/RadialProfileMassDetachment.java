/* 
 * Created on 10-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset.functions;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;
import nl.tudelft.bt.model.profiles1d.Radial2DMultigridProfile;

/**
 * Implements mass detachment that is proportional to the square of the inverse
 * of the distance to the center of the tube
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class RadialProfileMassDetachment extends DetachmentSpeedFunction {
	float _detachmentRateConstant;
	Radial2DMultigridProfile _profile;

	/**
	 * @param rate
	 */
	public RadialProfileMassDetachment(float rate, Radial2DMultigridProfile p) {
		_detachmentRateConstant = rate;
		_profile = p;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate,
	 *      float)
	 */
	public float getValue(ContinuousCoordinate c) {
		float dist = ((TubeBiomassParticleContainer) Model.model().biomassContainer)
				.rConvertToPolar(c.x, c.y);
		float density = Model.model().biomassContainer.getElementDensity(c);
		return _detachmentRateConstant / density * _profile.getValueAt(dist);
	}

	/**
	 * Set the detachment rate constant
	 * 
	 * @param r
	 *            The _detachmentRate to set.
	 */
	public void setDetachmentRateConstant(float r) {
		_detachmentRateConstant = r;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.detachment.DetachmentSpeedFunction#detachmentIsOff()
	 */
	public boolean detachmentIsOff() {
		return _detachmentRateConstant == 0;
	}
}