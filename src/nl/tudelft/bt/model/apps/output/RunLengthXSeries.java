/* 
 * Created on 9-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassParticleContainer;

/**
 * The time series of runlength of a granule along the X (vertical) direction
 * NOTE: may only be used with granule geometry
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class RunLengthXSeries extends VariableSeries {

	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public RunLengthXSeries() {
		super("RunLengthX", "Time [h]", "Run Length X [L]");
		setX(Model.model().getTimeSeries());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getLastY()
	 */
	public float getLastY() {
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			for (int i = sizeY; i < sizeX; i++) {
				getYArray()
						.add(
								((GranuleBiomassParticleContainer) Model.model().biomassContainer)
										.computeRunLengthX());
			}
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}
}