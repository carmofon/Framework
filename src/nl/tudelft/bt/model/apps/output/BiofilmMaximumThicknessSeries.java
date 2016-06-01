/* 
 * Created on 9-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;

/**
 * The time series
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiofilmMaximumThicknessSeries extends VariableSeries {

	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public BiofilmMaximumThicknessSeries() {
		super("Biofilm thickness", "Time [h]", "Biofilm thickness [um]");
		setX(Model.model().getTimeSeries());
	}

	
	/* (non-Javadoc)
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getLastY()
	 */
	public float getLastY() {
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		if (sizeY < sizeX) {
			for (int i = sizeY; i < sizeX; i++) {
				getYArray().add(Model.model().getCurrentBiofilmHeight());
			}
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

}
