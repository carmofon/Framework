/* 
 * Created on 9-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;
import nl.tudelft.bt.model.Model;
/**
 * The total biomass attached to the carrier
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class FixedTotalBiomassSeries extends VariableSeries {
	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public FixedTotalBiomassSeries() {
		super("Total biomass", "Time [h]",
				"Total biomass in biofilm [M]");
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
			getYArray().add(Model.model().getCurrentTotalBiomass());
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}
}
