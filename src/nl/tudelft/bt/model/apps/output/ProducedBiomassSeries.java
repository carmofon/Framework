/* 
 * Created on 9-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;
import nl.tudelft.bt.model.Model;
/**
 * The total produced biomass
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ProducedBiomassSeries extends VariableSeries {
	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public ProducedBiomassSeries() {
		super("Produced biomass", "Time [T]", "Produced biomass [M]");
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
		//
		if (sizeY < sizeX) {
			getYArray().add(Model.model().biomassContainer.getTotalProducedBiomass());
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}
}