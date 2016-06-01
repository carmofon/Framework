/* 
 * Created on 18-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;

/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConcentrationTimeSeries extends VariableSeries {
	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public ConcentrationTimeSeries(String name, String ylabel) {
		super(name, "Time [h]", ylabel);
		setX(Model.model().getTimeSeries());
	}	
}
