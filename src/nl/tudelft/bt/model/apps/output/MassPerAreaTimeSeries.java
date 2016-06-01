/* 
 * Created on 18-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;

/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MassPerAreaTimeSeries extends VariableSeries {
	/**
	 * @param name
	 * @param xlabel
	 * @param ylabel
	 */
	public MassPerAreaTimeSeries(String name) {
		super(name + " [M]", "Time [T]", name + " [M]");
		setX(Model.model().getTimeSeries());
	}
}
