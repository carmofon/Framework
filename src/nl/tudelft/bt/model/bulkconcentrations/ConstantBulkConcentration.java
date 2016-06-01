/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;

/**
 * Implements a constant bulk concentration manager class for multidimensional
 * biofilm modelling.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConstantBulkConcentration extends BulkConcentration {

	/**
	 * @param v
	 *            value for constant bulk concentration
	 */
	public ConstantBulkConcentration(float v) {
		super();
		setValue(v);
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		// do nothing, concentration is constant
		setValue(getValue());		
	}
}
