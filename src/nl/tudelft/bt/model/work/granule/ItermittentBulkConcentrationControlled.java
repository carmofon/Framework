/*
 * File created originally on Jan 17, 2005
 */
package nl.tudelft.bt.model.work.granule;

import nl.tudelft.bt.model.bulkconcentrations.SbrBulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.SbrOnOffBulkConcentration;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * An ittermitent bulk concentration with the diference that in the feast phase
 * bulk concentration is switched to zero when concetration of 2 other species
 * (speciesToWatch1 and speciesToWatch2) are lower than a threshold value
 * (tresholdValue)
 * 
 * @author jxavier
 */
public class ItermittentBulkConcentrationControlled
		extends
			SbrOnOffBulkConcentration {
	private SoluteSpecies _speciesToWatch1;
	private SoluteSpecies _speciesToWatch2;
	private float _thresholdValue;
	private float _onConcentration;
	private float _offConcentration;

	/**
	 * @param v value of the concentration (when on)
	 * @param timeToTurnOn
	 * @param cycle
	 * @param speciesToWatch1
	 * @param speciesToWatch2
	 * @param tresholdValue
	 */
	public ItermittentBulkConcentrationControlled(float v,
			float timeToTurnOn, SbrBulkConcentration.SbrCycle cycle,
			SoluteSpecies speciesToWatch1, SoluteSpecies speciesToWatch2,
			float tresholdValue) {
		super(v, timeToTurnOn, cycle);
		_speciesToWatch1 = speciesToWatch1;
		_speciesToWatch2 = speciesToWatch2;
		_thresholdValue = tresholdValue;
		_onConcentration = v;
		_offConcentration = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.bulkconcentrations.BulkConcentration#computeBulkConcentration(float)
	 */
	public void computeBulkConcentration(float tstep) {
		// let superclass compute if feast or famine phase
		super.computeBulkConcentration(tstep);
		// get the values of bulk concentrations from species to watch
		float sp1value = _speciesToWatch1.getBulkConcentration();
		float sp2value = _speciesToWatch2.getBulkConcentration();
		// if in famine phase, check if bulk concentration is bellow threshold.
		if ((getValue() == _onConcentration) & (sp1value <= _thresholdValue)
				& (sp2value <= _thresholdValue))
			//if so, turn of concentration to famine value
			setValue(_offConcentration);
	}
}