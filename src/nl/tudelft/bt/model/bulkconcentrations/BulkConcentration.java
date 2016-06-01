/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;
import java.io.Serializable;

import nl.tudelft.bt.model.apps.output.ConcentrationTimeSeries;
import nl.tudelft.bt.model.apps.output.MassPerAreaTimeSeries;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;
import nl.tudelft.bt.model.util.UnlimitedFloatArray;
/**
 * Abstract class for bulk concentration manager class. Extends VariableSeries
 * to keep track of the time series of bulk concentrations but also provides the
 * global concentration rates through time series
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class BulkConcentration implements Serializable {
	protected SoluteSpecies _species;
	// current value
	private float _value;
	private float _currentGlobalRateFromDiffusionReaction;
	private float _currentGlobalRateFromMassBalance;
	// time series
	private UnlimitedFloatArray _values;
	private UnlimitedFloatArray _globalRatesFromDiffusionReaction;
	private UnlimitedFloatArray _globalRatesFromMassBalances;
	/**
	 * Initializes values vector with 100 slots
	 */
	public BulkConcentration() {
		_value = 0;
		_currentGlobalRateFromDiffusionReaction = 0;
		_currentGlobalRateFromMassBalance = 0;
		_values = new UnlimitedFloatArray();
		_globalRatesFromDiffusionReaction = new UnlimitedFloatArray();
		_globalRatesFromMassBalances = new UnlimitedFloatArray();
	}
	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public abstract void computeBulkConcentration(float tstep);
	/**
	 * @return the value of the bulk concentration
	 */
	public final float getValue() {
		return _value;
	}
	/**
	 * Update the global rate series with a value
	 */
	public void updateGlobalRateFromDiffusionReactionSeries() {
		// compute the global rate and update the time series of rates
		_currentGlobalRateFromDiffusionReaction = _species
				.computeGlobalRateFromDiffusionReaction();
		_globalRatesFromDiffusionReaction
				.add(_currentGlobalRateFromDiffusionReaction);
	}
	/**
	 * Update the global rate series with a value
	 */
	public void updateGlobalRateFromMassBalance() {
		// compute the global rate and update the time series of rates
		_currentGlobalRateFromMassBalance = _species
				.computeGlobalRateFromMassBalance();
		_globalRatesFromMassBalances.add(_currentGlobalRateFromMassBalance);
	}
	/**
	 * @return the maximum time step in order not to have negative cocentrations
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		updateGlobalRateFromDiffusionReactionSeries();
		// default is positive infinity, which means no constraint is set by
		// default
		return new TimeStepConstraint();
	}
	/**
	 * Set the handle to the solute species that will use this BulkConcentration
	 * instance
	 * 
	 * @param s
	 *            The _species to set.
	 */
	public void setSpecies(SoluteSpecies s) {
		_species = s;
	}
	/**
	 * Updates the value and adds new value to the _values array for time series
	 * update.
	 * 
	 * @param v
	 *            The value to set.
	 */
	public void setValue(float v) {
		_value = v;
		_values.add(v);
	}
	/**
	 * Get the global convertion rates [M.L-2.T-1] time series
	 * 
	 * @return global rate time series
	 */
	public VariableSeries getRateTimeSeries() {
		MassPerAreaTimeSeries globalRatesSeries = new MassPerAreaTimeSeries(
				_species.getName());
		globalRatesSeries.setY(_globalRatesFromMassBalances);
		return globalRatesSeries;
	}
	/**
	 * Get the global convertion rates [M.L-2.T-1] time series
	 * 
	 * @return global rate time series
	 */
	public VariableSeries getBulkConcentrationTimeSeries() {
		String label = _species.getName() + " bulk concentration [M/L^3]";
		ConcentrationTimeSeries series = new ConcentrationTimeSeries(label,
				label);
		series.setY(_values);
		return series;
	}
	/**
	 * Get the the global rate for the whole computational volume without
	 * recomputing
	 * 
	 * @return the conversion rates for the entire computational volume
	 */
	public float getCurrentGlobalRateFromDiffusionReaction() {
		return _currentGlobalRateFromDiffusionReaction;
	}
	/**
	 * Get the the global rate for the whole computational volume without
	 * recomputing
	 * 
	 * @return the conversion rates for the entire computational volume
	 */
	public float getCurrentGlobalRateFromMassBalance() {
		return _currentGlobalRateFromMassBalance;
	}
	/**
	 * Gets the number of iterations already computed for this instance of
	 * BulkConcentration.
	 * 
	 * @return the number of iterations already passed
	 */
	public int getNumberOfIterationsComputed() {
		return _values.getSize();
	}
}