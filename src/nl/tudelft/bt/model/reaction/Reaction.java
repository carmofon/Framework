/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements a reaction rate
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Reaction implements Serializable {
	private String _name;

	private Species _catalyst;

	private float _constant;

	private ProcessFactor[] _factors;

	private int _addCounter;

	private float _globalReactionRate;

	private float _presentReactionRate;

	/**
	 * creates a new reaction with one or more factors
	 * 
	 * @param n
	 *            name of reaction
	 * @param bs
	 *            biomass species that catalyzes this reaction
	 * @param k
	 *            rate constant (for example umax)
	 * @param nfactors
	 *            number of factors in this reaction
	 */
	public Reaction(String n, Species bs, float k, int nfactors) {
		_name = n;
		_catalyst = bs;
		_constant = k;
		_factors = new ProcessFactor[nfactors];
		_addCounter = 0;
		_globalReactionRate = 0;
		// add the reaction thus created to the Model._reactions
		Model.model().addReaction(this);
	}

	/**
	 * Add a reaction factor to the array of factors in this recation
	 * 
	 * @param _f
	 */
	public void addFactor(ProcessFactor _f) {
		try {
			_factors[_addCounter++] = _f;
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ModelRuntimeException("Trying to add a " + _addCounter
					+ " factor to " + _name);
		}
	}

	/**
	 * Returns the rate at the present location
	 * 
	 * @return rate of reaction [g/um^3/h]
	 */
	public float getRate() {
		return getSpecificRateFactor() * _catalyst.getValue();
	}

	/**
	 * Pre-compute the mass-based growth rate and store it in atribute
	 * _presentReactionRate. Also, update the value of _globalReactionRate
	 * 
	 * @param mCatalyst
	 *            mass of catalyst species
	 */
	public void computeMassGrowthRateAndAddToGlobal(BiomassSpecies.Composition c) {
		float mCatalyst = c.getSpeciesMass((ParticulateSpecies) _catalyst);
		float rateFactor = getSpecificRateFactor(c);
		// Check if the values are valid
		if ((mCatalyst != mCatalyst) | (rateFactor != rateFactor)) {
			String st = "reaction " + _name + "produced NaN (mCatalyst = "
					+ mCatalyst + " of " + _catalyst.getName()
					+ ", rateFactor = " + rateFactor + ")";
			throw new ModelRuntimeException(st);
		}
		_presentReactionRate = rateFactor * mCatalyst;
		// to update the value of the global reaction rate
		// it must be taken into account the numeric integration
		// used for Composition.growWithIntegratedPrecision (NINT)
		_globalReactionRate += _presentReactionRate / BiomassSpecies.NINT;
	}

	/**
	 * Get the value of _presentReactionRate.
	 * 
	 * @return the mass growth rate [g/h]
	 */
	public float getPreComputedMassGrowthRate() {
		return _presentReactionRate;
	}

	/**
	 * Returns the specific rate at the present location
	 * 
	 * @return specific rate of reaction [h^-1]
	 */
	public float getSpecificRateFactor() {
		float r = _constant;
		for (int i = 0; i < _addCounter; i++) {
			r *= _factors[i].getValue();
		}
		return r;
	}

	/**
	 * Returns the specific rate at the present location
	 * 
	 * @param c
	 *            composition of biomass particle to compute specific rate
	 * @return specific rate of reaction [h^-1]
	 */
	public float getSpecificRateFactor(BiomassSpecies.Composition c) {
		float r = _constant;
		for (int i = 0; i < _addCounter; i++) {
			r *= _factors[i].getValue(c);
		}
		return r;
	}

	/**
	 * Returns the maximum specific rate in the system
	 * 
	 * @return maximum specific rate in the system [h^-1]
	 */
	public float getSpecificRateMaximum() {
		float r = _constant;
		for (int i = 0; i < _addCounter; i++) {
			r *= _factors[i].getMaximumValue();
		}
		return r;
	}

	/**
	 * Returns the derivative of the rate of a given reaction at the present
	 * location
	 * 
	 * @param c
	 *            chemical species to derivate rate to
	 * @return derivative of rate of reaction [g/um^3/h]
	 */
	public float getRateDerivative(SoluteSpecies c) {
		float dr = 0;
		// implements derivative of multiplying functions
		for (int i = 0; i < _addCounter; i++) {
			float v = 1;
			v *= _factors[i].getDerivative(c);
			for (int j = 0; j < _addCounter; j++)
				if (j != i)
					v *= _factors[j].getValue();
			dr += v;
		}
		if (c == _catalyst) {
			// in case the solute species is the catalyst of the reaciton
			float r = _constant;
			for (int i = 0; i < _addCounter; i++)
				r *= _factors[i].getValue();
			return dr * _constant * _catalyst.getValue() + r;
		}
		// if the species is not the catalyst
		return dr * _constant * _catalyst.getValue();
	}

	/**
	 * Update the array with values of the rate and rate derivative
	 * 
	 * @param c
	 *            the solute species for which this is being computed
	 * @param rDr
	 *            [rate, rateDerivative]
	 */
	public void updateValuesForRateAndRateDerivative(SoluteSpecies c,
			float[] rDr) {
		rDr[0] = _constant;
		rDr[1] = 0;
		float catalystValue = _catalyst.getValue();
		//compute the rates and store the values of the rates for each
		//factor for later use
		for (int i = 0; i < _addCounter; i++) {
			_factors[i].preComputedValue = _factors[i].getValue();
			rDr[0] *= _factors[i].preComputedValue;
		}
		rDr[0] *= catalystValue; //rate is finished computing
		// Compute the derivative implemented by multiplying functions
		for (int i = 0; i < _addCounter; i++) {
			float v = 1;
			v *= _factors[i].getDerivative(c);
			for (int j = 0; j < _addCounter; j++)
				if (j != i)
					v *= _factors[j].preComputedValue;
			rDr[1] += v;
		}
		rDr[1] *= _constant * catalystValue;
		// if the species is the catalyst		
		if (c == _catalyst) {
			// in case the solute species is the catalyst of the reaciton
			float r = _constant;
			for (int i = 0; i < _addCounter; i++)
				r *= _factors[i].preComputedValue;
			rDr[1] += r;
			return;
		}
	}

	/**
	 * @return Returns the _catalyst.
	 */
	public Species getCatalystSpecies() {
		return _catalyst;
	}

	/**
	 * Set the reaction constant rate
	 * 
	 * @param _constant
	 *            The _constant to set.
	 */
	public void setConstant(float constant) {
		this._constant = constant;
	}

	/**
	 * Set the value of the global reaction rate to 0. This method is called at
	 * the begining of the Model.performGrowthAndDivision()
	 */
	public void resetGlobalReactionRate() {
		_globalReactionRate = 0;
	}

	/**
	 * Get the the current value of the global reaction rate
	 * 
	 * @return Returns the _globalReactionRate.
	 */
	public float getGlobalReactionRate() {
		return _globalReactionRate;
	}
}