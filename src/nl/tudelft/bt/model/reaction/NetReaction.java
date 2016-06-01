/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Implements a container for processes and their sthoichiometric coeficients
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class NetReaction implements Serializable {
	private Reaction[] _reactions;

	private float[] _coeficients;

	private int _addCounter;
	
	//auxiliary variable
	static final float[] _reactionRDr = new float[2];

	/**
	 * creates a new instance o ReactionStoichiometry with space for n reactions
	 * 
	 * @param n
	 *            number of reactions
	 */
	public NetReaction(int n) {
		_reactions = new Reaction[n];
		_coeficients = new float[n];
		_addCounter = 0;
	}

	/**
	 * add a new reaction to the stoichiometry container
	 * 
	 * @param r
	 *            reaction to add
	 * @param coef
	 *            stoichiometry coefficient
	 */
	public void addReaction(Reaction r, float coef) {
		try {
			_reactions[_addCounter] = r;
			_coeficients[_addCounter++] = coef;
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ModelRuntimeException("Trying to add a "
					+ (_addCounter + 1) + " reaction to stoichiometry");
		}
	}

	/**
	 * Computes the net sum of processe factors and their stoichiometric
	 * coeficients
	 * 
	 * @return net sum of rate factors [1/h]
	 */
	public float getSpecificRate() {
		float u = 0;
		for (int i = 0; i < _addCounter; i++) {
			u += _coeficients[i] * _reactions[i].getSpecificRateFactor();
		}
		return u;
	}

	/**
	 * Computes the net sum of processe factors and their stoichiometric
	 * coeficients (for bulk conentration of chemicals)
	 * 
	 * @return the maximum of the factor net sum
	 */
	public float getSpecificRateMaximum() {
		float u = 0;
		for (int i = 0; i < _addCounter; i++) {
			// NOTE: if rate is inhibitting (coeficient
			// less than 0) factor is considered 0 for purposes of colloring
			// particles with colors proportional to growth rate
			if (_coeficients[i] > 0)
				u += _coeficients[i] * _reactions[i].getSpecificRateMaximum();
		}
		return u;
	}

	/**
	 * Get the mass growth rate for a particle from stoichiopmetry calculations
	 * of pre-computed reaction rates
	 * 
	 * @return the mass rate for this stoichiometry [g/h]
	 */
	public float computeMassRateFromPreComputedReactionRates() {
		float r = 0;
		for (int i = 0; i < _addCounter; i++) {
			r += _coeficients[i] * _reactions[i].getPreComputedMassGrowthRate();
		}
		return r;
	}
	/**
	 * Computes the net sum of process rates and their stoichiometric
	 * coeficients
	 * 
	 * @return net sum of rate factors [g/um^3/h]
	 */
	public float getRate() {
		float r = 0;
		for (int i = 0; i < _addCounter; i++) {
			r += _coeficients[i] * _reactions[i].getRate();
		}
		return r;
	}
	/**
	 * @return the current Global rate for this net reaction
	 */
	public float getCurrentGlobalRate() {
		float r = 0;
		for (int i = 0; i < _addCounter; i++) {
			r += _coeficients[i] * _reactions[i].getGlobalReactionRate();
		}
		return r;
	}

	/**
	 * Computes the derivative of net sum of process rates and their
	 * stoichiometric coeficients in respect to variable
	 * 
	 * @param c
	 *            the chemical species to derivate rate
	 * @return net sum of rate factors [g/um^3/h]
	 */
	public float getRateDerivative(SoluteSpecies c) {
		float dr = 0;
		for (int i = 0; i < _addCounter; i++) {
			dr += _coeficients[i] * _reactions[i].getRateDerivative(c);
		}
		return dr;
	}

	/**
	 * Update the array with values of the rate and rate derivative
	 * 
	 * @param c the solute species for which this is being computed
	 * @param rDr [rate, rateDerivative]
	 */
	public void updateValuesForRateAndRateDerivative(SoluteSpecies c, float [] rDr) {
		rDr[0] = 0;
		rDr[1] = 0;
		for (int i = 0; i < _addCounter; i++) {
			_reactions[i].updateValuesForRateAndRateDerivative(c, _reactionRDr);
			rDr[0] += _coeficients[i] * _reactionRDr[0];
			rDr[1] += _coeficients[i] * _reactionRDr[1];
		}
	}
	
	/**
	 * @return the reactions as a List
	 */
	public List getReactionsAsList() {
		List l = Arrays.asList(_reactions);
		return l;
	}
}