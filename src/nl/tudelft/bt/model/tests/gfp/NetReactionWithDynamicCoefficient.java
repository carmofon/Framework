/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.tests.gfp;

import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.multigrid.Species;
import nl.tudelft.bt.model.reaction.NetReaction;
import nl.tudelft.bt.model.reaction.Reaction;

/**
 * Implements a container for processes and their sthoichiometric coeficients
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class NetReactionWithDynamicCoefficient extends NetReaction {
	private Reaction[] _reactions;

	private DynamicCoefficient[] _coeficients;

	private int _addCounter;
	/**
	 * Implements a reaction coefficient with dynamic behaviour
	 */
	private class DynamicCoefficient implements Serializable{
		float _value;

		float _t;

		float _duration;

		/**
		 * Create a new dynamic coefficient
		 * 
		 * @param v
		 */
		private DynamicCoefficient(float v) {
			_value = v;
			_t = 0;
			_duration = Float.POSITIVE_INFINITY;
		}

		/**
		 * Create a new dynamic coefficient
		 * 
		 * @param v
		 * @param t
		 * @param duration
		 */
		private DynamicCoefficient(float v, float t, float duration) {
			_value = v;
			_t = t;
			_duration = duration;
		}

		/**
		 * @return the current value of the coefficient
		 */
		private float getValue() {
			if ((Model.model().getTime() >= _t)
					& (Model.model().getTime() <= (_t + _duration))) {
				return _value;
			}
			return 0;
		}
	}

	/**
	 * creates a new instance o ReactionStoichiometry with space for n reactions
	 * 
	 * @param n
	 *            number of reactions
	 */
	public NetReactionWithDynamicCoefficient(int n) {
		super(n);
		_reactions = new Reaction[n];
		_coeficients = new DynamicCoefficient[n];
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
			_coeficients[_addCounter++] = new DynamicCoefficient(coef);
		} catch (ArrayIndexOutOfBoundsException e) {
			throw new ModelRuntimeException("Trying to add a "
					+ (_addCounter + 1) + " reaction to stoichiometry");
		}
	}

	/**
	 * add a new reaction to the stoichiometry container
	 * 
	 * @param r
	 *            reaction to add
	 * @param coef
	 *            stoichiometry coefficient
	 * @param t
	 *            time of begining of value coef
	 * @param duration
	 *            duration
	 */
	public void addReaction(Reaction r, float coef, float t, float duration) {
		try {
			_reactions[_addCounter] = r;
			_coeficients[_addCounter++] = new DynamicCoefficient(coef, t,
					duration);
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
			u += _coeficients[i].getValue()
					* _reactions[i].getSpecificRateFactor();
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
			if (_coeficients[i].getValue() > 0)
				u += _coeficients[i].getValue()
						* _reactions[i].getSpecificRateMaximum();
		}
		return u;
	}

	/**
	 * get the mass growth rate for a particle of composition c
	 * 
	 * @param c
	 *            composition of biomass particle
	 * @param speciesReacting
	 *            The species that is reacting
	 * 
	 * @return the mass rate for this stoichiometry [g/h]
	 */
	public float computeMassRateFromPreComputedReactionRates() {
		float r = 0;
		for (int i = 0; i < _addCounter; i++) {
			r += _coeficients[i].getValue()
					* _reactions[i].getPreComputedMassGrowthRate();
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
			r += _coeficients[i].getValue() * _reactions[i].getRate();
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
			dr += _coeficients[i].getValue()
					* _reactions[i].getRateDerivative(c);
		}
		return dr;
	}
}