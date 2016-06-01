package nl.tudelft.bt.model.multigrid;

import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements the class for solute species which react and diffuse
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SoluteSpecies extends Species {
	private BulkConcentration _bulkConcentration;

	protected MultigridVariable rhs;

	private float _diffusivity;

	private NetReaction _processes;

	protected float truncationError; // used in MG

	/**
	 * Create a chemical species
	 * 
	 * @param n
	 *            name of species
	 * @param d
	 *            diffusivity of species
	 * @throws MultigridSystemNotSetException
	 */
	public SoluteSpecies(String n, float d)
			throws MultigridSystemNotSetException {
		super();
		_name = n;
		_diffusivity = d;
		rhs = new MultigridVariable();
	}

	/**
	 * Set the value for bulk concentration of this chemical species
	 * 
	 * @param bc
	 *            the instance of bulk concentration
	 */
	public void setBulkConcentration(BulkConcentration bc) {
		_bulkConcentration = bc;
		_bulkConcentration.setSpecies(this);
	}

	/**
	 * @return the time step constraint
	 */
	public TimeStepConstraint getBulkSpeciesTimeStepConstraint() {
		return _bulkConcentration.getMaximumTimeStep();
	}

	/**
	 * Initialize the rate series by adding the first value
	 */
	public void initializeGlobalRateSeries() {
		_bulkConcentration.updateGlobalRateFromDiffusionReactionSeries();
		_bulkConcentration.updateGlobalRateFromMassBalance();
	}

	/**
	 * Update the value of the bulk concentration by evoking
	 * _bulkConcentration.computeBulkConcentration
	 * 
	 * @param t
	 *            the time step for this iteration
	 */
	public void updateBulkConcentrationAndRateSeries(float t) {
		// update the rate value
		_bulkConcentration.updateGlobalRateFromMassBalance();
		// Compute new concentration based on mass balance
		_bulkConcentration.computeBulkConcentration(t);
	}

	/**
	 * Creates a new rate time series which dfd
	 * 
	 * @return array of the global rate time series
	 */
	public VariableSeries getRateTimeSeries() {
		return _bulkConcentration.getRateTimeSeries();
	}

	/**
	 * @return the value of the bulk concentration for this solute
	 */
	public float getBulkConcentration() {
		return _bulkConcentration.getValue();
	}

	/**
	 * @return the series of bulk concentrations through the simulation
	 *         iterations.
	 */
	public VariableSeries getBulkConcentrationSeries() {
		return _bulkConcentration.getBulkConcentrationTimeSeries();
	}

	/**
	 * Set all values in the _mg matrices to the bulk
	 */
	void resetMultigridCopies() {
		for (int i = _order - 1; i > COARSEST; i--) {
			MultigridUtils.setValues(_mg[i], _bulkConcentration.getValue());
		}
	}

	/**
	 * Set value of coarsest grid to the bulk concnetration
	 */
	void setValueCoarsestToBulk() {
		MultigridUtils.setValues(_mg[COARSEST], _bulkConcentration.getValue());
	}

	/**
	 * @param stoichiometry
	 */
	public void setProcesses(NetReaction stoichiometry) {
		_processes = stoichiometry;
	}

	/**
	 * get the value of the chemical diffusvity
	 * 
	 * @return molecular diffussivity
	 */
	public float getDiffusivity() {
		return _diffusivity;
	}


	/**
	 * Compute the global biofilm convertion rate, if not yet computed
	 * 
	 * @return the integral of the rate for the whole computational volume
	 *         [MT-1]
	 */
	public float computeGlobalRateFromDiffusionReaction() {
		float r = 0;
		_g = _order - 1;
		for (_k = 1; _k <= _l; _k++) {
			for (_i = _n; _i >= 1; _i--) {
				for (_j = 1; _j <= _m; _j++) {
					r += getRate();
				}
			}
		}
		// multiplication by voxel volume is made in the end to save
		// computation time
		return r * _voxelVolume;
	}

	/**
	 * Compute the global biofilm convertion rate, if not yet computed
	 * 
	 * @return the integral of the rate for the whole computational volume
	 *         [MT-1/CV]
	 */
	public float computeGlobalRateFromMassBalance() {
		float r = _processes.getCurrentGlobalRate();
		// multiplication by voxel volume is made in the end to save
		// computation time
		return r;
	}

	/**
	 * Get the net rate at a the present position
	 * 
	 * @return net rate [g/um^3/h]
	 */
	public float getRate() {
		return _processes.getRate();
	}


	/**
	 * Update the array with values of the rate and rate derivative
	 * 
	 * @param rDr [rate, rateDerivative]
	 */
	public void updateValuesForRateAndRateDerivative(float [] rDr) {
		_processes.updateValuesForRateAndRateDerivative(this, rDr);
	}
	
	/**
	 * Return a formatted string with rate for this chemical
	 * 
	 * @return
	 */
	private String finestRateToString() {
		_g = _order - 1;
		StringBuffer out = new StringBuffer();
		for (_k = 1; _k <= _l; _k++) {
			for (_i = _n; _i >= 1; _i--) {
				for (_j = 1; _j <= _m; _j++) {
					out.append(getRate());
					out.append(", ");
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	public void finestRateToFile(String s) {
		java.io.File f = new java.io.File(s + _name + "Rate.txt");
		try {
			java.io.FileWriter fr = new java.io.FileWriter(f);
			fr.write(finestRateToString());
			fr.close();
		} catch (java.io.IOException e) {
			System.out.println(e);
			System.exit(-1);
		}
	}
}