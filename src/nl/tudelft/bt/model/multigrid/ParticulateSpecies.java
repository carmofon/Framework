package nl.tudelft.bt.model.multigrid;

import java.awt.Color;
import java.util.List;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.apps.output.MassPerAreaTimeSeries;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.reaction.*;

/**
 * Biomass species multigrid variable
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ParticulateSpecies extends Species {
	private float _biomassDensity;

	private NetReaction _processes;

	private float _currentMaximumSpecificGrowthRate;

	private Color _color;
	
	private float _detachedBiomass;

	/**
	 * Creates a new fixed species
	 * 
	 * @param n
	 *            name of species
	 * @param dens
	 *            maximum biomass density [g/l = 10^-15g/um^3]
	 * @param c
	 *            the color to represent this species
	 * @throws MultigridSystemNotSetException
	 */
	public ParticulateSpecies(String n, float dens, Color c)
			throws MultigridSystemNotSetException {
		super();
		_name = n;
		_biomassDensity = dens;
		_color = c;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return _name;
	}

	/**
	 * @return the biomass density
	 */
	public float getDensity() {
		return _biomassDensity;
	}

	/**
	 * Get the net growth rate at a given position
	 * 
	 * @param c
	 *            position to get rate
	 * @return net gorwth rate at position [1/h]
	 */
	public float getSpecificRate(ContinuousCoordinate c) {
		setCurrentLocation(c);
		return _processes.getSpecificRate();
	}

	/**
	 * @return the mass growth rate from precomputed reaction rates [g/h]
	 */
	public float getMassRate() {
		return _processes.computeMassRateFromPreComputedReactionRates();
	}

	/**
	 * Compute the maximum possible growth rate (considering values of bulk
	 * conentration of chemicals involved)
	 * 
	 * sets the value of the maximum growth rate (considering maximum
	 * concentrations in system) [1/h]
	 */
	public void computeMaximumSpecificGrowthRate() {
		_currentMaximumSpecificGrowthRate = _processes.getSpecificRateMaximum();
	}

	/**
	 * Sets the processes and pre-computes the maximum rate for this species
	 * 
	 * @param stoichiometry
	 */
	public void setProcesses(NetReaction stoichiometry) {
		_processes = stoichiometry;
	}

	/**
	 * Set all values in the concentration matrix to 0
	 */
	public void resetDiscreteMatrix() {
		setValue(0.0f);
	}

	/**
	 * Refresh the data padding for boudary consitions
	 */
	public void refreshBoundaryConditions() {
		_boundaryConditions.refreshBoundaryConditions(_mg[_order - 1]);
	}

	/**
	 * Add contribution of a bacterium to the discrete grid of biomass
	 * concentration
	 * 
	 * @param c
	 *            position
	 * @param mass
	 *            to add
	 */
	public void addContributionToDiscreteData(ContinuousCoordinate c, float mass) {
		incrementValueAt(c, mass / _voxelVolume);
	}

	/**
	 * Return a formatted string with rate for this fixed species
	 * 
	 * @return
	 */
	private String finestRateToString() {
		_g = _order - 1;
		StringBuffer out = new StringBuffer();
		for (_k = 1; _k <= _l; _k++) {
			for (_i = _n; _i >= 1; _i--) {
				for (_j = 1; _j <= _m; _j++) {
					out.append(_processes.getSpecificRate());
					out.append(", ");
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	/**
	 * Save the rate for this bacteria species at the finest grid. If file is
	 * not written, process exits with error.
	 * 
	 * @param s
	 *            filename
	 */
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

	/**
	 * @return Returns the _currentMaximumSpecificGrowthRate.
	 */
	public float getCurrentMaximumSpecificGrowthRate() {
		return _currentMaximumSpecificGrowthRate;
	}

	/**
	 * Computes the total mass of this particle species in the computational
	 * volume
	 * 
	 * @return the total mass of this particle species per area [10-15g]
	 */
	public float getTotalMass() {
		//
		float concentrationSum = MultigridUtils.computeSum(_mg[_order - 1]);
		//
		return _voxelVolume * concentrationSum;
	}

	/**
	 * Returns a variable series
	 * 
	 * @return the variable series with the total mass in the biofilm
	 */
	public VariableSeries getTotalMassSeries() {
		// creates an intance of MassPerAreaTimeSeries overloading
		// getLastY()
		MassPerAreaTimeSeries m = new MassPerAreaTimeSeries(_name) {
			/*
			 * This method updates the Y array each time that last Y is called
			 * and X is longer than Y
			 */
			public float getLastY() {
				if (getYArray().getSize() < (getXArray().getSize() - 1)) {
					// Throw exception if time array is more than one entry
					// larger than the mass array
					throw new ModelRuntimeException(
							"Unsynchronized total mass series");
				} else if (getYArray().getSize() == (getXArray().getSize() - 1)) {
					getYArray().add(getTotalMass());
				}
				return super.getLastY();
			}
		};
		return m;
	}

	/**
	 * @return get the detached biomass accumulate in this iteration
	 */
	public float getDetachedBiomass() {
		return _detachedBiomass;
	}

	/**
	 * @return set the detached biomass accumulate in this iteration to 0
	 */
	public void resetDetachedBiomass() {
		_detachedBiomass = 0;
	}

	/**
	 * @return increment the detached biomass accumulate in this iteration by v
	 */
	public void addToDetachedBiomass(float v) {
		_detachedBiomass += v;
	}

	/**
	 * Get the reactions in which the species is involved by querying _processes
	 * attribute
	 * 
	 * @return the reactions involved as an List
	 */
	public List getReactionsInvlovedAsArrayList() {
		return _processes.getReactionsAsList();
	}

	/**
	 * @return Returns the _color representig this particulate species.
	 */
	public Color getColor() {
		return _color;
	}

	/**
	 * Read the discrete biomass concentrations from a file in a results
	 * directory
	 * 
	 * @param dir
	 *            the directory where results are located
	 * @param iteration
	 */
	public void readDiscreteConcentrationFromDir(String dir, int iteration) {
		// create the name of the file to read
		String f = dir + "/" + _name + (_order - 1) + ".txt";
		// read it
		resetDiscreteMatrix();
		// copy the 2D matrix read from file into the finest grid
		float[][] m = MultigridUtils.readSquareMatrixFromFile(f);
		for (int i = _n; i >= 1; i--) {
			for (int j = 1; j <= _m; j++) {
				_mg[_order - 1][i][j][1] = m[i - 1][j - 1];
			}
		}
		//
		refreshBoundaryConditions();
	}
}