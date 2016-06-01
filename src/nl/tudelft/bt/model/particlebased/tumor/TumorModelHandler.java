/*
 * Created on May 5, 2011
 */
package nl.tudelft.bt.model.particlebased.tumor;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.DilationFillHolesBoundaryLayer;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.util.ExtraMath;

import java.util.Collection;
import java.util.Iterator;
import java.util.Random;

/**
 * @author jxavier
 */
public abstract class TumorModelHandler extends GranuleModelHandler {

	/**
	 * Builds a new model handler for tumor simulations based on the granule
	 * model handler
	 */
	public TumorModelHandler() {
		super();
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		// uncomment the type of boundary layer wanted
		// _boundaryLayer = new FrontBoundaryLayer(h);
		_boundaryLayer = new DilationFillHolesBoundaryLayer(h);
		// create the boundary conditions
		MultigridVariable
				.setBoundaryConditions(new BiofilmBoundaryConditions());
	}

	/**
	 * Places cells at center of computational volume. Uses arrays to allow for
	 * any number of biomass species
	 * 
	 * 
	 * @param nCells
	 *            array of cells to place
	 * @param cancerBiomassSpecies
	 *            biomass species array (length must match that of nCells)
	 */
	protected void placeCellsAtCenter(int[] nCells,
			BiomassSpecies[] cancerBiomassSpecies) {
		int nspecies = nCells.length;
		if (nspecies != cancerBiomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + nCells.length
							+ " when it has " + cancerBiomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[nCells.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = nCells[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(Model.model().getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) Model.model().biomassContainer)
							.placeInoculumParticle(cancerBiomassSpecies[r]);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completely for homogenizing
		// the innoculum
		spreadByShovingCompletely();

	}

	/**
	 * Picks one cell of type "source" and switches by cell of type "target" 
	 * 
	 * @param source
	 * @param target
	 */
	protected void switchOneRandomnlyPickedCell(BiomassSpecies source,
			BiomassSpecies target) {
		// get the number of cells of type "source"
		int nSource = Model.model().biomassContainer
				.getNumberOfParticles(source);
		// pick one at random to switch
		int numberOfCellToSwitch = (int) Math
				.ceil(Model.model().getRandom() * (nSource-1));
		//System.out.println("switch cell #" + numberOfCellToSwitch + " of " + source + " by " + target);
		// get the biomass particle collection
		Collection<BiomassParticle> particleList = Model.model().biomassContainer
				.getBiomassAsBiomassParticleCollection();
		// iterate until finding the cell to switch
		int n = 0;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle cellToSwitch = (BiomassParticle) iter.next();
			if (cellToSwitch.isOfSpecies(source)) {
				n++;
				if (n == numberOfCellToSwitch) {
					System.out.println("switch cell #" + numberOfCellToSwitch + " of " + source + " by " + target);
					// get location of cell to switch
					float x = cellToSwitch.getCenterX();
					float y = cellToSwitch.getCenterY();
					float z = cellToSwitch.getCenterZ();
					// make a new cell of "target"
					Model.model().biomassContainer.placeBiomassAt(
							target, x, y, z);
					// remove source cell
					((BiomassParticleContainer)Model.model().biomassContainer).remove(cellToSwitch);
					break;

				}
			}
		}
	}

	/**
	 * Places cells throughout the entire computational volume. Uses arrays to
	 * allow for any number of biomass species
	 * 
	 * 
	 * @param nCells
	 *            array of cells to place
	 * @param macrophageBiomassSpecies
	 *            biomass species array (length must match that of nCells)
	 */
	protected void placeCellsRandomly(int[] nCells,
			BiomassSpecies[] macrophageBiomassSpecies) {
		if (nCells.length != macrophageBiomassSpecies.length) {
			new ModelRuntimeException(
					"number of cells must match number of biomass species");
		}
		// copy the n array
		int[] numberOfCellsToPlace = new int[nCells.length];
		int nspecies = macrophageBiomassSpecies.length;
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = nCells[i];
		//
		while (true) {
			// get the random value
			int r = (int) Math.floor(Model.model().getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float x = Model.model().getRandom()
							* Model.model().systemSize.x;
					float y = Model.model().getRandom()
							* Model.model().systemSize.y;
					float z = Model.model().getRandom()
							* Model.model().systemSize.z;
					// add to the bacteria list
					Model.model().biomassContainer.placeBiomassAt(
							macrophageBiomassSpecies[r], x, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	protected void placeCellsRandomlybetweenRadius(int[] nCells,
			BiomassSpecies[] macrophageBiomassSpecies, float rada, float radb) {
		if (nCells.length != macrophageBiomassSpecies.length) {
			new ModelRuntimeException(
					"number of cells must match number of biomass species");
		}
		// copy the n array
		int[] numberOfCellsToPlace = new int[nCells.length];
		int nspecies = macrophageBiomassSpecies.length;
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = nCells[i];
		//
		while (true) {
			// get the random value
			int r = (int) Math.floor(Model.model().getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					float radius = Model.model().getRandom() * (rada - radb)
							+ radb;
					((GranuleBiomassParticleContainer) Model.model().biomassContainer)
							.placeInoculumParticleAtRadius(
									macrophageBiomassSpecies[r], radius);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	protected void placeCellsAtRadius(int[] nCells,
			BiomassSpecies[] macrophageBiomassSpecies, float rad) {
		if (nCells.length != macrophageBiomassSpecies.length) {
			new ModelRuntimeException(
					"number of cells must match number of biomass species");
		}
		// copy the n array
		int[] numberOfCellsToPlace = new int[nCells.length];
		int nspecies = macrophageBiomassSpecies.length;
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = nCells[i];
		//
		while (true) {
			// get the random value
			int r = (int) Math.floor(Model.model().getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) Model.model().biomassContainer)
							.placeInoculumParticleAtRadius(
									macrophageBiomassSpecies[r], rad);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	protected void placeCellsAtRadiusandAngle(int[] nCells,
			BiomassSpecies[] macrophageBiomassSpecies, float rad, float a) {
		if (nCells.length != macrophageBiomassSpecies.length) {
			new ModelRuntimeException(
					"number of cells must match number of biomass species");
		}
		// copy the n array
		int[] numberOfCellsToPlace = new int[nCells.length];
		int nspecies = macrophageBiomassSpecies.length;
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = nCells[i];
		//
		while (true) {
			// get the random value
			int r = (int) Math.floor(Model.model().getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) Model.model().biomassContainer)
							.placeInoculumParticleAtRadiusandAngle(
									macrophageBiomassSpecies[r], rad, a);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

}