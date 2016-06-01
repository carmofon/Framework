package nl.tudelft.bt.model;

import java.io.IOException;
import java.io.Serializable;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Abstract class that defines all biomass spreading operations
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class BiomassContainer implements Serializable {
	/**
	 * Grow each biomass particle according to a time step
	 * 
	 * @param t
	 *            time step
	 */
	public abstract void grow(float t);

	/**
	 * Iterate through bacteria and perform division when necessary.
	 */
	public abstract void divideAndExcreteEPS();

	/**
	 * Relax completly by shoving
	 */
	public abstract void spread();

	public abstract void spreadByShoving();

	public abstract void performSpreadingByPulling();

	/**
	 * Perform a shoving iteration
	 * 
	 * @return number of shoved particles in this iteration
	 */
	public abstract int performSpreadingStep();

	/**
	 * Remove all cells located above the system limits
	 */
	public abstract void removeDetachedBiomass() throws ModelException;

	/**
	 * Iterate through all bacteria in the system removing "dead" cells.
	 */
	public abstract void removeDeadBiomass();

	/**
	 * Return the number of particles of a species
	 * 
	 * @param s
	 *            species
	 * @return number of particles presently in the container
	 */
	public abstract int getNumberOfParticles(BiomassSpecies s);

	/**
	 * Get the height of biofilm top
	 * 
	 * @return biofilm height
	 */
	public abstract float getTopBiofilmHeight();

	/**
	 * Get the total biomass in the fixed biofilm
	 * 
	 * @return total biofilm biomass
	 */
	public abstract float getTotalBiomass();

	/**
	 * computes a time constraint that prevents the biofilm to increase in size
	 * relative to a function of the maximum particle radius
	 * 
	 * @return the biofilm matrix quantity contraint
	 */
	public abstract TimeStepConstraint getGlobalSizeConstraint();

	/**
	 * Get the time constraint for the particle scale
	 * 
	 * @return time constraint for the particle scale
	 */
	public abstract TimeStepConstraint getParticleTimeConstraint();

	/**
	 * Add contribution of biomass particles to the biomass discrete data
	 * matrices
	 */
	public abstract void addContributionsToSpeciesDiscreteData();

	/**
	 * Add a biomass particle of species s at a defined position with a mass set
	 * by the biomass species default
	 * 
	 * @param s
	 *            species of biomass particle to add
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @return the biomass particle placed
	 */
	public abstract BiomassParticle placeBiomassAt(BiomassSpecies s, float x,
			float y, float z);

	/**
	 * Return the biomass as a collection of bacteria (spherical particles)
	 * 
	 * @return collection of biomass particles
	 */
	public abstract Collection getBiomassAsBiomassParticleCollection();

	/**
	 * @return the total produced biomass [10^-15 g] in this iteration
	 */
	public abstract float getTotalProducedBiomass();

	/**
	 * Write all bacteria (as spherical particles of varying color) information
	 * to file
	 * 
	 * @param fn
	 *            file name
	 */
	public void bacteriaToFile(String fn) throws IOException {
		java.io.File f = new java.io.File(fn);
		java.io.FileWriter fr = new java.io.FileWriter(f);
		for (Iterator iter = getBiomassAsBiomassParticleCollection().iterator(); iter
				.hasNext();) {
			// get the current bacterium
			BiomassParticle b = (BiomassParticle) iter.next();
			fr.write(b.toString());
		}
		fr.close();
	}

	/**
	 * Writes the detachment level set matrix to file before detachment is
	 * carried out
	 * 
	 * @param f
	 * @throws IOException
	 *             if write to file fails
	 */
	public abstract void writeDetachmentLevelSet(String f) throws IOException;

	/**
	 * writes the biofilm front as boolean matrix
	 * 
	 * @param f
	 * @throws IOException
	 */
	public abstract void writeBiofilmFront(String f) throws IOException;

	/**
	 * Returns the density at location c
	 * 
	 * @param c
	 * @return Biomass density at location c
	 */
	public abstract float getElementDensity(ContinuousCoordinate c);

	/**
	 * Returns the concentration of s particulate species at location c
	 * 
	 * @param c
	 * @return Biomass density at location c
	 */
	public abstract float getConcentrationInElement(ContinuousCoordinate c,
			ParticulateSpecies s);

	/**
	 * Returns the local concentration of EPS
	 * 
	 * @param c
	 * @return concentration of EPS at location c
	 */
	public abstract float getEpsInElement(ContinuousCoordinate c);

	/**
	 * Compute the biovolume from the mapping representation of the biomass
	 * 
	 * @return the total volume occupied by the biofilm (biovolume)
	 */
	public abstract float computeBiovolume();
	
	
	/**
	 * Turns off the removal of disconnected cells
	 */
	public abstract void turnSloughingOff();
}