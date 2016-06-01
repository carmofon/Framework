
// CHANGE NAME OF PACKAGE TO MATCH THE LOCAITON IN THE FRAMEWORK
package nl.tudelft.bt.model.work.carlos;

import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassParticle;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassSpecies;

public class ImotileBiomassSpecies extends GranuleBiomassSpecies {

	public ImotileBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		// TODO Auto-generated constructor stub
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new ImotileBiomassParticle(this);
	}

}
