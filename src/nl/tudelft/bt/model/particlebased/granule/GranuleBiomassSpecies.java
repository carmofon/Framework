package nl.tudelft.bt.model.particlebased.granule;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements cyclic borders overall
 * 
 * 
 * @author xavierj
 *
 */
public class GranuleBiomassSpecies extends BiomassSpecies {

	public GranuleBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
	}

	
	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new GranuleBiomassParticle(this);
	}

}
