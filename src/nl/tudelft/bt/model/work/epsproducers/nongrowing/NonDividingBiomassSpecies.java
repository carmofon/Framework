package nl.tudelft.bt.model.work.epsproducers.nongrowing;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements particles that do not divide for EPS producers paper
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) 18 October 2006
 */
public class NonDividingBiomassSpecies extends BiomassSpecies {
	public NonDividingBiomassSpecies(String name, ParticulateSpecies[] species,
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
		return new NonDividingBiomassParticle(this);
	}
}
