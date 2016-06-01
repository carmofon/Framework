package nl.tudelft.bt.model.work.epsproducers.distance;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class BiomassSpeciesTrackDistances extends BiomassSpecies {
	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @throws NonMatchingNumberException
	 */
	public BiomassSpeciesTrackDistances(String name,
			ParticulateSpecies[] species, float[] fractionalCompositionInVolume)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new ParticleTrackingDistances(this);
	}
}
