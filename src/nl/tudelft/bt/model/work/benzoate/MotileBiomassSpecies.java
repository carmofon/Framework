package nl.tudelft.bt.model.work.benzoate;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) 12 October 2006
 */
public class MotileBiomassSpecies extends BiomassSpecies {
	// the motility diffusivity
	private float _diffusivity;

	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param motilityDiffusivity
	 * @throws NonMatchingNumberException
	 */
	public MotileBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, float motilityDiffusivity)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_diffusivity = motilityDiffusivity;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new MotileBiomassParticle(this);
	}

	public float getDiffusivity() {
		return _diffusivity;
	}
}
