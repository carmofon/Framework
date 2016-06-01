package nl.tudelft.bt.model.work.epsproducers.mutation;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class BiomassSpeciesForKeepingDistances extends BiomassSpecies {
	protected ParticulateSpecies _wildTypeTrait;

	protected ParticulateSpecies _mutantTrait;

	protected float _mutationTime;

	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param wildTypeTrait
	 * @param mutantTrait
	 * @param mutationTime
	 * @throws NonMatchingNumberException
	 */
	public BiomassSpeciesForKeepingDistances(String name,
			ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume,
			ParticulateSpecies wildTypeTrait, ParticulateSpecies mutantTrait,
			float mutationTime) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_wildTypeTrait = wildTypeTrait;
		_mutantTrait = mutantTrait;
		_mutationTime = mutationTime;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new MutatorBiomassParticleKeepDistances(this);
	}
}
