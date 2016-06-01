package nl.tudelft.bt.model.work.epsproducers.mutation;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class MutatorBiomassSpeciesForImpact extends BiomassSpecies {
	protected ParticulateSpecies _wildTypeTrait;

	protected ParticulateSpecies _mutantTrait;

	protected ParticulateSpecies _neighborsTrait;

	protected float _neighborhoodRadius;

	protected float _mutationTime;

	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param wildTypeTrait
	 * @param mutantTrait
	 * @param neighborsTrait
	 * @param mutationTime
	 * @param neighborhoodRadius
	 * @throws NonMatchingNumberException
	 */
	public MutatorBiomassSpeciesForImpact(String name,
			ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume,
			ParticulateSpecies wildTypeTrait, ParticulateSpecies mutantTrait,
			ParticulateSpecies neighborsTrait, float mutationTime,
			float neighborhoodRadius) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_wildTypeTrait = wildTypeTrait;
		_mutantTrait = mutantTrait;
		_neighborsTrait = neighborsTrait;
		_mutationTime = mutationTime;
		_neighborhoodRadius = neighborhoodRadius;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new MutatorBiomassParticle(this);
	}
}
