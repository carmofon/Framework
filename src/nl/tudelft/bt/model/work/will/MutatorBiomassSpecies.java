package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class MutatorBiomassSpecies extends BiomassSpecies {
	protected ParticulateSpecies _wildTypeTrait;

	protected ParticulateSpecies _mutantTrait;

	protected float _mutationRate;

	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param wildTypeTrait
	 * @param mutantTrait
	 * @param neighborsTrait
	 * @param mutationRate
	 * @param neighborhoodRadius
	 * @throws NonMatchingNumberException
	 */
	public MutatorBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume,
			ParticulateSpecies wildTypeTrait, ParticulateSpecies mutantTrait,
			float mutationRate) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_wildTypeTrait = wildTypeTrait;
		_mutantTrait = mutantTrait;
		_mutationRate = mutationRate;
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
