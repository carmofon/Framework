package nl.tudelft.bt.model.work.vanni;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.reaction.ProcessFactor;

public class MutatorBiomassSpecies extends BiomassSpecies {
	
	protected ParticulateSpecies[] arraySpecies;
	
	protected float _mutationRate;
	protected int _numberOfLevels;
	protected ProcessFactor _m0;
	
	/**
	 * @param name
	 * @param species
	 * @param fractionalCompositionInVolume
	 * @param mutationRate
	 * @param numberOfLevels the number of mutation levels
	 * @param wildTypeTrait
	 * @param mutantTrait
	 * @param neighborsTrait
	 * @param neighborhoodRadius
	 * @throws NonMatchingNumberException
	 */
	public MutatorBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume,
			float mutationRate, int numberOfLevels) throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		arraySpecies = species;
	    _numberOfLevels = numberOfLevels;
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
